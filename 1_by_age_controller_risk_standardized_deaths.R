##########################################################
# Date: 01/23/19
# Description: Controls mortality amenable to healthcare and avertable burden analyses
##########################################################
## DEFINE ROOT AND LIBRARIES
rm(list=ls())

##########################################################
## DEFINE FUNCTIONS
source("FILEPATH")

##########################################################
## ACTION TOGGLES
sweep.standardized <- F

##########################################################
## PARAMETERIZATION
vanilla_data_dir <- "FILEPATH"
prog_dir <- "FILEPATH"
vanilla_prog_dir <- "FILEPATH"
scale_ceiling <- 90
lsid <- 35
csid <- 3
cycle <- 6
d_step <- "step4"
cc_vers <- 135
mid <- 1
paf_vers <- "283"
yids <- c(1990, 2000, 2010, 2015, 2019)

##########################################################
## BUILD FILE SYSTEM
data_dir <- paste0("FILEPATH", cc_vers,"_paf_", paf_vers, "_amenable")
dir.create(data_dir)
dir.create(paste0(data_dir, "/draws"))

# Risk-standardized
if (sweep.standardized) {
  print("About to delete risk-standardized directory...")
  Sys.sleep(10)
  print("...")
  unlink(paste0(data_dir, "/draws/standardized"), recursive=TRUE)
  print("Finished.")
}
dir.create(paste0(data_dir, "/draws/standardized"))
for (yid in yids) {
  dir.create(paste0(data_dir, "/draws/standardized/", yid))
}

##########################################################
## SAVE COMMON DATA BASE FETCHES
if (!file.exists(paste0(data_dir, "locsdf_", lsid, ".csv"))) {
  locsdf <- get_location_metadata(location_set_id = lsid, gbd_round_id = cycle)
  write.csv(locsdf, file = paste0(data_dir, "/locsdf_", lsid, ".csv"), row.names = F) 
} else locsdf <- fread(paste0(data_dir, "locsdf_", lsid, ".csv"))

if (!file.exists(paste0(data_dir, "popsdf_", lsid, ".csv"))) {
  popsdf <- get_population(location_id = 'all', year_id = yids, sex_id = 'all', age_group_id = 'all', gbd_round_id = cycle, decomp_step = d_step, status = "best")
  write.csv(popsdf, file = paste0(data_dir, "/popsdf_", lsid, ".csv"), row.names = F) 
} else popsdf <-fread(paste0(data_dir, "popsdf_", lsid, ".csv"))

if (!file.exists(paste0(data_dir, "causesdf.csv"))) {
  causesdf <- get_cause_metadata(cause_set_id = csid, gbd_round_id = cycle)
  write.csv(causesdf, file = paste0(data_dir, "/causesdf_", lsid, ".csv"), row.names = F)
} else causesdf <- fread(paste0(data_dir, "causesdf_35.csv"))

if(!file.exists(paste0(data_dir, "/sdidf_", lsid, ".RData"))) {
  sdidf <- get_covariate_estimates(covariate_id = 881, year_id = yids, gbd_round_id = cycle, decomp_step = d_step)
  setnames(sdidf, "mean_value", "sdi")
  save(sdidf,
       file = paste0(data_dir, "/sdidf.RData"))
}

##########################################################
## SUBMISSION BY LOCATION SET
# Risk standardize by measure (only deaths)

# Risk-standardized mortality calculator, on Prod cluster

for (yid in yids) {
  for (lid in locsdf$location_id) {
    if (!file.exists(paste0(data_dir, "/draws/standardized/", yid, "/", lid, ".csv"))) {
      qsub_fair(
        jobname = paste0("by_age_rsd_", lid, "_year_", yid),
        shell = paste0(vanilla_prog_dir, "/r_shell.sh"),
        code = paste0(prog_dir, "/1a_by_age_risk_standardizer.R"),
        logloc = '/share/temp/sgeoutput/jyear',
        queue = "all.q",
        project = "health_sys",
        runtime = 40:00,
        threads = 2,
        args = list(
          lid,
          yid,
          prog_dir,
          vanilla_data_dir,
          data_dir,
          d_step,
          cycle,
          lsid),
        mem = 8
      )
      
    }
  }
}

# Check files
for (yid in yids) {
  for (lid in locsdf$location_id) {
    while (!file.exists(paste0(data_dir, "/draws/standardized/", yid, "/", lid, ".csv"))) {
      print(paste0("Waiting for risk-standardized deaths to be stored (location_id: ", lid, "; year_id: ", yid, ") -- ", Sys.time()))
      Sys.sleep(3)
    }
  }
}

##########################################################