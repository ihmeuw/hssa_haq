##########################################################
# Date: 04/23/2020
# Description: Qsubs age-standardized MI ratios 
##########################################################
## DEFINE ROOT AND LIBRARIES
rm(list=ls())

##########################################################
## DEFINE FUNCTIONS
source("FILEPATH/utilities.R")

##########################################################
## ACTION TOGGLES
sweep.mi <- F

##########################################################
## PARAMETERIZATION
prog_dir <- "FILEPATH"
vanilla_prog_dir <- "FILEPATH"
data_dir <- "FILEPATH"
lsid <- 35
yids <- c(1990, 2000, 2010, 2015, 2019)
csid <- 3
cycle <- 6
d_step <- "step4"
cc_vers <- 135
como_vers <- 470

##########################################################
##########################################################
## BUILD FILE SYSTEM
if (sweep.mi) {
  print("About to delete mi directory...")
  Sys.sleep(10)
  print("...")
  unlink(paste0(data_dir, "/draws/mi"), recursive= T)
  print("Finished.")
}

if(!file.exists(data_dir)){
  dir.create(data_dir)
  dir.create(paste0(data_dir, "draws"))
}

# Inputs
if(!file.exists(paste0(data_dir, "draws/mi"))){
  dir.create(paste0(data_dir, "draws/mi"))
  for (yid in yids) {
    dir.create(paste0(data_dir, "draws/mi/", yid))
    
  }
}
##########################################################
#########################################################
## SAVE COMMON DATA BASE FETCHES
if (!file.exists(paste0(data_dir, "locsdf_", lsid, ".csv"))) {
  locsdf <- get_location_metadata(location_set_id = lsid, gbd_round_id = cycle)
  write.csv(locsdf, file = paste0(data_dir, "locsdf_", lsid, ".csv"), row.names = F) 
} else locsdf <- fread(paste0(data_dir, "locsdf_", lsid, ".csv"))

if (!file.exists(paste0(data_dir, "popsdf_", lsid, ".csv"))) {
  popsdf <- get_population(location_id = 'all', year_id = yids, sex_id = 'all', age_group_id = 'all', gbd_round_id = cycle, decomp_step = d_step, status = "best")
  write.csv(popsdf, file = paste0(data_dir, "popsdf_", lsid, ".csv"), row.names = F) 
} else popsdf <-fread(paste0(data_dir, "popsdf_", lsid, ".csv"))

if (!file.exists(paste0(data_dir, "causesdf_", lsid, ".csv"))) {
  causesdf <- get_cause_metadata(cause_set_id = csid, gbd_round_id = cycle)
  write.csv(causesdf, file = paste0(data_dir, "causesdf_", lsid, ".csv"), row.names = F)
} else causesdf <- fread(paste0(data_dir, "causesdf_35.csv"))

if(!file.exists(paste0(data_dir, "/sdidf.RData"))) {
  sdidf <- get_covariate_estimates(covariate_id = 881, year_id = yids, gbd_round_id = cycle, decomp_step = d_step)
  setnames(sdidf, "mean_value", "sdi")
  save(sdidf,
       file = paste0(data_dir, "/sdidf.RData"))
}

##########################################################
# Submit Jobs
# MI ratio qsubmitter 
for (yid in yids) {
  for (lid in unique(locsdf$location_id)) {
    if (!file.exists(paste0(data_dir, "draws/mi/", yid, "/", lid, ".csv"))) {
      qsub_fair(
        jobname = paste0("by_age_mi_", lid , "_", yid),
        shell = paste0(vanilla_prog_dir, "/r_shell.sh"),
        code = paste0(prog_dir, "/2a_by_age_mi_ratios.R"),
        logloc = 'FILEPATH',
        queue = "all.q",
        project = "health_sys",
        runtime = 40,
        threads = 2,
        args = list(lid, yid, prog_dir, data_dir, cc_vers, como_vers, d_step, cycle, lsid),
        mem = 8)
      
    }
  }
}

# Check files
for (yid in yids) {
  for (lid in locsdf$location_id) {
    while (!file.exists(paste0(data_dir, "/draws/mi/", yid, "/", lid, ".csv"))) {
      print(paste0("Waiting for MIRs to be stored (location_id: ", lid, "; year_id: ", yid, ") -- ", Sys.time()))
      Sys.sleep(30)
    }
  }
}

##########################################################