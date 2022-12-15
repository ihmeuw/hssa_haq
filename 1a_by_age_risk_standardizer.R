##########################################################
# Date: 01/23/2019
# Description: Risk standardization
##########################################################
## DEFINE ROOT AND LIBRARIES
rm(list=ls())

library(data.table)
##########################################################
## PARSE ARGUMENTS
print("reading args")
arg <- commandArgs(trailingOnly = T)
if (length(arg)==0) {
  #toggle for arguments
  arg <- list(lid = 102,
              yid = 2010,
              prog_dir = "FILEPATH",
              vanilla_data_dir = "FILEPATH",
              data_dir = "FILEPATH",
              d_step = "step4",
              cycle = 6,
              lsid = 35,
              age = "young")
}
lid <- arg[[1]]
yid <- arg[[2]]
prog_dir <- arg[[3]]
vanilla_data_dir <- arg[[4]]
data_dir <- arg[[5]]
d_step <- arg[[6]]
cycle <- arg[[7]]
lsid <- arg[[8]]


if (length(arg) == 9){
  age <- arg[[9]]
}

years <- c(1990, 2000, 2010, 2015, 2019)
agid <- 12
mid <- 1
prefix <- "PAF"
scale_ceiling <- 90

popsdf <- fread(paste0(data_dir, "/popsdf_35.csv"))

##########################################################
## DEFINE FUNCTIONS
source("FILEPATH/utilities.R")

joinGlobal <- function(df, data_dir, lsid, mid, prefix) {
  # Load all years, and use average
  globdf <-
    rbindlist(lapply(as.list(years), function(yid, data_dir, prefix) {
      load(paste0(
        vanilla_data_dir,
        "/draws/inputs/",
        yid,
        "/",
        prefix,
        "_",
        mid,
        "_1.RData"
      ))
      return(inputdf)
    },
    data_dir, prefix))
  if (prefix == "PAF") {
    # Aggregate attributable burden over all years and join
    globdf <- globdf[, .(glob_paf = sum(paf * val) / sum(val), val = sum(val)), by = .(age_group_id, sex_id, measure_id, cause_id, draw)]
    globdf <- globdf[val == 0, glob_paf := 0]
    globdf$val <- NULL
    df <- merge(df,
                globdf,
                by = c("age_group_id", "sex_id", "measure_id", "cause_id", "draw"))
  }
  return(df)
}

rescale_make_age_weights <- function(cause_list, age){
  
  if (age == "young"){
    start <- 0
    end <- 15
  } else if ( age == "middle"){
    start <- 15
    end <- 65
  } else if ( age == "old"){
    start <- 65
    end <- 75
  }
  
  ageweightdf <- get_age_metadata(agid, cycle)
  cs_ageweightdf <- data.table()
  for (cid in cause_list$cause_id) {
    cause_weights <- data.table(ageweightdf)
    cause_weights$cause_id <- cid
    cause_weights <- cause_weights[age_group_years_start >= start & age_group_years_start < end]
    cause_weights <- cause_weights[, age_group_weight_value := age_group_weight_value / sum(age_group_weight_value)]
    cs_ageweightdf <- rbind(cs_ageweightdf,
                            cause_weights[, c("cause_id", "age_group_id", "age_group_weight_value"), with = FALSE])
  }
  return(cs_ageweightdf)
}

riskStandardize <- function(prefix, data_dir, yid, lid, mid, lsid, scale_ceiling, age) {
  
  # Specifying for by age analysis
  if (age == "young"){
    age_groups <- 2:7
    bin_age_group_id <- 39
  } else if ( age == "middle"){
    age_groups <- 8:17
    bin_age_group_id <- 200
  } else if ( age == "old"){
    age_groups <- 18:19
    bin_age_group_id <- 232
  }
  
  # Read in PAFs and CoDCorrect for each most-detailed location
  load(paste0(vanilla_data_dir, "/draws/inputs/", yid, "/", prefix, "_", mid, "_", lid, ".RData"))
  
  inputdf <- inputdf[age_group_id %in% age_groups]
  
  # Load location hierarchy data frame
  locsdf <- fread(paste0(data_dir, "/locsdf_", lsid, ".csv"))
  
  # Use deaths to aggregate PAF up to global
  inputdf <- joinGlobal(inputdf, data_dir, lsid, mid, prefix)
  
  if (prefix == "PAF") {
    # Apply scalars
    load(paste0(vanilla_data_dir, "/draws/inputs/PAF_", mid, "_scalar_", scale_ceiling, ".RData"))
    inputdf <- merge(inputdf,
                     scalardf,
                     by = c("age_group_id", "sex_id", "measure_id", "cause_id"))
    
    inputdf <- inputdf[, c("paf", "glob_paf") := .(paf * pafscalar, glob_paf * pafscalar)][, pafscalar := NULL]
    
    # Produce risk-standardized deaths
    inputdf <- inputdf[, rsval := val * (1 - paf) * (1 / (1 - glob_paf))]
    inputdf <- inputdf[is.na(rsval), rsval := 0]
    
    # For PAFs of one and for diarrhea, use observed
    inputdf <- inputdf[paf == 1 | cause_id == 302, rsval := val]
  }
  
  # Only keep Nolte/McKee causes and ages
  cause_list <- fread("FILEPATH/amenable_cause_list_GBD.csv")
  causesdf <- prepCauseHierarchy(data_dir)
  
  
  # Aggregate GBD causes and append amenable
  inputdf <- merge(inputdf,
                   causesdf,
                   by = "cause_id")
  
  inputdf <- rbindlist(lapply(as.list(c("root", paste0("L", seq(1, max(causesdf$level))))), hierarchyAgg, inputdf, "cause_id", c("rsval", "val")))
  
  # Aggregate to Nolte/McKee list, attach to main data frame
  inputdf <- merge(inputdf,
                   cause_list[, c("cause_id", "age_group_id_start", "age_group_id_end"), with = FALSE],
                   by = "cause_id")
  inputdf <- inputdf[age_group_id >= age_group_id_start & age_group_id <= age_group_id_end,]
  
  # put into rate space using the female and male specific populations for each cause
  popsdf <- fread(paste0(data_dir, "/popsdf_35.csv"))
  inputdf <- merge(inputdf,
                   popsdf,
                   by = c("location_id", "year_id", "age_group_id", "sex_id"))
  
  inputdf <- inputdf[, .(rsval = sum((rsval / population)),
                         val = sum((val / population))), 
                     by = .(location_id, year_id, age_group_id, sex_id, measure_id, cause_id, draw)]
  
  # Aggregate sexes
  inputdf <- inputdf[, bothsex := 3]
  inputdf <- hierarchyAgg("bothsex", inputdf, "sex_id", c("rsval", "val"))
  
  
  # Age-standardize 
  rescale_ageweightdf <- rescale_make_age_weights(cause_list, age)
  
  inputdf <- merge(inputdf,
                   rescale_ageweightdf,
                   all.x = T,
                   by = c("cause_id", "age_group_id"))

  if (nrow(inputdf[is.na(age_group_weight_value)]) > 0) stop("Problem with age weight merge")
  
  # Sum all risk-standardized deaths to see how the number of risk-standarized deaths varies by age
  
  inputdf <- inputdf[, age_group_id := bin_age_group_id]
  inputdf <- inputdf[, .(rsval = sum(rsval * age_group_weight_value),
                         val = sum(val * age_group_weight_value)), 
                     by = .(location_id, year_id, age_group_id, sex_id, cause_id, draw)]
  
  # Write data frame
  write.csv(inputdf[, c("location_id", "year_id", "age_group_id", "sex_id", "cause_id", "draw", "rsval", "val"), with = FALSE],
            file = paste0(data_dir, "/draws/standardized/", yid, "/", lid, "_", age, ".csv"),
            row.names = FALSE)
}

##########################################################
## RUN PROGRAM
# Risk standardize (should be changed from rsdeathsdf to rsdf, missed that in the update)
riskStandardize("PAF", data_dir, yid, lid, mid, lsid, scale_ceiling, "young")
riskStandardize("PAF", data_dir, yid, lid, mid, lsid, scale_ceiling, "middle")
riskStandardize("PAF", data_dir, yid, lid, mid, lsid, scale_ceiling, "old")
