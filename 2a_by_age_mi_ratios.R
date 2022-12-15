##########################################################
# Date: 04/23/2020
# Description: Make MI ratios for all causes
##########################################################
## Define Root
rm(list=ls())

library(data.table)
library(ggplot2)

#############################################################
## PARSE ARGUMENTS
print("reading args")
arg <- commandArgs(trailingOnly = T)
if (length(arg)==0) {
  #toggle for arguments
  arg <- list(lid = 102,
              yid = 2019,
              prog_dir = "FILEPATH",
              data_dir = "FILEPATH",
              cc_vers = 135,
              como_vers = 470,
              d_step = "step4",
              cycle = 6,
              lsid = 35,
              age = "young"
              
  )
}
lid <- arg[[1]]
yid <- arg[[2]]
prog_dir <- arg[[3]]
data_dir <- arg[[4]]
cc_vers <- arg[[5]]
como_vers <- arg[[6]]
d_step <- arg[[7]]
cycle <- arg[[8]]
lsid <- arg[[9]]

# if not running controller script set age to young
if (length(arg) == 10){
  age <- arg[[10]]
}

agid <- 12
causes <- c(297,302,322,328,338,339,340,341,366,380,429,432,435,441,468,484,487,
            492,493,494,498,508,527,529,531,534,545,587,589,643,708,849)

print(arg)
print("args read")

##########################################################
## DEFINE FUNCTIONS
source("FILEPATH/utilities.R")

##############################################################
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

makeMIratios <- function(data_dir, yid, lid, age) {
  
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
  
  print("pulling deaths")
  
  # Pull raw counts of deaths at draw level, for each cause
  deathsdf <- get_draws(rep("cause_id",length(causes)), gbd_id = causes, source = 'codcorrect', measure_id = 1, age_group_id = age_groups,
                      location_id = lid, year_id = yid, sex_id = c(1,2), gbd_round_id = cycle, decomp_step = d_step, metric_id = 1, version_id = cc_vers)
  
  print("deaths complete")
  
  # Make it long
  deathsdf <- melt(deathsdf, id.vars = c('location_id','year_id','sex_id','age_group_id','cause_id', 'measure_id', 'metric_id'),
                 value.name = 'deaths', variable.name = 'draw')
  deathsdf <- deathsdf[, draw := as.numeric(gsub("draw_", "", draw))]
  
  # Pull incidence, estimates do not exist for Hypertensive heart disorder (cause id = 498) 
  incidencedf <- get_draws(rep("cause_id",length(causes)), gbd_id = causes, source = 'como', measure_id = 6, age_group_id = age_groups,
                         location_id = lid, year_id = yid, sex_id = c(1,2), gbd_round_id = cycle, decomp_step = d_step, metric_id = 3, version_id = como_vers)
  # Make it long
  incidencedf <- melt(incidencedf, id.vars = c('location_id','year_id','sex_id','age_group_id','cause_id', 'measure_id', 'metric_id'), value.name = 'incidence', variable.name = 'draw')
  incidencedf <- incidencedf[, draw := as.numeric(gsub("draw_", "", draw))]
  
  
  #Merge dataframes. 380 neonatal, 498 hypertensive , 643 congenital heart anomalies will not be present. 
  #COMO assigns incidence for 380 and 643 as "birth". Codcorrect assigns deaths starting with age group id 2 = 0-6 days 
  inputdf <- merge(incidencedf[,list(location_id, year_id, sex_id, age_group_id, cause_id, draw, incidence)],
                   deathsdf[,list(location_id, year_id, sex_id, age_group_id, cause_id, draw, deaths)],
                   by = c("location_id", "year_id", "sex_id", "age_group_id", "cause_id", "draw"))
  
  #Put deaths in rate space by sex specific populations
  popsdf <- fread(paste0(data_dir, "/popsdf_35.csv"))
  inputdf <- merge(inputdf, popsdf, by = c("location_id", "age_group_id", "sex_id", "year_id"))
  inputdf <- inputdf[, .(deaths = (deaths / population), incidence), by = .(location_id, year_id, age_group_id, sex_id, draw, cause_id)]
  
  if (age == "young"){
  inputdf <- inputdf[!c(cause_id == 484 & age_group_id %in% c(2,3))]
  }
  
  # For age groups that there are relevant
  # Only keep Nolte/McKee causes and ages
  cause_list <- fread("FILEPATH/amenable_cause_list_GBD.csv")
  
  # Aggregate to Nolte/McKee list, attach to main data frame
  inputdf <- merge(inputdf,
                   cause_list[, c("cause_id", "age_group_id_start", "age_group_id_end"), with = FALSE],
                   by = "cause_id")
  inputdf <- inputdf[age_group_id >= age_group_id_start & age_group_id <= age_group_id_end,]

  print("age-standardizing")
  
  # Age-standardize
  rescale_ageweightdf <- rescale_make_age_weights(cause_list, age = age)
  inputdf <- merge(inputdf,
                   rescale_ageweightdf,
                   by = c("cause_id", "age_group_id"),
                   all.x = TRUE)
  
  if (nrow(inputdf[is.na(age_group_weight_value)]) > 0) stop("Problem with age weight merge")
  inputdf <- inputdf[, age_group_id := bin_age_group_id]
  inputdf <- inputdf[, sex_id := 3]
  inputdf <- inputdf[, .(deaths = sum(deaths * age_group_weight_value),
                           incidence = sum(incidence  * age_group_weight_value)),
                       by = .(location_id, year_id, age_group_id, sex_id, cause_id, draw)]
  
  ## Compute aggregated, all ages both sex MI ratio
  inputdf[,rsval := deaths/incidence]
  
  na_draws <- nrow(inputdf[which(!is.finite(rsval))])
  na_causes <- unique(inputdf[which(!is.finite(rsval)),cause_id])
  if (na_draws > 0) print(paste0("Number of NA draws in df is ", na_draws))
  if (na_draws > 0) print(paste0("Cause id ", na_causes, " has NA values"))
  
  
  inputdf[,val := NA] #Needed to keep shape of dfs uniform
  inputdf[, draw := as.numeric(draw)]
  inputdf<-inputdf[order(cause_id,draw)] 

  # Changed congenital heart anomalies and Hypertenisve Heart diseases to 0. Using RSD anyway
  inputdf[which(!is.finite(rsval)), rsval := 0] 
  
  # Return data frame
  write.csv(inputdf[, c("location_id", "year_id", "age_group_id", "sex_id", "cause_id", "draw", "rsval", "val"), with = FALSE],
            file = paste0(data_dir, "/draws/mi/", yid, "/", lid, "_", age, ".csv"),
            row.names = FALSE)
}

##########################################################
## RUN PROGRAM

makeMIratios(data_dir, yid, lid, age = "young")
makeMIratios(data_dir, yid, lid, age = "middle")
makeMIratios(data_dir, yid, lid, age = "old")
