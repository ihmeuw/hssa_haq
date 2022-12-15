##########################################################
# Date: 01/29/2019
# Description: Create by age index. This is a very memory intensive script.  
# Memory: 200 gb, 20 > threads
# Runtime: ~ 1.5 hours or more to run all three indices
##########################################################
## DEFINE ROOT AND LIBRARIES
rm(list=ls())
library(data.table)
library(feather)
library(broom)
library(parallel)

##########################################################
## List Args

yids <- c(1990, 2000, 2010, 2015, 2019)
data_dir <- "FILEPATH"
composite_data_dir <- "FILEPATH"
lsid <- 35
version <- "v2"
composite_version <- "v2"
use_composite_scale <- F

causes <- fread("FILEPATH/amenable_cause_list_GBD.csv")
mi_causes <- causes[measure == "MIR", cause_id]

##########################################################
## DEFINE FUNCTIONS
source("FILEPATH/utilities.R")

#Read in files for deaths and MI ratios

compileAmenable <- function(lids, yids, data_dir, age) {
  
  if (age == "young"){
    age <- "young"
  } else if (age == "middle"){
    age <- "middle"
  } else if (age == "old"){
    age <- "old"
  }
  
  # Load all specified risk-standardized mortality summaries
  # ncores <- round(detectCores()/2, 0)
  args <- expand.grid(location_id = lids, year_id = yids)
  amendf <- rbindlist(mcmapply(function(data_dir, yid, lid)  {
    if (file.exists(paste0(data_dir, yid, "/", lid, "_", age,".csv"))) {
      amenabledf<-fread(paste0(data_dir, yid, "/", lid, "_", age, ".csv"))
    }
  },
  yid = args$year_id,
  lid = args$location_id,
  MoreArgs = list(data_dir = data_dir),
  SIMPLIFY = FALSE,
  mc.cores = 8))
  return(amendf)
}



calcIndex <- function(base_lsid, yids, data_dir, age) {
  # Create directory
  dir.create(paste0(data_dir, '/'), recursive = T)
  
  # Get requisite datasets
  print("Getting locations")
  locsdf <- fread(paste0(data_dir, "/locsdf_", lsid, ".csv"))
  
  # Load age-standardized/risk standardized values, and MI ratios
  print("Compiling ASD Risk-Standardized Deaths")
  amenabledf_rsd <- compileAmenable(lids = locsdf$location_id, yids, data_dir = paste0(data_dir, "/"), age) 
  amenabledf_rsd <- amenabledf_rsd[!cause_id %in% mi_causes] # exclude non-relevant causes

  if ((length(unique(locsdf$location_id)) * length(unique(yids))* length(unique(amenabledf_rsd$cause_id)) * 1000) != nrow(amenabledf_rsd)) stop("Check for duplicates or missing values")
      # Check to see if correlation with SDI is negative for each cause
      load(paste0(data_dir, "/sdidf.RData"))
      
      rsd_sdi_df <- merge(sdidf[,.(location_id, year_id, sdi)], amenabledf_rsd, by = c("location_id", "year_id"))
      rsd_sdi_cor <- data.table()
      for (cause in unique(rsd_sdi_df$cause_id)) {
        
        cause_cor <- cor(rsd_sdi_df[cause_id == cause, sdi], rsd_sdi_df[cause_id == cause, rsval])
        colname <- paste0("correlation_", cause)
        row_value <- data.table(cause_id = cause, correlation = cause_cor, measure = "rsd")
        rsd_sdi_cor <- rbind(row_value, rsd_sdi_cor)  
        
      }
      
      if (age == "old"){
        # Remove congentinal heart anomalies. SDI correlation is in reverse direction
        rsd_sdi_cor <- rsd_sdi_cor[!cause_id == 643]
      } else {
      if (any(rsd_sdi_cor[, correlation] > 0)) stop("Check RSD Data, SDI correlation opposite of what's expexcted")
      }
      
  gc()
    
  # Compile MI Ratios
  print("Compiling ASD MI Ratios")
  amenabledf_MI_ratios <- compileAmenable(lids = locsdf$location_id, yids, paste0(data_dir, "FILEPATH"), age)
  amenabledf_MI_ratios <- amenabledf_MI_ratios[cause_id %in% mi_causes]
  
  if ((length(unique(locsdf$location_id)) * length(unique(yids))* length(unique(amenabledf_MI_ratios$cause_id)) * 1000) != nrow(amenabledf_MI_ratios)) stop("Check for duplicates or missing values")
  
      # Check to see if correlation with SDI is negative for each cause
      mi_sdi_df <- merge(sdidf[,.(location_id, year_id, sdi)], amenabledf_MI_ratios, by = c("location_id", "year_id"))
      mi_sdi_cor <- data.table()
      for (cause in unique(mi_sdi_df$cause_id)) {
        
        cause_cor <- cor(mi_sdi_df[cause_id == cause, sdi], mi_sdi_df[cause_id == cause, rsval])
        colname <- paste0("correlation_", cause)
        row_value <- data.table(cause_id = cause, correlation = cause_cor, measure = "mir")
        mi_sdi_cor <- rbind(row_value, mi_sdi_cor)  
        
      }
      
      if (any(mi_sdi_cor[,correlation] > 0)) stop("Check MIR Data, SDI correlation opposite of what's expexcted")
      
      # Combine correlations table's save for diagnostics
      all_cor <- rbind(mi_sdi_cor, rsd_sdi_cor)
      write.csv(all_cor, paste0(data_dir, "unscaled_sdi_correlations_", version, ".csv"), row.names = F)
      
      rm(mi_sdi_cor, rsd_sdi_cor)
      
  gc()
  
  
  amenabledf <- rbindlist(list(amenabledf_rsd, amenabledf_MI_ratios), use.names = T)
  rm(amenabledf_rsd, amenabledf_MI_ratios)
  
  # Check if rsval's greater than 1, likely there will be but still should know which causes are potentailly problematic 
  high_vals <- amenabledf[rsval > 1] # Note all of the causes with rsval > 1 were M/I
  write.csv(high_vals, paste0(data_dir, "high_risk_standardized_vals_",age,"_", version, ".csv"), row.names = F)
  
  if (nrow(amenabledf[!is.finite(rsval)]) > 0) stop("Non-finite values present")
  gc()
  
  if (use_composite_scale == T){
    
    minmaxdf <- fread(paste0(composite_data_dir, "all_age_both_sex_min_max_", composite_version, ".csv"))
  
    } else {
    
    # Calculate cause-specific index values 
    minmaxdf <- data.table(amenabledf)
    minmaxdf <- minmaxdf[location_id %in% locsdf$location_id[locsdf$level == 3]][, rsval := (rsval + 1e-6)] # use location 3 for setting min/max, log offset with 1 death per million
    minmaxdf <- minmaxdf[, .(min_rsval = quantile(rsval, .01), max_rsval = quantile(rsval, .99)), by = .(cause_id, draw)]
  }
  
  gc()

  # Scale amenabledf to min/max
  print("Combining MI Ratios/RSD Deaths then scaling to upper/lower bounds")
  scaleddf <- merge(amenabledf, minmaxdf,
                    by = c("cause_id", "draw"))
  scaleddf <- scaleddf[rsval < min_rsval, rsval := min_rsval][rsval > max_rsval, rsval := max_rsval]
  scaleddf <- scaleddf[, log_index_cause := (1 - ((log(rsval) - log(min_rsval)) / (log(max_rsval) - log(min_rsval)))) * 100]
  scaleddf <- scaleddf[log_index_cause < 0, log_index_cause := 0][log_index_cause > 100, log_index_cause := 100]
  
  gc()
  
  
  # Save both draws and summaries for unscaled data
  print("Saving unscaled draws")
  write_feather(amenabledf,
                path = paste0(data_dir, "amenabledf_", age,"_", version, ".feather"))
  amenabledf <- amenabledf[, .(rs_mval = mean(rsval), rs_lval = quantile(rsval, 0.025), rs_uval = quantile(rsval, 0.975)), by = .(location_id, year_id, age_group_id, sex_id, cause_id)]
  
  write_feather(amenabledf,
                path = paste0(data_dir, "amenabledf_", age,"_", version, ".feather"))
  
  rm(amenabledf)
  
  gc()
  
  # Calculate the composite index and label cause id 100
  aggdf <- scaleddf[, .(log_index_cause_mean = mean(log_index_cause),
                  log_index_cause_geom_mean = exp(mean(log(log_index_cause + 1)))), by = .(location_id, year_id, age_group_id, sex_id, draw)]
  aggdf[, cause_id := 100]
  indexdf <- rbind(scaleddf[, c("location_id", "year_id", "age_group_id", "sex_id", "cause_id", "draw", "log_index_cause"), with = FALSE],
              aggdf[, c("location_id", "year_id", "age_group_id", "sex_id", "cause_id", "draw", "log_index_cause_mean", "log_index_cause_geom_mean"), with = FALSE],
              fill = TRUE)
  
  # Enforce limits of 0-100 scale. Shouldn't be needed but here as a double-check.
  if (nrow(indexdf[log_index_cause < 0 | log_index_cause > 100]) > 0) stop(" Values outside of 0-100")
  
  # Save draws of index
  print("Saving scaled draws")
  write_feather(indexdf,
                path = paste0(data_dir, "indexdf_", age, "_", version, ".feather"))
  
  # Find mean, upper, lower of scaled draws
  # Have to do separately for non-pca weighted indices because they only exist for cause id 100
  indexdf_alt <- indexdf[cause_id == 100, 
                         .(index_mean = mean(log_index_cause_mean),
                           index_mean_lval = quantile(log_index_cause_mean, 0.025),
                           index_mean_uval = quantile(log_index_cause_mean, 0.975),
                           index_geom = mean(log_index_cause_geom_mean),
                           index_geom_lval = quantile(log_index_cause_geom_mean, 0.025),
                           index_geom_uval = quantile(log_index_cause_geom_mean, 0.975)),
                         by = .(location_id, year_id, age_group_id, sex_id, cause_id)] 
  

  # Can do it for all pca-weighted causes since values present for all causes including 100  
  indexdf <- indexdf[cause_id != 100, .(index = mean(log_index_cause),
                         index_lval = quantile(log_index_cause, 0.025),
                         index_uval = quantile(log_index_cause, 0.975)),
                     by = .(location_id, year_id, age_group_id, sex_id, cause_id)]
  
  # Merge non-pca indices onto data set
  indexdf <- rbindlist(list(indexdf, indexdf_alt), use.names = T, fill = T)
  
  
  print("Saving summmary of scaled draws")
  write_feather(indexdf,
                path = paste0(data_dir, "indexdf_",age,"_", version, ".feather"))
  }

##########################################################
## RUN PROGRAM
indexdf_young <- calcIndex(base_lsid, yids, data_dir, age = "young")
rm(indexdf_young)
indexdf_middle <- calcIndex(base_lsid, yids, data_dir, age = "middle" )
rm(indexdf_middle)
indexdf_old <- calcIndex(base_lsid, yids, data_dir, age = "old")
