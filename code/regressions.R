
##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
neededPackages = c("plyr", "dplyr", "here", #"tibble", "data.table",
                   "foreign", # "readxl",
                   # "raster", "rgdal",  "sp", "sf", # "spdep",
                   "DataCombine",
                   "knitr", "kableExtra",
                   "fixest" #,"msm", "car",  "sandwich", "lmtest", "boot", "multcomp",
                   #"ggplot2", "leaflet", "htmltools"
                   )
# "pglm", "multiwayvcov", "clusterSEs", "alpaca", "clubSandwich",

# Install them in their project-specific versions
renv::restore(packages = neededPackages)

# Load them
lapply(neededPackages, library, character.only = TRUE)

# /!\/!\ IF renv::restore(neededPackages) FAILS TO INSTALL SOME PACKAGES /!\/!\ 

# For instance sf could cause trouble https://github.com/r-spatial/sf/issues/921 
# or magick, as a dependency of raster and rgdal. 

# FOLLOW THESE STEPS:
# 1. Remove these package names from neededPackages above, and rerun renv::restore(packages = neededPackages)
# 2. Write them in troublePackages below, uncomment, and run the following code chunk: 

# # /!\ THIS BREAKS THE PROJECT REPRODUCIBILITY GUARANTY /!\
# troublePackages <- c() 
# # Attempt to load packages from user's default libraries.
# lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

# 3. If the troubling packages could not be loaded ("there is no package called ...") 
#   you should try to install them, preferably in their versions stated in the renv.lock file. 
#   see in particular https://rstudio.github.io/renv/articles/renv.html 


# # # /!\ THIS BREAKS THE PROJECT REPRODUCIBILITY GUARANTY /!\
# troublePackages <- c("leaflet", "leaflet.providers", "png")
# # Attempt to load packages from user's default libraries.
# lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing

### NEW FOLDERS USED IN THIS SCRIPT 
dir.create("temp_data/reg_results")

### SET NUMBER OF THREADS USED BY {FIXEST} TO ONE (TO REDUCE R SESSION CRASH)
getFixest_nthreads()





### READ ALL POSSIBLE DATASETS HERE 
# but not all together because of memory issues. 
d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_gaez_long_final.Rdata"))

# d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_gaez_long_final.Rdata"))
# 
# d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_gaez_long_final.Rdata"))

prices <- readRDS(here("temp_data", "prepared_prices.Rdata"))

dataset = "glass"
start_year = 1983
end_year = 2015
crop_j = "Oilpalm"
j_soy = "Soybean_oil"
price_k = "Soybean_oil"
standardized_si = TRUE
price_lag = 1
SjPj = TRUE
SkPk = TRUE
fe = "grid_id + country_year"
distribution = "quasipoisson"
output_full = FALSE

make_base_reg <- function(dataset, 
                          start_year = 2002, 
                          end_year = 2015, 
                          crop_j = "Oilpalm", # in GAEZ spelling
                          j_soy = "Soybean_oil", # in case crop_j is Soybean, which price should be used? One of "Soybeans", "Soybean_oil", "Soybean_meal".
                          price_k = "Soybean_oil", # in prices spelling
                          standardized_si = TRUE,
                          price_lag = 1, 
                          SjPj = TRUE,
                          SkPk = TRUE,
                          fe = "grid_id + country_year", 
                          distribution = "quasipoisson", 
                          output_full = FALSE){

   

  
  
  ### SPECIFICATIONS  
  
  ## Outcome variable
  if(dataset=="glass"){
    # d <- glass
    outcome_variable <- "first_loss"} #   "sbqt_direct_lu"      "sbqt_mode_lu" 
  if(dataset=="fl8320"){
    # d <- fl8320
    outcome_variable <- "firstloss_glassgfc"}
  if(dataset=="phtfl"){
    # d <- phtfl
    outcome_variable <- "phtf_loss"}
  
  # We need to create the variables we will need, based on:
  # crop_j, price_k, whether the price is lagged, whether the suitability is standardized, and the controls we want 
  
  ## Identify the price of j based on crop_j (GAEZ spelling)
  if(crop_j == "Pasture"){price_j <- "Beef"} 
  if(crop_j == "Oilpalm"){price_j <- "Palm_oil"} 
  if(crop_j == "Soybean"){price_j <- j_soy} 
  if(crop_j == "Cocoa"){price_j <- "Cocoa"} 
  if(crop_j == "Coffee"){price_j <- "Coffee"} 
  
  # and lag it
  if(price_lag != 0){price_j <- paste0(price_j,"_lag",price_lag)}
  
  ## Do the revert for k: from price of k to the GAEZ crop 
  if(grepl(pattern = "Banana", x = price_k, ignore.case = TRUE)){crop_k <- "Banana"} 
  if(grepl(pattern = "Barley", x = price_k, ignore.case = TRUE)){crop_k <- "Barley"} 
  if(grepl(pattern = "Beef", x = price_k, ignore.case = TRUE)){crop_k <- "Pasture"} # currently does not exist
  if(grepl(pattern = "Orange", x = price_k, ignore.case = TRUE)){crop_k <- "Citrus"} 
  if(grepl(pattern = "Cocoa", x = price_k, ignore.case = TRUE)){crop_k <- "Cocoa"} 
  if(grepl(pattern = "Coconut", x = price_k, ignore.case = TRUE)){crop_k <- "Coconut"} 
  if(grepl(pattern = "Coffee", x = price_k, ignore.case = TRUE)){crop_k <- "Coffee"} 
  if(grepl(pattern = "Groundnuts", x = price_k, ignore.case = TRUE)){crop_k <- "Groundnut"} 
  if(grepl(pattern = "Maize", x = price_k, ignore.case = TRUE)){crop_k <- "Maize"} 
  if(grepl(pattern = "Oat", x = price_k, ignore.case = TRUE)){crop_k <- "Oat"} # aggregate different crops there ? like rye... 
  if(grepl(pattern = "Olive", x = price_k, ignore.case = TRUE)){crop_k <- "Olive"} 
  if(grepl(pattern = "Palm", x = price_k, ignore.case = TRUE)){crop_k <- "Oilpalm"} 
  if(grepl(pattern = "Rapeseed", x = price_k, ignore.case = TRUE)){crop_k <- "Rapeseed"} 
  if(grepl(pattern = "Rice", x = price_k, ignore.case = TRUE)){crop_k <- "Rice"} # currently does not exist
  if(grepl(pattern = "Sorghum", x = price_k, ignore.case = TRUE)){crop_k <- "Sorghum"} 
  if(grepl(pattern = "Soy", x = price_k, ignore.case = TRUE)){crop_k <- "Soybean"} 
  if(grepl(pattern = "Sugar", x = price_k, ignore.case = TRUE)){crop_k <- "sugar_crops"} 
  if(grepl(pattern = "Sunflower", x = price_k, ignore.case = TRUE)){crop_k <- "Sunflower"} 
  if(grepl(pattern = "Tea", x = price_k, ignore.case = TRUE)){crop_k <- "Tea"} 
  if(grepl(pattern = "Tobacco", x = price_k, ignore.case = TRUE)){crop_k <- "Tobacco"} 
  if(grepl(pattern = "Wheat", x = price_k, ignore.case = TRUE)){crop_k <- "Wheat"} # aggregate different wheat crops there ? 
  
  # Construct the price of k
  if(price_lag != 0){price_k <- paste0(price_k,"_lag",price_lag)}
  
  
  # don't condition that, as we won't use non logged prices a priori. 
  price_j <- paste0("ln_", price_j)
  price_k <- paste0("ln_", price_k)
  
  # Identify the suitability index needed - don't condition that either for now, as the non stded indexes are not in the main dataset now (for memory issues) we'll see later how to conduct the rob check
  # if(standardized_si){
    suitability_j <- paste0(crop_j, "_std") 
    suitability_k <- paste0(crop_k, "_std") 
  # }else{
    # suitability_j <- crop_j
    # suitability_k <- price_j
  # }
  
  # Keep only in data the useful variables 
  d <- dplyr::select(d, all_of(c("grid_id", "year", "country_name", "country_year",
                               outcome_variable, suitability_j, suitability_k))) 
    
    
  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", price_j, price_k)], by = "year")
  
  # Create the variables
  d <- mutate(d, 
              SjPk = !!as.symbol(suitability_j)*!!as.symbol(price_k), 
              SkPk = !!as.symbol(suitability_k)*!!as.symbol(price_k), 
              SjPj = !!as.symbol(suitability_j)*!!as.symbol(price_j))
  
  
  ### KEEP OBSERVATIONS THAT: 
  
  # - are in study period 
  d <- dplyr::filter(d, year >= start_year)
  d <- dplyr::filter(d, year <= end_year)
  
  used_vars <- names(d)
  
  # - have no NA nor INF on any of the variables used (otherwise they get removed by {fixest})
  # for instance, there are some NAs in the suitability index (places in water that we kept while processing other variables...) 
  usable <- lapply(used_vars, FUN = function(var){is.finite(d[,var]) | is.character(d[,var])})
  # used_vars[!(used_vars%in%names(d))]
  names(usable) <- used_vars            
  usable <- bind_cols(usable)
  filter_vec <- base::rowSums(usable)
  filter_vec <- filter_vec == length(used_vars)
  d_nona <- d[filter_vec, c(used_vars)]
  if(anyNA(d_nona)){stop()}
  rm(filter_vec, usable)
  
  # is.na(d$Oilpalm_std) %>% sum()
  
  
  # - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove
  d_clean <- d_nona[-obs2remove(fml = as.formula(paste0(outcome_variable, " ~ ", fe)),
                                d_nona, 
                                family = "poisson"),]
  
  
  ### REGRESSIONS
  controls <- c()
  if(SkPk){controls <- c(controls, "SkPk")}
  if(SjPj){controls <- c(controls, "SjPj")}
  
  
  # Model specification
  if(length(controls) > 0){
    fe_model <- as.formula(paste0(outcome_variable,
                                  " ~ SjPk",
                                  " + ",
                                  paste0(controls, collapse = "+"),
                                  " | ",
                                  fe))
  }else{
    fe_model <- as.formula(paste0(outcome_variable,
                                  " ~ SjPk",
                                  " | ",
                                  fe))
  }

  
  if(distribution != "negbin"){ # i.e. if it's poisson or quasipoisson or gaussian
    reg_res <- fixest::feglm(fe_model,
                             data = d_clean, 
                             family = distribution, 
                             #glm.iter = 100,
                             #fixef.iter = 100000,
                             notes = TRUE)
  }
  
  if(output_full){
    toreturn <- list(reg_res, d_clean, d)
  }else{
    toreturn <- list(reg_res, d_clean)
  }
  
  rm(d, d_nona, d_clean, reg_res)
  return(toreturn)
  
}


K_commodities <- c("ln_Soybean_oil", "ln_Olive_oil", "ln_Rapeseed_oil", "ln_Sunflower_oil", "ln_Coconut_oil", "ln_Sugar", "ln_Crude_oil")
for(commodity in K_commodities){
  
}




