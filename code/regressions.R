
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
glass <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_gaez_long_country_nf.Rdata"))

fl8320 <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_gaez_long_country_nf.Rdata"))

phtfl <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_gaez_long_country_nf.Rdata"))

prices <- readRDS(here("temp_data", "prepared_prices.Rdata"))

make_base_reg <- function(dataset, 
                          start_year = 2002, 
                          end_year = 2015, 
                          commo_j = "Oilpalm", # in GAEZ spelling
                          commo_k = c("ln_Soybean_oil", "ln_Olive_oil", "ln_Rapeseed_oil", "ln_Sunflower_oil", "ln_Coconut_oil", "ln_Sugar", "ln_Crude_oil"),
                          price_lag = 1, 
                          fe = "grid_id + country_year"){
  
  
  
  if(dataset=="glass"){
    d <- glass
    outcome_variable <- "first_loss"}
  if(dataset=="fl8320"){
    d <- fl8320
    outcome_variable <- "firstloss_glassgfc"}
  if(dataset=="phtfl"){
    d <- phtfl
    outcome_variable <- "phtf_loss"}
  
  
  ### SPECIFICATIONS  
  
  
  ### KEEP OBSERVATIONS THAT: 
  
  # - are in study period 
  d <- dplyr::filter(d, year >= start_year)
  d <- dplyr::filter(d, year <= end_year)
  
  used_vars <- c(outcome_variable, regressors,
                 "grid_id",  "year", "lat", "lon", 
                 "country", "country_year")
  
  # - have no NA nor INF on any of the variables used (otherwise they get removed by {fixest})
  usable <- lapply(used_vars, FUN = function(var){is.finite(d[,var]) | is.character(d[,var])})
  # used_vars[!(used_vars%in%names(d))]
  names(usable) <- used_vars            
  usable <- bind_cols(usable)
  filter_vec <- base::rowSums(usable)
  filter_vec <- filter_vec == length(used_vars)
  d_nona <- d[filter_vec, c(used_vars)]
  if(anyNA(d_nona)){stop()}
  rm(filter_vec, usable)
}








