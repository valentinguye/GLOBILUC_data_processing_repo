

##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
neededPackages = c("Matrix",
                   "plyr", "dplyr", "here", #"tibble", "data.table",
                   "foreign", # "readxl",
                   "raster", "rgdal",  "sp", "sf", # "spdep",
                   "DataCombine",
                   "knitr", "kableExtra",
                   "fixest", "boot",#,"msm", "car",  "sandwich", "lmtest",  "multcomp",
                   "ggplot2", "dotwhisker", "viridis", "hrbrthemes", "rasterVis", #"tmap",# "leaflet", "htmltools"
                   "foreach", "parallel"
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
dir.create(here("temp_data","reg_results", "alpha"))
dir.create(here("temp_data","reg_results", "beta"))
dir.create(here("temp_data","reg_results", "rfs"))

### SET NUMBER OF THREADS USED BY {FIXEST} TO ONE (TO REDUCE R SESSION CRASH)
getFixest_nthreads()

setFixest_dict(c(grid_id = "grid cell",
                 country_name = "country",
                 country_year = "country*year", 
                 driven_loss = "Forest loss",
                 first_loss = "First forest loss",
                 Soybean_meal = "Soybean meal",
                 Soybean_oil = "Soybean oil",
                 Olive_oil = "Olive oil",
                 Rapeseed_oil = "Rapeseed oil",
                 Sunflower_oil = "Sunflower oil",
                 Coconut_oil = "Coconut oil", 
                 Fodder_X_Beef = "Beef", Fodder_X_Soybean = "Soybean", Fodder_X_Maize = "Maize", Fodder_X_Sugar = "Sugar", Fodder_X_Wheat = "Wheat", Fodder_X_Barley = "Barley", Fodder_X_Oat = "Oat", Fodder_X_Chicken = "Chicken", Fodder_X_Pork = "Pork", Fodder_X_Palm_oil = "Palm oil",Fodder_X_Rapeseed_oil = "Rapeseed oil",Fodder_X_Sunflower_oil = "Sunflower oil", 
                 Oilpalm_X_Beef = "Beef", Oilpalm_X_Soybean = "Soybean", Oilpalm_X_Maize = "Maize", Oilpalm_X_Wheat = "Wheat", Oilpalm_X_Barley = "Barley", Oilpalm_X_Oat = "Oat", Oilpalm_X_Chicken = "Chicken", Oilpalm_X_Pork = "Pork", Oilpalm_X_Palm_oil = "Palm oil", Oilpalm_X_Sugar = "Sugar",Oilpalm_X_Rapeseed_oil = "Rapeseed oil",Oilpalm_X_Sunflower_oil = "Sunflower oil",
                 Soybean_X_Beef = "Beef", Soybean_X_Soybean = "Soybean", Soybean_X_Maize = "Maize", Soybean_X_Wheat = "Wheat", Soybean_X_Barley = "Barley", Soybean_X_Oat = "Oat", Soybean_X_Chicken = "Chicken", Soybean_X_Pork = "Pork", Soybean_X_Palm_oil = "Palm oil", Soybean_X_Sugar = "Sugar",Soybean_X_Rapeseed_oil = "Rapeseed oil",Soybean_X_Sunflower_oil = "Sunflower oil",
                 Cocoa_X_Cocoa = "Cocoa", Cocoa_X_Coffee = "Coffee", Cocoa_X_Beef = "Beef", Cocoa_X_Soybean = "Soybean", Cocoa_X_Maize = "Maize", Cocoa_X_Wheat = "Wheat", Cocoa_X_Barley = "Barley", Cocoa_X_Oat = "Oat", Cocoa_X_Chicken = "Chicken", Cocoa_X_Pork = "Pork", Cocoa_X_Palm_oil = "Palm oil", Cocoa_X_Sugar = "Sugar",Cocoa_X_Rapeseed_oil = "Rapeseed oil",Cocoa_X_Sunflower_oil = "Sunflower oil",
                 Coffee_X_Coffee = "Coffee", Coffee_X_Cocoa = "Cocoa", Coffee_X_Tea = "Tea", Coffee_X_Tobacco = "Tobacco",  Coffee_X_Beef = "Beef", Coffee_X_Soybean = "Soybean", Coffee_X_Maize = "Maize", Coffee_X_Wheat = "Wheat", Coffee_X_Barley = "Barley", Coffee_X_Oat = "Oat", Coffee_X_Chicken = "Chicken", Coffee_X_Pork = "Pork", Coffee_X_Palm_oil = "Palm oil", Coffee_X_Sugar = "Sugar", Coffee_X_Rapeseed_oil = "Rapeseed oil", Coffee_X_Sunflower_oil = "Sunflower oil",
                 Rubber_X_Rubber = "Rubber", Rubber_X_Beef = "Beef", Rubber_X_Soybean = "Soybean", Rubber_X_Maize = "Maize", Rubber_X_Wheat = "Wheat", Rubber_X_Barley = "Barley", Rubber_X_Oat = "Oat", Rubber_X_Chicken = "Chicken", Rubber_X_Pork = "Pork", Rubber_X_Palm_oil = "Palm oil", Rubber_X_Sugar = "Sugar", Rubber_X_Rapeseed_oil = "Rapeseed oil", Rubber_X_Sunflower_oil = "Sunflower oil"
))


### GLOBAL CRS USED ### 
mercator_world_crs <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "


# this is mapmat for eaear 
eaear_mapmat_data <- c(
  "Banana","eaear_Banana",
  "Crude_oil", "eaear_Biomass",
  #"Barley", "eaear_Barley",
  "Cereals", "eaear_Cereals",
  "Orange", "eaear_Citrus",
  #"Cocoa", "eaear_Cocoa",
  "Cocoa_Coffee", "eaear_Cocoa_Coffee",
  "Coconut_oil", "eaear_Coconut",
  #"Coffee", "eaear_Coffee",
  "Cotton", "eaear_Cotton",
  #"Groundnuts", "eaear_Groundnut",
  "Beef", "eaear_Fodder", 
  "Maize", "eaear_Maizegrain",
  "Oilfeed_crops", "eaear_Oilfeed_crops",
  # "Oat", "eaear_Oat",
  # "Olive_oil", "eaear_Olive",  
  "Palm_oil", "eaear_Oilpalm",
  #"Rapeseed_oil", "eaear_Rapeseed",
  "Rice", "eaear_Rice",
  "Rubber", "eaear_Rubber",
  #"Sorghum", "eaear_Sorghum2", 
  "Soy_index", "eaear_Soy_compo",
  "Sugar", "eaear_Sugarcane", 
  #"Sunflower_oil", "eaear_Sunflower",
  "Tea", "eaear_Tea",
  "Tobacco", "eaear_Tobacco" 
  #"Wheat", "eaear_Wheat"
)

eaear_mapmat <- matrix(data = eaear_mapmat_data, 
                       nrow = length(eaear_mapmat_data)/2,
                       ncol = 2, 
                       byrow = TRUE)

colnames(eaear_mapmat) <- c("Prices", "Crops")

agri_crops <- eaear_mapmat[,"Crops"][!(eaear_mapmat[,"Crops"] %in% c("eaear_Oilpalm", "eaear_Rubber"))]
plantation_crops <- c("eaear_Oilpalm", "eaear_Rubber")

limited_crops1 <- c("eaear_Maizegrain", "eaear_Fodder", "eaear_Cereals", "eaear_Soy_compo", 
                    "eaear_Oilfeed_crops", "eaear_Sugarcane", "eaear_Oilpalm", "eaear_Tobacco",
                    "eaear_Cocoa_Coffee")


### MAIN DATA SET ### 
brazil_data <- readRDS(here("temp_data", "merged_datasets", "brazil_aoi", "loss_commodity_aeay_long_final.Rdata"))
southam_data <- readRDS(here("temp_data", "merged_datasets", "southam_aoi", "loss_commodity_aeay_long_final.Rdata"))
seasia_data <- readRDS(here("temp_data", "merged_datasets", "SEAsia_aoi", "loss_commodity_aeay_long_final.Rdata"))

### PRICE DATA ###
prices <- readRDS(here("temp_data", "prepared_international_prices.Rdata"))

# renewable fuel standards, as from https://www.epa.gov/renewable-fuel-standard-program/renewable-fuel-annual-standards
# remember: biodiesel are nested within advanced biofuels, which are themeselves nested within total
# conventional biofuels are the part of the total that is not "advanced"
rfs <- data.frame(year = min(prices$year):2022, 
                  statute_total = 0, final_total = 0, 
                  statute_advanced = 0, final_advanced = 0, 
                  statute_biodiesel = 0, final_biodiesel = 0, 
                  statute_conv_earliest = 0)

rfs <- dplyr::arrange(rfs, year)
# this is if RFS2 takes over as soon as 2008 (9bgal)
rfs[rfs$year >= 2006 & rfs$year <= 2022, c("statute_total")] <- c(4, 4.7, 9, 11.1, 12.95,	13.95,	15.2,	16.55,	18.15,	20.5,	22.25,	24.0,	26.0,	28.0, 30.0, 33.0, 36.0)
rfs[rfs$year >= 2006 & rfs$year <= 2022, c("final_total")] <- c(4, 4, 9, 11.1, 12.95,	13.95,	15.2,	16.55,	16.28,	16.93,	18.11,	19.28,	19.29,	19.92, 20.09, NA, NA)

rfs[rfs$year >= 2006 & rfs$year <= 2022, c("statute_advanced")] <- c(0, 0, 0, 0.6, 0.95,	1.35,	2.0,	2.75,	3.75,	5.5,	7.25,	9.0,	11.0,	13.0, 15.0, 18.0, 21.0)
rfs[rfs$year >= 2006 & rfs$year <= 2022, c("final_advanced")] <- c(0, 0, 0, 0.6, 0.95,	1.35,	2.0,	2.75,	2.67,	2.88,	3.61,	4.28,	4.29,	4.92, 5.09, NA, NA)

rfs[rfs$year >= 2009 & rfs$year <= 2022, c("final_biodiesel")] <- c(0.5, 1.15, 0.8, 1.0, 1.28, 1.63, 1.73, 1.9, 2.0, 2.1, 2.1, 2.43, 2.43, NA)
rfs[rfs$year >= 2009 & rfs$year <= 2022, c("statute_biodiesel")] <- c(0.5, 0.65, 0.8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

rfs <- mutate(rfs, statute_conv = statute_total - statute_advanced)
rfs <- mutate(rfs, final_conv = final_total - final_advanced)
# RFS mandates with for each year, the earliest set target (by RFS1 up to 2012 included, and by RFS2 afterwards)
rfs[rfs$year >= 2006 & rfs$year <= 2022, c("statute_conv_earliest")] <- c(4, 4.7, 5.4, 6.1, 6.8, 7.4, 7.5, 13.8, 14.4, 15, 15, 15, 15, 15, 15, 15, 15)
#rfs[rfs$year >= 2006 & rfs$year <= 2020, c("statute_conv")] <- c(4, 4, 9, 10.5, 12, 12.6, 13.2, 13.8, 14.4, 15, 15, 15, 15, 15, 15)

# rfs <- mutate(rfs, amendments = final_conv - statute_conv) # in 2014 it ramped up to 13.61 instead of 14.4, which is thus a negative shock on demand of .79 bgal
# rfs <- round(rfs, 2)
# rfs$blendwall <- 0
# rfs[rfs$year == 2014, "blendwall"] <- -0.79

# make random permutations of main rfs vars as placebos
# rfs <- mutate(rfs, statute_conv_placebo = sample(statute_conv, length(statute_conv)), 
#               statute_biodiesel_placebo = sample(statute_biodiesel, length(statute_biodiesel)))

# merge with price time series
prices <- full_join(prices, rfs, by = "year")

# make lags and leads
rfs_vars <- c("statute_conv", "statute_conv_earliest")#, "amendments", "blendwall"
for(voi in rfs_vars){
  for(lead in c(1:6)){
    prices <- dplyr::arrange(prices, year)
    prices <- DataCombine::slide(prices,
                                 Var = voi, 
                                 TimeVar = "year",
                                 NewVar = paste0(voi,"_lead",lead),
                                 slideBy = lead, 
                                 keepInvalid = FALSE)
    prices <- dplyr::arrange(prices, year)
    
    # up to 2005, leads are actually 0 because mandates were not known/disclosed
    prices[prices$year <= 2005, paste0(voi,"_lead", lead)] <- 0
    # and in 2006 the expected mandates are those of the RFS1 (the EISA was discussed in Senate as of January 2007)
    prices[prices$year == 2006, paste0(voi,"_lead", lead)] <- prices[prices$year == 2006+lead, "statute_conv_earliest"]
  }  
  
  for(lag in c(1:6)){
    prices <- dplyr::arrange(prices, year)
    prices <- DataCombine::slide(prices,
                                 Var = voi, 
                                 TimeVar = "year",
                                 NewVar = paste0(voi,"_lag",lag),
                                 slideBy = -lag, 
                                 keepInvalid = FALSE)
    prices <- dplyr::arrange(prices, year)
  }  
  
  
  for(py in c(1:6)){
    ## Future-year averages (1, 2, 3, and 4 years, after and EXCLUDING contemporaneous)  
    # note that we DON'T add voi column (not lagged) in the row mean
    # prices$newv <- rowMeans(x = prices[,c(voi, paste0(voi,"_lead",c(1:py)))], na.rm = FALSE)
    # prices[is.nan(prices$newv),"newv"] <- NA
    # colnames(prices)[colnames(prices)=="newv"] <- paste0(voi,"_",py,"fya")
    prices <- mutate(prices, 
                     !!as.symbol(paste0(voi,"_",py,"fya")) := round(rowMeans(across(.cols = any_of(paste0(voi,"_lead",c(1:py)))), 
                                                                             na.rm = TRUE), 2) 
    )
    ## Past-year averages (1, 2, 3, and 4 years, before and INCLUDING contemporaneous)  
    # note that we DO add voi column (not lagged) in the row mean
    prices <- mutate(prices, 
                     !!as.symbol(paste0(voi,"_",py,"pya")) := round(rowMeans(across(.cols = any_of(c(voi, paste0(voi,"_lag",c(1:py))))), 
                                                                             na.rm = TRUE), 2)
    )
    
    
    # AS A CONSEQUENCE, IT IS NORMAL THAT ABS. EXPOSURE CONTROLS ARE INTERACTED WITH 1fya AND 2pya: BOTH ARE AVERAGES OF TWO YEARS
  }
}

# prices[,c("year", grep(names(prices), pattern = "statute_conv", value = TRUE))] %>% tail(25)


#### REGRESSION FUNCTION #### 

### TEMPORARY OBJECTS 
outcome_variable = "loss_commodity" # "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
endo_vars = "unidir_soy"
start_year = 2011
end_year = 2019
region = "southam"

pasture_shares <- FALSE
standardization = ""

rfs_rando <- ""
original_rfs_treatments <- c("statute_conv")
rfs_lag <- 2
rfs_lead <- 3
rfs_fya <-  0
rfs_pya <- 0
lag_controls = NULL
aggr_dyn <- TRUE
exposure_rfs <- eaear_mapmat[,"Crops"]#"eaear_Oilpalm"
group_exposure_rfs <- FALSE
control_all_absolute_rfs <- FALSE
most_correlated_only = FALSE
annual_rfs_controls <- FALSE

control_pasture <- FALSE
pasture_trend <- FALSE

fc_trend <- FALSE
s_trend <- FALSE
s_trend_loga <- FALSE
fc_s_trend <- FALSE

sjpos <- FALSE # should the sample be restricted to cells where sj is positive? 

fe = "grid_id + country_year" #   
preclean_level = "FE"
distribution_1st = "quasipoisson"
distribution_2nd = "quasipoisson"
invhypsin = TRUE
conley_cutoff <- 100
clustering = "oneway" # either "oneway" or "twoway". If oneway, it clusters on cluster_var1. 
cluster_var1 = "grid_id_10" 
cluster_var2 = "grid_id_5_year"
output = "est_object"
glm_iter <- 25 

stddev = FALSE # if TRUE, the PEs are computed for a one standard deviation (after removing variation in the fixed-effect dimensions)
rel_lu_change = 0.01
abs_lu_change = 1 
rounding = 2

rm(outcome_variable, start_year, end_year, pasture_shares, 
   standardization, group_prices, biofuel_focus, aggregate_K, control_interact_sj, control_interact_Pk, reference_crop, control_direct, 
   rfs_reg, original_rfs_treatments, rfs_lead, rfs_lag, original_exposure_rfs, group_exposure_rfs, control_absolute_rfs, control_all_absolute_rfs, sjpos, fe, 
   distribution_1st, distribution_2nd, invhypsin, conley_cutoff, se, boot_cluster, coefs_to_aggregate, 
   output, glm_iter)




make_main_reg <- function(region = "southam", # one of "brazil", "southam", or "seasia" 
                          endo_vars = "unidir_soy", 
                          
                          outcome_variable = "loss_commodity", # one of "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
                          start_year = 2011, 
                          end_year = 2019, 
                          continent = "all", # one of "Africa", "America", "Asia", or "all"
                          
                          rfs_rando = "", # either "between", "within", or any other string. If one of the former two, randomization inference of the type is performed
                          original_rfs_treatments = c("statute_conv"),
                          rfs_lead = 0,
                          rfs_lag = 0,
                          rfs_fya = 0, 
                          rfs_pya = 0,
                          lag_controls = NULL, # which of the lags specified above should be CONTROLLED for rather than counted in the cumulative effect
                          aggr_dyn = TRUE, # whether to report aggregate coefficients of all leads and lags ("all") or leads and lags separately ("leadlag"), or no aggregate (any other string)
                          exposure_rfs = "eaear_Soy_compo",
                          all_exposures_rfs = eaear_mapmat[,"Crops"],
                          group_exposure_rfs = FALSE, # This is deprecated, as it was relevant only when rfs exposures where standardized and could thus be added. If exposure_rfs is length 1, it does not matter whether this option is TRUE or FALSE.
                          control_all_absolute_rfs = TRUE, # there is not really an alternative where this could be FALSE and have sense. 
                          most_correlated_only = FALSE, # but this restricts the controls to only interactions with the most correlated crop. 
                          annual_rfs_controls = FALSE,
                          
                          # second stage dynamics
                          annual = FALSE, 
                          outcome_lags = 0, # 0 means current deforestation is regressed on current expansion. 1 means deforestation one year ago is regressed on current expansion, etc. 
                          
                          control_pasture = FALSE,
                          pasture_trend = FALSE,
                          remaining = FALSE, # should remaining forest be controlled for STOP DOING THIS BECAUSE IT INTRODUCES NICKELL BIAS
                          s_trend = TRUE,
                          s_trend_loga = FALSE,
                          fc_trend = FALSE,
                          fc_s_trend = FALSE,
                          
                          sjpos = FALSE, # should the sample be restricted to cells where sj is positive? 
                          
                          fe = "grid_id + country_year", 
                          preclean_level = "FE", 
                          distribution_1st = "quasipoisson",#  "quasipoisson", 
                          distribution_2nd = "quasipoisson",#  "quasipoisson", 
                          invhypsin = TRUE, # if distribution is gaussian, should the dep. var. be transformed to inverse hyperbolic sine?
                          
                          clustering = "oneway", # either "oneway" or "twoway". If oneway, it clusters on cluster_var1. 
                          cluster_var1 = "grid_id_10", 
                          cluster_var2 = "grid_id_10",

                          glm_iter = 25,
                          # dyn_tests = FALSE, # should the Fisher-type panel unit root test be returned, instead of the regressions, for the outcome_variable and the first regressor
                          boot_rep = 500, # number of bootstrap replicates to compute APE SE. 

                          # parameters for APE computation 
                          stddev = FALSE, # if TRUE, the PEs are computed for a one standard deviation (after removing variation in the fixed-effect dimensions)
                          rel_lu_change = 0.01, 
                          abs_lu_change = 1, 
                          rounding = 2,
                          
                          output = "coef_table" # one of "data", est_object, or "coef_table" 
){
  
  # Define the outcome_variable based on the crop under study (if we are not in the placebo case)
  if(outcome_variable == "tmf_agri" & ("eaear_Oilpalm" %in% exposure_rfs | "eaear_Rubber" %in% exposure_rfs)){
    outcome_variable <- "tmf_plantation"
  }
  
  #### PREPARE NEEDED VARIABLE NAMES
  # this does not involve data, just arguments of the make_reg function
  # original_ names are used to get generic covariate names in coefficient tables (irrespective of modelling choices)
  
  
  all_rfs_treatments <- grep(pattern = original_rfs_treatments, names(prices), value = TRUE)
  
  rfs_treatments <- original_rfs_treatments
  if(rfs_lead >= 1){
    for(rfs_var in original_rfs_treatments){
      rfs_treatments <- c(rfs_treatments, paste0(rfs_var, "_lead", 1:rfs_lead))
    }
    #  for(rfs_var in original_rfs_treatments){
    #   rfs_treatments <- paste0(rfs_var, "_lead", rfs_lead)
    # }
  }
  if(rfs_lag >= 1){
    for(rfs_var in original_rfs_treatments){
      rfs_treatments <- c(rfs_treatments, paste0(rfs_var, "_lag", 1:rfs_lag))
    }
    # for(rfs_var in original_rfs_treatments){
    #   rfs_treatments <- paste0(rfs_var, "_lag", rfs_lag)
    # }
  }
  
  # Or they are past or future values averaged over a certain amount of time
  if(rfs_fya >= 1 & rfs_pya == 0){
    for(rfs_var in original_rfs_treatments){
      rfs_treatments <- paste0(rfs_var, "_", rfs_fya, "fya") 
    }
  }  
  if(rfs_pya >= 1 & rfs_fya == 0){
    for(rfs_var in original_rfs_treatments){
      rfs_treatments <- paste0(rfs_var, "_", rfs_pya, "pya") 
    }
  }
  if(rfs_fya >=1 & rfs_pya >=1){
    for(rfs_var in original_rfs_treatments){
      # ORDER MATTERS
      rfs_treatments <- c(paste0(rfs_var, "_", rfs_fya, "fya"), paste0(rfs_var, "_", rfs_pya, "pya") )
    }
  }
  
  # or they are a combination of averaged past treatments, and annual future ones
  if(rfs_pya > 0 & rfs_fya == 0 & rfs_lead > 0 & rfs_lag == 0){
    for(rfs_var in original_rfs_treatments){
      rfs_treatments <- c(paste0(rfs_var, "_lead", 1:rfs_lead), paste0(rfs_var, "_", rfs_pya, "pya") )
    }
  }
  
  #### MAKE THE VARIABLES NEEDED IN THE DATA
  # manipulate a different data set so that original one can be provided to all functions and not read again every time. 
  if(region == "brazil"){
    d <- brazil_data
  }
  
  if(region == "southam"){
    d <- southam_data
  }
  
  if(region == "seasia"){
    d <- seasia_data
  }
  
  if(grepl("tmf_", outcome_variable)){
    d  <- readRDS(here("temp_data", "merged_datasets", "tmf_aoi", "tmf_aeay_pantrop_long_final_1990_2020.Rdata"))
    # release some memory upfront
    d <- dplyr::filter(d, year >= 2008, year <= 2019)
    
    d <- dplyr::mutate(d, tmf_deforestation = tmf_agri + tmf_plantation)
  }
  
  # code below should work whether d is from tmf or losscommo 
  # Keep only in data the useful variables 
  d <- dplyr::select(d, all_of(unique(c("grid_id", "year", "lat", "lon", "continent_name", "country_name", "country_year",
                                        # "fc_2000", "fc_2009", "remaining_fc", # "accu_defo_since2k",
                                        "grid_id_5", "grid_id_10", "grid_id_20", "grid_id_5_year", "grid_id_10_year", "grid_id_20_year",
                                        outcome_variable, # "tmf_agri", "tmf_flood", "tmf_plantation",
                                        endo_vars,
                                        "pasture_share_2000",
                                        eaear_mapmat[,"Crops"], exposure_rfs )))) #sj, 
  
  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", unique(c(rfs_treatments, all_rfs_treatments)))], by = c("year"))#, all_treatments
  
  
  # If we want the exposure to deforestation for pasture to be proxied with the share of pasture area in 2000, rather than suitability index, 
  if(pasture_shares){
    d[,grepl("Fodder", names(d))] <- d$pasture_share_2000
  }
  
  if((distribution_2nd == "gaussian") & invhypsin){
    # transform dependent variable, if gaussian GLM 
    d <- dplyr::mutate(d, !!as.symbol(outcome_variable) := asinh(!!as.symbol(outcome_variable)))
  }
 
  
  #### INSTRUMENTS #### 
  
  instruments <- c()
  for(rfs_var in rfs_treatments){
    if(group_exposure_rfs){
      # group (sum) exposures
      d <- dplyr::mutate(d, grp_exp = rowSums(across(.cols = (any_of(exposure_rfs)))))
      # make the regressor of interest
      varname <- paste0("grp_exp_X_", rfs_var)
      instruments <- c(instruments, varname)
      d <- mutate(d, !!as.symbol(varname) := grp_exp * !!as.symbol(rfs_var))
    }else{
      for(exp_rfs in exposure_rfs){
        # make instruments of interest
        varname <- paste0(exp_rfs, "_X_", rfs_var)
        instruments <- c(instruments, varname)
        d <- mutate(d, !!as.symbol(varname) := !!as.symbol(exp_rfs) * !!as.symbol(rfs_var))
      }
    }
  }
  
  ### CONTROLS ####
  # it's important that this is not conditioned on anything so these objects exist
  controls <- c()
  
  controls <- c(controls, instruments)
  
  # remove lags of the treatment of interest from instruments (i.e. transfer them to controls)
  if(length(lag_controls)>0){
    lags_to_transfer <- c()
    for(rfs_var in original_rfs_treatments){
      lags_to_transfer <- c(lags_to_transfer, paste0(exposure_rfs, "_X_", rfs_var, "_lag", lag_controls))
    }
    
    instruments <- instruments[instruments != lags_to_transfer]
    
    # controls <- c(controls, lags_to_transfer) # not in IV context
    # no need to change their names, they are not captured by cumulative effects as long as they are in controls. 
    
    # old_controls <- c()
    # for(old_rc in lags_to_transfer){
    #   varname <- paste0("jexp_X_", old_rc)
    #   controls <- c(controls, varname)
    #   d <- mutate(d, !!as.symbol(varname) := !!as.symbol(exposure_rfs) * !!as.symbol(old_rc))
    # }
  }
  
  # IT IS NORMAL THAT ABS. EXPOSURE CONTROLS ARE INTERACTED WITH 2fya AND 1pya: BOTH ARE AVERAGES OF TWO YEARS. 
  # 1pya is the avg of current and lag1 values, which are grouped together because they reflect the same kind of mechanism. 
  if(control_all_absolute_rfs){
    all_abs <- eaear_mapmat[, "Crops"][!(eaear_mapmat[, "Crops"] %in% original_exposure_rfs)]
    
    # add the share of pasture in grid cell in 2000 as an exposure to pastures (if the exposure of interest is not pasture)
    if(control_pasture & !grepl("Fodder", original_exposure_rfs)){
      all_abs <- c(all_abs, "pasture_share_2000")
    }
    
    for(abs in all_abs){
      # if there are several rfs treatments due to lags and leads, then we interact exposures to other crops with an average of those only, not with each
      if(rfs_lead>=1 & !annual_rfs_controls){
        rfs_control <- paste0(original_rfs_treatments, "_", rfs_lead, "fya")  
        varname <- paste0(abs, "_X_", rfs_control)
        controls <- c(controls, varname)
        d <- mutate(d, !!as.symbol(varname) := !!as.symbol(abs) * !!as.symbol(rfs_control))
      }
      if(rfs_lag>=1 & !annual_rfs_controls){
        rfs_control <- paste0(original_rfs_treatments, "_", rfs_lag, "pya")  
        varname <- paste0(abs, "_X_", rfs_control)
        controls <- c(controls, varname)
        d <- mutate(d, !!as.symbol(varname) := !!as.symbol(abs) * !!as.symbol(rfs_control))
      }
      if(rfs_lag==0 & rfs_lead>0 & rfs_pya > 0 & rfs_fya ==0){
        rfs_control <- paste0(original_rfs_treatments, "_", rfs_pya, "pya")  
        varname <- paste0(abs, "_X_", rfs_control)
        controls <- c(controls, varname)
        d <- mutate(d, !!as.symbol(varname) := !!as.symbol(abs) * !!as.symbol(rfs_control))
      }
      # note the use of the original_rfs_treatments vector: if no lag is specified, we still want to control for interactions with contemporaneous RFS
      # but this should never happen
      if(rfs_lag==0 & rfs_pya == 0 & !annual_rfs_controls){
        for(ogn_rfs_var in original_rfs_treatments){
          varname <- paste0(abs, "_X_", ogn_rfs_var)
          controls <- c(controls, varname)
          d <- mutate(d, !!as.symbol(varname) := !!as.symbol(abs) * !!as.symbol(ogn_rfs_var))
        }
      }
      if(rfs_lead==0 & rfs_lag==0 & !annual_rfs_controls){
        for(rfs_var in rfs_treatments){
          varname <- paste0(abs, "_X_", rfs_var)
          controls <- c(controls, varname)
          d <- mutate(d, !!as.symbol(varname) := !!as.symbol(abs) * !!as.symbol(rfs_var))
        }    
      }
      
      
      # unless we specify that we want the full range of interactions as controls, i.e. all exposures with all leads and lags
      if(annual_rfs_controls){
        for(rfs_var in rfs_treatments){
          varname <- paste0(abs, "_X_", rfs_var)
          controls <- c(controls, varname)
          d <- mutate(d, !!as.symbol(varname) := !!as.symbol(abs) * !!as.symbol(rfs_var))
        }
      }
    }  
  }
  
  
  # add pasture share trend
  if(pasture_trend){
    d <- mutate(d, pasture_share_trend = pasture_share_2000 * year)
    controls <- c(controls, "pasture_share_trend")
  }
  
  if(fc_trend){
    # the conditions are just for the sake of computation efficiency
    if(start_year != 2001 & start_year != 2010){
      fc_year <- d[d$year == start_year - 1, c("grid_id", "remaining_fc")]
      names(fc_year) <- c("grid_id", paste0("fc_",start_year-1))
      d <- left_join(d, fc_year, by = "grid_id")
      d <- mutate(d, forest_cover_trend := !!as.symbol(paste0("fc_",start_year-1))  * year)
    }
    if(start_year == 2001){
      d <- mutate(d, forest_cover_trend = fc_2000 * year)
    }
    if(start_year == 2010){
      d <- mutate(d, forest_cover_trend = fc_2009 * year)
    }
    
    controls <- c(controls, "forest_cover_trend")
  }
  
  if(s_trend){
    # Suitability time trend
    # d <- mutate(d, suitability_trend := !!as.symbol(sj)*(year))#-2000 # doesn't change anything that it's multiplied by 1...19 or 2001...2019
    # controls <- c(controls, "suitability_trend")
    for(eaear_exp_rfs in exposure_rfs){
      varname <- paste0(eaear_exp_rfs, "_trend")
      controls <- c(controls, varname)
      d <- mutate(d, !!as.symbol(varname) := !!as.symbol(eaear_exp_rfs) * (year-2000))
    }    
  }
  
  if(fc_s_trend){
    # d <- mutate(d, forest_cover_suitability_trend := !!as.symbol(sj)*fc_2000*(year))#-2000 # doesn't change anything that it's multiplied by 1...19 or 2001...2019
    # controls <- c(controls, "forest_cover_suitability_trend")
    for(eaear_exp_rfs in exposure_rfs){
      varname <- paste0(eaear_exp_rfs, "_fc_trend")
      controls <- c(controls, varname)
      d <- mutate(d, !!as.symbol(varname) := !!as.symbol(eaear_exp_rfs) * fc_2000 * (year-2000))
    }  
  }
  

  ### KEEP OBSERVATIONS THAT: ####
  
  # - are suitable to crop j 
  if(sjpos){
    d <- dplyr::filter(d, !!as.symbol(sj) > 0)  
  }
  
  # - are in study period 
  if(start_year != 2001 | end_year != 2019){
    d <- dplyr::filter(d, year >= start_year)
    d <- dplyr::filter(d, year <= end_year)
  }

  # remove units with no country name (in the case of all_drivers data set currently, because nearest_feature function has not been used in this case, see add_variables.R)
  d <- dplyr::filter(d, !is.na(country_name))

  used_vars <- unique(c("grid_id", "year", "lat", "lon","continent_name",  "country_year",  #"country_name",  "remaining_fc", "accu_defo_since2k", # "sj_year",
                        "grid_id_5", "grid_id_10", "grid_id_20", "grid_id_5_year", "grid_id_10_year", "grid_id_20_year",
                        outcome_variable,# "tmf_agri", "tmf_flood", "tmf_plantation",
                        endo_vars, instruments, controls, 
                        exposure_rfs, original_rfs_treatments, rfs_treatments))    # this is necessary to reconstruct variables in randomization inference processes

  # - have no NA nor INF on any of the variables used (otherwise they get removed by {fixest})
  # for instance, there are some NAs in the suitability index (places in water that we kept while processing other variables...) 
  usable <- lapply(used_vars, FUN = function(var){is.finite(d[,var]) | is.character(d[,var])})
  # used_vars[!(used_vars%in%names(d))]
  names(usable) <- used_vars            
  usable <- bind_cols(usable)
  filter_vec <- base::rowSums(usable)
  filter_vec <- filter_vec == length(used_vars)
  d <- d[filter_vec, c(used_vars)]
  if(anyNA(d)){stop()}
  rm(filter_vec, usable)
  
  # is.na(d$Oilpalm) 
  
  # experience deforestation at least once (automatically removed if distribution is quasipoisson, but not if it's gaussian. 
  # Though we want identical samples to compare distributional assumptions, or FE specifications
  if(preclean_level == "FE"){
    preclean_level <- fe
  }
  temp_est <- feglm(fml = as.formula(paste0(outcome_variable, " ~ 1 | ", preclean_level)),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    d_clean <- d[unlist(temp_est$obs_selection),]
  }  else { 
    d_clean <- d
  }
  
  # if we don't remove obs that have no remaining forest, and if we don't remove always 0 outcome obs. in previous step, d IS A BALANCED PANEL !!
  # this is a matter only if we do beta process, with cluster bootstrap
  # if(!nrow(d) == length(unique(d$grid_id))*length(unique(d$year))){
  #   warning("data is not balanced")
  # }
  rm(d)
  
  # save these for later
  avg_defo_ha <- temp_est$fitted.values %>% mean() %>% round(1)
  
  G <- length(unique(d_clean[,cluster_var1]))
  
  d_clean_save <- d_clean
  
  ## Infrastructure to store outputs from estimations below
  toreturn <- list(firstg_summary = NA, 
                   boot_info = NA,
                   APE_mat = NA)
  output_list <- list()
  
  ### FORMULAE ####
  
  # INTRUMENTS ARE IN CONTROLS   
  # store first stage formulae in this list 
  fml_1st_list <- list()
  classic_iv_fixest_fml <- list()
  for(f in 1:length(endo_vars)){
    fml_1st_list[[f]] <- as.formula(paste0(endo_vars[f],
                                           " ~ ", 
                                           paste0(controls, collapse = "+"), #  instruments are in controls
                                           " | ", 
                                           fe))
    
    # this is not for within bootstraps, but to fit in feols-iv in fixest package to extract 1st stage stats readily
    if(length(controls[!(controls %in% instruments)] > 0)){
      excluded_controls <- controls[!(controls %in% instruments)]
    }else { 
      excluded_controls <- 1 
    }
    
    classic_iv_fixest_fml[[f]] <- as.formula(paste0(outcome_variable,
                                                    " ~ ", 
                                                    paste0(excluded_controls, collapse = "+"), # only "included" controls, i.e. z1, i.e. not instruments.
                                                    " | ",
                                                    fe, 
                                                    " | ",
                                                    endo_vars[f], 
                                                    " ~ ",
                                                    paste0(instruments, collapse = "+")))
    
  }
  
  
  # need to specify the names of first stages' residuals now
  errors_1stg <- paste0("est_res_1st_", c(1:length(endo_vars)))
  
  fml_2nd <- as.formula(paste0(outcome_variable,
                               " ~ ", 
                               paste0(endo_vars, collapse = "+"),
                               " + ", 
                               paste0(errors_1stg, collapse = "+"),
                               " + ",
                               paste0(excluded_controls, collapse = "+"), # only "included" controls, i.e. z1, i.e. not instruments.
                               " | ",
                               fe))
  ## ESTIMATION 
  se <- as.formula(paste0("~ ", paste0(c(cluster_var1, cluster_var2), collapse = "+")))
  
  if(clustering =="oneway"){
    se <- as.formula(paste0("~ ", cluster_var1))
  }
  
  # We run the first stage outside the bootstrapping process, in order to extract related information
  est_1st_list <- list()
  Ftest_list <- list()
  Waldtest_list <- list()
  for(f in 1:length(endo_vars)){
    tsls <- fixest::feols(classic_iv_fixest_fml[[f]], 
                          data = d_clean)
    
    # Extract first stage information
    # est_1st <- summary(tsls, stage = 1)
    # est_1st_list[[f]] <- summary(est_1st, .vcov = vcov(est_1st, cluster = CLUSTER))#se = "cluster", cluster = CLUSTER
    
    Ftest_list[[f]] <- fitstat(tsls, type = "ivf1")
    Waldtest_list[[f]] <- fitstat(tsls, type = "ivwald1")
    
    est_1st_list[[f]] <- summary(fixest::feols(fml_1st_list[[f]], d_clean), cluster = se)
    
    # save residuals 
    # d_clean[,errors_1stg[f]] <- summary(est_1st, stage = 1)$residuals 
  }

  toreturn[["firstg_summary"]] <- est_1st_list
  
  Fstat <- unlist(Ftest_list)[[1]] %>% round(2)
  Waldstat <- unlist(Waldtest_list)[[1]] %>% as.numeric() %>% round(2)

  myfun_data <- d_clean
  fsf_list <- fml_1st_list
  ssf <- fml_2nd
  num_iter <- 25
  
ctrl_fun_endo <- function(myfun_data, fsf_list, ssf, num_iter = 25){
  
  statistics <- tryCatch(expr = {
    # As many first stages as there are endogenous variables
    for(f in 1:length(endo_vars)){
      BS_est_1st <- fixest::feols(fsf_list[[f]], 
                                  data = myfun_data)
      
      # save residuals 
      myfun_data[,errors_1stg[f]] <- BS_est_1st$residuals
      rm(BS_est_1st)
    }
    
    # 2nd stage
    # increase glm algorithm number of iterations... iteratively 
    # i.e. ask for more glm iterations only if it did not converge with fewer. 
    # assign("warning", NULL)
    # actually, this took more than 24h to run 500 bootstraps for all 9 plantation groups. 
    # SO we rather impose low # of glm.iter, catch out all those that don't converge, and perform more bootstraps instead. 
    
    #while(num_iter <= 25 ) {
    
    BS_est_2nd <-   fixest::feglm(ssf, 
                                  data = myfun_data, 
                                  family = distribution_2nd, 
                                  glm.iter = num_iter,
                                  notes = TRUE) 
    
    
    # # if feglm did not converge, try with more iterations
    # if(class(BS_est_2nd) == "warnings"){
    #   num_iter <- num_iter + 500
    # 
    # } else { # else, finish the estimation
    #   num_iter <- 3000
    # }
    #}
    
    
    ## MAKE APE 
    # Unlike in the non-IV case, APEs are not estimated with the delta Method, because their SEs are only computed asymptotically, by bootstrap. 
    # we don't need to compute inference statistics on them, as inference is handled by bootstrap (on APE).
    
    # select coefficient of interest 
    roi <- BS_est_2nd$coefficients[endo_vars]
    
    boot_output <- c()
    i <- 1
    for(r in roi){
      # the APE is different depending on the regressor of interest being in the log scale or not. 
      if(grepl("ln_",endo_vars[1])){
        boot_output[i] <- ((1+rel_lu_change)^(r) - 1)*100
      } else{
        boot_output[i] <- (exp(r*abs_lu_change) - 1)*100 
      } 
      i <- i+1
    }
    
    # control function variables 
    cfreg <- BS_est_2nd$coefficients[errors_1stg]
    
    # (i is already taken to the +1 level, so index i is empty at this point)
    boot_output[i] <- sum(cfreg)
    
    boot_output
    
  }, # closes tryCatch expr argument  
  warning = function(w) print("warning in feglm captured"), #
  error = function(e) print("error in feglm captured"), 
  finally = "0" ) # ça sert à rien, mais il faut passer quelque chose 
  
  # if, for this boot sample, fixest returns an error in either stage, or the second stage returns a warning, 
  # then statistics is resp. an error or NULL, in which case we coerce it to NA 
  if(!is.numeric(statistics)){
    # this weird formula makes final_stats the right length, either 1 if annual = FALSE, 
    # and, in the other case, the number of years the annual effects are estimated plus one for the cumulative effect, plus one for the 1st stage residual
    final_stats <- rep(NA, annual*(outcome_lags) + 1 + 1)
  } else {
    # return both cumulative and annual APE as statistics to bootstrap, as well as the (cumulative) coeffs on control function terms
    # which are in statistics
    if(annual){ # PAS POSSIBLE DE FAIRE DU ANNUAL SI LE REFERENTIEL EST ENDO_VAR. 
      final_stats <- c(sum(statistics[1:length(endo_vars)]), statistics) 
    } else { 
      final_stats <- statistics
    }        
  }
  rm(statistics)      
  
  # statistics we want to evaluate the variance of:
  return(final_stats)
}



  # We run the first stage outside the bootstrapping process, in order to extract related information
  est_1st <- fixest::feglm(fml_1st, 
                           data = d_clean, 
                           vcov = se,
                           family = "quasipoisson",
                           fixef.rm = "perfect",
                           nthreads = 3, 
                           notes = TRUE)
  # this is the data actually used in the estimation of the first stage, not necessarily the same as input data d_clean because of units having always zero endogenous var. 
  d_clean_1st <- d_clean[est_1st$obs_selection[[1]], ]
  
  # est_1stgauss <- fixest::feglm(fml_1st, 
  #                            data = d_clean_1st, 
  #                            vcov = se,
  #                            family = "gaussian",
  #                            fixef.rm = "perfect",
  #                            nthreads = 3, 
  #                            notes = TRUE)
  # 
  # # compare fits 
  # summary(est_1st$fitted.values)
  # summary(est_1stgauss$fitted.values)
  # summary(d_clean_1st[,endo_vars])
  # # quasipoisson does not predict negative values, contrary to gaussian
  # # all means are the same. Maximum is slightly better predicted by gaussian
  
  # those two lines necessary to compute aggregated effects below
  df_res_1st <- summary(est_1st)$coeftable
  
  fixest_df_1st <- degrees_freedom(est_1st, type = "t")
  
  ## Extract first stage information
  if(aggr_dyn & rfs_lead > 0 & rfs_lag > 0 & rfs_fya == 0 & rfs_pya == 0){
    # In this case, we are interested in LEAD AND LAG effects, aggregated separately, and all together.
    df_res_1st <- rbind(rep(NA, ncol(df_res_1st)), rep(NA, ncol(df_res_1st)), df_res_1st)
    
    # ORDER MATTERS
    aggr_names <- paste0(original_exposure_rfs, c("_X_aggrleads", "_X_aggrlags"))
    row.names(df_res_1st)[1:2] <- aggr_names
    # instruments of interest
    base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
    # Contemporaneous value is NOT in leads
    lead_roi <- c(grep(pattern = paste0(base_reg_name,"_lead"), 
                       instruments, value = TRUE))
    # Contemporaneous value is in lags
    lag_roi <- c(base_reg_name, 
                 grep(pattern = paste0(base_reg_name,"_lag"), 
                      instruments, value = TRUE))
    
    df_res_1st[aggr_names[1],"Estimate"] <- est_1st$coefficients[lead_roi] %>% sum()
    df_res_1st[aggr_names[2],"Estimate"] <- est_1st$coefficients[lag_roi] %>% sum()
    
    # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
    # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
    df_res_1st[aggr_names[1],"Std. Error"] <- est_1st$cov.scaled[lead_roi, lead_roi] %>% as.matrix() %>% sum() %>% sqrt()
    df_res_1st[aggr_names[2],"Std. Error"] <- est_1st$cov.scaled[lag_roi, lag_roi] %>% as.matrix() %>% sum() %>% sqrt()
    
    df_res_1st[aggr_names[1],"t value"]  <- (df_res_1st[aggr_names[1],"Estimate"] - 0)/(df_res_1st[aggr_names[1],"Std. Error"])
    df_res_1st[aggr_names[2],"t value"]  <- (df_res_1st[aggr_names[2],"Estimate"] - 0)/(df_res_1st[aggr_names[2],"Std. Error"])
    
    # use t distribution with degrees of freedom equal to that used by fixest, i.e. after two way cluster adjustment.
    # does not make a significant difference given sample size
    df_res_1st[aggr_names[1],"Pr(>|t|)"]  <- (2*pt(abs(df_res_1st[aggr_names[1],"t value"]), 
                                                   lower.tail = FALSE, 
                                                   df = fixest_df_1st)) 
    df_res_1st[aggr_names[2],"Pr(>|t|)"]  <- (2*pt(abs(df_res_1st[aggr_names[2],"t value"]), 
                                                   lower.tail = FALSE, 
                                                   df = fixest_df_1st)) 
    
    
    ## Then (on top of dataframe), aggregate all leads and lags together
    # it's important that overall aggregate comes after, for at least two reasons in current code:
    # 1. because selection of estimate to plot is based on order: it takes the first row of df_res_1st 
    # 2. so that aggr_names corresponds to *_X_aggrall in randomization inference below
    
    df_res_1st <- rbind(rep(NA, ncol(df_res_1st)), df_res_1st)
    
    # ORDER MATTERS
    aggr_names <- paste0(original_exposure_rfs, c("_X_aggrall"))
    row.names(df_res_1st)[1] <- aggr_names
    # instruments of interest
    base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
    
    # Contemporaneous, lead, and lag values
    all_roi <- c(base_reg_name, 
                 grep(pattern = paste0(base_reg_name,"_lead"), 
                      instruments, value = TRUE),
                 grep(pattern = paste0(base_reg_name,"_lag"), 
                      instruments, value = TRUE))
    
    df_res_1st[aggr_names[1],"Estimate"] <- est_1st$coefficients[all_roi] %>% sum()
    
    # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
    # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
    df_res_1st[aggr_names[1],"Std. Error"] <- est_1st$cov.scaled[all_roi, all_roi] %>% as.matrix() %>% sum() %>% sqrt()
    
    df_res_1st[aggr_names[1],"t value"]  <- (df_res_1st[aggr_names[1],"Estimate"] - 0)/(df_res_1st[aggr_names[1],"Std. Error"])
    
    # use t distribution with degrees of freedom equal to that used by fixest, i.e. after two way cluster adjustment.
    # does not make a significant difference given sample size
    df_res_1st[aggr_names[1],"Pr(>|t|)"]  <- (2*pt(abs(df_res_1st[aggr_names[1],"t value"]), 
                                                   lower.tail = FALSE, 
                                                   df = fixest_df_1st)) 
  }
  if(aggr_dyn & rfs_lead == 0 & rfs_lag == 0 & rfs_fya > 0 & rfs_pya > 0){
    # In this case, we are interested in AVERAGE LEAD AND LAG effects, separately and aggregated
    df_res_1st <- rbind(rep(NA, ncol(df_res_1st)), df_res_1st)
    
    # ORDER MATTERS
    aggr_names <- paste0(original_exposure_rfs, c("_X_aggrall"))
    row.names(df_res_1st)[1] <- aggr_names
    # instruments of interest
    base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
    
    # Contemporaneous, lead, and lag values
    all_roi <- c(grep(pattern = paste0(base_reg_name,"_",rfs_fya,"fya"), 
                      instruments, value = TRUE),
                 grep(pattern = paste0(base_reg_name,"_",rfs_pya,"pya"), 
                      instruments, value = TRUE))
    
    df_res_1st[aggr_names[1],"Estimate"] <- est_1st$coefficients[all_roi] %>% sum()
    
    # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
    # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
    df_res_1st[aggr_names[1],"Std. Error"] <- est_1st$cov.scaled[all_roi, all_roi] %>% as.matrix() %>% sum() %>% sqrt()
    
    df_res_1st[aggr_names[1],"t value"]  <- (df_res_1st[aggr_names[1],"Estimate"] - 0)/(df_res_1st[aggr_names[1],"Std. Error"])
    
    # use t distribution with degrees of freedom equal to that used by fixest, i.e. after two way cluster adjustment.
    # does not make a significant difference given sample size
    df_res_1st[aggr_names[1],"Pr(>|t|)"]  <- (2*pt(abs(df_res_1st[aggr_names[1],"t value"]), 
                                                   lower.tail = FALSE, 
                                                   df = fixest_df_1st)) 
  }
  if(aggr_dyn & rfs_lead > 0 & rfs_lag == 0 & rfs_fya == 0 & rfs_pya > 0){
    # In this case, we are interested in aggregated LEAD effects, and all together.
    
    # first aggregate leads
    df_res_1st <- rbind(rep(NA, ncol(df_res_1st)), df_res_1st)
    
    # ORDER MATTERS
    aggr_names <- paste0(original_exposure_rfs, c("_X_aggrleads"))
    row.names(df_res_1st)[1] <- aggr_names
    # instruments of interest
    base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
    # Contemporaneous value is NOT in leads
    lead_roi <- c(grep(pattern = paste0(base_reg_name,"_lead"), 
                       instruments, value = TRUE))
    
    df_res_1st[aggr_names[1],"Estimate"] <- est_1st$coefficients[lead_roi] %>% sum()
    
    # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
    # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
    df_res_1st[aggr_names[1],"Std. Error"] <- est_1st$cov.scaled[lead_roi, lead_roi] %>% as.matrix() %>% sum() %>% sqrt()
    
    df_res_1st[aggr_names[1],"t value"]  <- (df_res_1st[aggr_names[1],"Estimate"] - 0)/(df_res_1st[aggr_names[1],"Std. Error"])
    
    # use t distribution with degrees of freedom equal to that used by fixest, i.e. after two way cluster adjustment.
    # does not make a significant difference given sample size
    df_res_1st[aggr_names[1],"Pr(>|t|)"]  <- (2*pt(abs(df_res_1st[aggr_names[1],"t value"]), 
                                                   lower.tail = FALSE, 
                                                   df = fixest_df_1st)) 
    
    ## Then (on top of dataframe), aggregate all leads and lags together
    # it's important that overall aggregate comes after, for at least two reasons in current code:
    # 1. because selection of estimate to plot is based on order: it takes the first row of df_res_1st 
    # 2. so that aggr_names corresponds to *_X_aggrall in randomization inference below
    
    df_res_1st <- rbind(rep(NA, ncol(df_res_1st)), df_res_1st)
    
    # ORDER MATTERS
    aggr_names <- paste0(original_exposure_rfs, c("_X_aggrall"))
    row.names(df_res_1st)[1] <- aggr_names
    # instruments of interest
    base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
    
    # Contemporaneous, lead, and lag values
    all_roi <- c(grep(pattern = paste0(base_reg_name,"_lead"), 
                      instruments, value = TRUE),
                 grep(pattern = paste0(base_reg_name,"_",rfs_pya,"pya"), 
                      instruments, value = TRUE))
    
    df_res_1st[aggr_names[1],"Estimate"] <- est_1st$coefficients[all_roi] %>% sum()
    
    # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
    # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
    df_res_1st[aggr_names[1],"Std. Error"] <- est_1st$cov.scaled[all_roi, all_roi] %>% as.matrix() %>% sum() %>% sqrt()
    
    df_res_1st[aggr_names[1],"t value"]  <- (df_res_1st[aggr_names[1],"Estimate"] - 0)/(df_res_1st[aggr_names[1],"Std. Error"])
    
    # use t distribution with degrees of freedom equal to that used by fixest, i.e. after two way cluster adjustment.
    # does not make a significant difference given sample size
    df_res_1st[aggr_names[1],"Pr(>|t|)"]  <- (2*pt(abs(df_res_1st[aggr_names[1],"t value"]), 
                                                   lower.tail = FALSE, 
                                                   df = fixest_df_1st)) 
  }
  
  # output wanted
  
  toreturn <- list(cummulative_annual_effects_1st = df_res_1st[grepl(pattern = exposure_rfs, x = row.names(df_res_1st)),], 
                   wald_jointnull_1st = fitstat(est_1st, "wald"), # F-test not available for GLM. 
                   sum_fv_1st = sum(est_1st$fitted.values), # leave it in bare unit and convert/scale post estimation 
                   APE_estimand = NA # this will be filled below
  )
  
  # # 2nd stage: necessary to compute degrees of freedom, to compute p-value of APE 
  # don't take DoF from 1st stage: there are much more DoF there, as the initial number of observation is higher before removing always null deforestation cells
  d_clean_1st$est_res_1st <- est_1st$residuals
  est_2nd <- fixest::feglm(fml_2nd,
                           data = d_clean_1st,
                           family = distribution_2nd,
                           vcov = se,
                           fixef.rm = "perfect",
                           nthreads = 3,
                           glm.iter = glm_iter,
                           notes = TRUE)
  
  fixest_df_2nd <- degrees_freedom(est_2nd, type = "t")
  
  ## Redefine changes in endogenous variable (first stage outcome) to one standard deviation if asked 
  if(stddev){
    # remove fixed effect variations from the regressor
    reg_sd <- fixest::feols(fml = as.formula(paste0(
      endo_vars,
      " ~ 1 | ", 
      paste0(est_1st$fixef_vars, collapse = "+"))),
      data = d_clean)
    
    # and take the remaining standard deviation
    rel_lu_change <- sd(reg_sd$residuals)
    abs_lu_change <- sd(reg_sd$residuals)
  }
  
  

  # myfun_data <- d_clean
  # fsf <- fml_1st
  # ssf <- fml_2nd
  
  # function applied to each bootstrap replicate
  ctrl_fun_endo <- function(myfun_data, fsf, ssf){
    
    est_1st <- fixest::feglm(fsf, 
                             data = myfun_data, 
                             # vcov = se,
                             family = distribution_1st,
                             fixef.rm = "perfect",
                             nthreads = 3, 
                             notes = TRUE)
    # how the vcov is computed does not change the residuals, so it does not affect the second stage. 
    # However, it affects fit statistics of the first stage, but we look at this outside the bootstrap process
    
    myfun_data_1st <- myfun_data[est_1st$obs_selection[[1]], ]
    
    # save estimated residuals
    myfun_data_1st$est_res_1st <- est_1st$residuals
    
    # 2nd stage
    est_2nd <- fixest::feglm(ssf, 
                             data = myfun_data_1st, 
                             family = distribution_2nd, 
                             fixef.rm = "perfect",
                             nthreads = 3,
                             glm.iter = glm_iter,
                             notes = TRUE)
    # regarding vcov, we don't need to extract SEs from this function, so let the argument be its default here too
    
    ## MAKE APE 
    # identify the nature of different variables
    coeff_names <- names(coef(est_2nd))
    
    k <- which(coeff_names == endo_vars)
    beta_k <- coef(est_2nd)[k]
    
    # the APE is different depending on the regressor of interest being in the log scale or not. 
    if(grepl("ln_",coeff_names[k])){
      ape <- ((1+rel_lu_change)^(beta_k) - 1)*100
    } else{
      ape <- (exp(beta_k*abs_lu_change) - 1)*100 
    } 
    
    # statistics we want to evaluate the variance of:
    return(ape)
}
  
  ## BOOTSTRAP
  
  # helper function
  ran.gen_cluster_blc <- function(original_data, arg_list){
    # to store 
    cl_boot_dat <- list()
    
    # such that we don't have to call it from the list every time
    cluster_var <- arg_list[["cluster_variable"]]
    
    # non-unique names of clusters (repeated when there is more than one obs. in a cluster) 
    nu_cl_names <- as.character(original_data[,cluster_var]) 
    
    # sample, in the vector of names of clusters as many draws as there are clusters, with replacement
    sample_cl <- sample(x = arg_list[["cluster_names"]], 
                        size = arg_list[["number_clusters"]], 
                        replace = TRUE) 
    
    # because of replacement, some names are sampled more than once
    # we need to give them a new cluster identifier, otherwise a cluster sampled more than once 
    # will be "incorrectly treated as one large cluster rather than two distinct clusters" (by the fixed effects) (Cameron and Miller, 2015)    
    sample_cl_tab <- table(sample_cl)
    
    for(n in 1:max(sample_cl_tab)){ # from 1 to the max number of times a cluster was sampled bc of replacement
      # vector to select obs. that are within the clusters sampled n times. 
      # seems slightly faster to construct the names_n object beforehand 
      names_n <- names(sample_cl_tab[sample_cl_tab == n])
      sel <- nu_cl_names %in% names_n
      
      # select data accordingly to the cluster sampling (duplicating n times observations from clusters sampled n times)
      clda <- original_data[sel,][rep(seq_len(sum(sel)), n), ]
      
      #identify row names without periods, and add ".0" 
      row.names(clda)[grep("\\.", row.names(clda), invert = TRUE)] <- paste0(grep("\\.", row.names(clda), invert = TRUE, value = TRUE),".0")
      
      # add the suffix due to the repetition after the existing cluster identifier. 
      clda[,cluster_var] <- paste0(clda[,cluster_var], sub(".*\\.","_",row.names(clda)))
      
      # stack the bootstrap samples iteratively 
      cl_boot_dat[[n]] <- clda
    }
    return(bind_rows(cl_boot_dat))
  }
  # names and numbers of clusters for cluster var 1 
  par_list_dim1 <- list(cluster_variable = cluster_var1, 
                         cluster_names = unique(d_clean[,cluster_var1]),
                         number_clusters = length(unique(d_clean[,cluster_var1])))
  
  ## RUN BOOTSTRAP PROCEDURE
  # store Standard Errors here
  SEs <- c(cluster_var1 = NA, cluster_var2 = NA, cluster_var12 = NA)
  
  # bootstrap on d_clean, no specific subset, because bootstrapping WILL lead to different data sets that have different 
  # patterns of always zero units. 
  bootstraped_1 <- boot(data = d_clean, 
                      statistic = ctrl_fun_endo, # 2 first arguments do not need to be called.
                      # the first one, arbitrarily called "myfun_data" is passed the previous "data" argument 
                      fsf = fml_1st,
                      ssf = fml_2nd,
                      ran.gen = ran.gen_cluster_blc,
                      mle = par_list_dim1,
                      sim = "parametric",
                      R = boot_rep)
  
  SEs[cluster_var1] <- sd(bootstraped_1$t)
  SEfinal <- SEs[cluster_var1]
  
  if(clustering=="twoway"){
    ## Supplementary parameters needed
    # names and numbers of clusters for cluster var 2 
    par_list_dim2 <- list(cluster_variable = cluster_var2, 
                          cluster_names = unique(d_clean[,cluster_var2]),
                          number_clusters = length(unique(d_clean[,cluster_var2])))
    
    # names and numbers of clusters for cluster var 12 
    d_clean <- dplyr::mutate(d_clean, cluster_var12 := paste0(!!as.symbol(cluster_var1), "_", !!as.symbol(cluster_var2)))
    # length(unique(d_clean$cluster_var12)) == nrow(d_clean)
    cluster_var12 <- paste0(cluster_var1, "_", cluster_var2)
    par_list_dim12 <- list(cluster_variable = cluster_var12, 
                           cluster_names = unique(d_clean[,cluster_var12]),
                           number_clusters = length(unique(d_clean[,cluster_var12])))
    
    ## Supplmentary bootstrapping 
    bootstraped_2 <- boot(data = d_clean, 
                        statistic = ctrl_fun_endo, # 2 first arguments do not need to be called.
                        # the first one, arbitrarily called "myfun_data" is passed the previous "data" argument 
                        fsf = fml_1st,
                        ssf = fml_2nd,
                        ran.gen = ran.gen_cluster_blc,
                        mle = par_list_dim2,
                        sim = "parametric",
                        R = boot_rep)
    SEs[cluster_var2] <- sd(bootstraped_2$t)
    
    bootstraped_12 <- boot(data = d_clean, 
                        statistic = ctrl_fun_endo, # 2 first arguments do not need to be called.
                        # the first one, arbitrarily called "myfun_data" is passed the previous "data" argument 
                        fsf = fml_1st,
                        ssf = fml_2nd,
                        ran.gen = ran.gen_cluster_blc,
                        mle = par_list_dim12,
                        sim = "parametric",
                        R = boot_rep)
    SEs[cluster_var12] <- sd(bootstraped_12$t)
    
    ## Final standard errors 
    SEfinal <- sqrt( (SEs[cluster_var1])^2 + (SEs[cluster_var2])^2 - (SEs[cluster_var12])^2 )
    # would need to weight each by G/G-1 maybe ? not sure, and if it was the case, would not change much
  }

  APE_pt_est <- bootstraped_1$t0
  names(APE_pt_est) <- NULL
  names(SEfinal) <- NULL
  APE_estimand <- c(Estimate = APE_pt_est, `Std. Error` = SEfinal, `t value` = (APE_pt_est - 0)/SEfinal, `Pr(>|t|)` = NA)
  
  APE_estimand["Pr(>|t|)"] <- (2*pt(abs(APE_estimand["t value"]), 
                                    lower.tail = FALSE, 
                                    df = fixest_df_2nd)) 
  
  toreturn[["APE_estimand"]] <- APE_estimand
  
  rm(d_clean, d_clean_1st)
  return(toreturn)
  rm(toreturn)
}

