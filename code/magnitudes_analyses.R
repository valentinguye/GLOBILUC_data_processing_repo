

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
  "Sugar", "eaear_Sugar", 
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


### MAIN DATA SET ### 
brazil_data <- readRDS(here("temp_data", "merged_datasets", "brazil_aoi", "driverloss_all_aeay_long_final.Rdata"))
southam_data <- readRDS(here("temp_data", "merged_datasets", "southam_aoi", "driverloss_all_aeay_long_final.Rdata"))

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


### HELPER FUNCTION TO BOOTSTRAP CLUSTER SE 
# will tell boot::boot how to sample data at each replicate of statistic

# without clustering
ran.gen_blc <- function(original_data, arg_list){
  
  new_rowids <- sample(x = row.names(original_data), 
                       size = arg_list[["size"]], 
                       replace = TRUE)
  
  return(original_data[new_rowids,])
}



#### REGRESSION FUNCTION #### 

### TEMPORARY OBJECTS 
outcome_variable = "driven_loss_commodity" # "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
endo_vars = "pastures"
start_year = 2011
end_year = 2019
region = "brazil"

pasture_shares <- FALSE
standardization = ""

estimated_effect = "alpha"
rfs_reg <- TRUE
rfs_rando <- ""
original_rfs_treatments <- c("statute_conv")
rfs_lead <- 3
rfs_lag <- 3 
rfs_fya <-  0
rfs_pya <- 0
aggr_dyn <- TRUE
original_exposure_rfs <- "eaear_Fodder"
group_exposure_rfs <- FALSE
control_all_absolute_rfs <- TRUE
annual_rfs_controls <- FALSE

control_pasture <- FALSE
pasture_trend <- FALSE

fc_trend <- FALSE
s_trend <- TRUE
fc_s_trend <- FALSE

sjpos <- FALSE
fe = "grid_id + country_year" #   
distribution <- "quasipoisson"
invhypsin = TRUE
conley_cutoff <- 100
se = "twoway"
output = "est_object"
glm_iter <- 25 




rm(outcome_variable, start_year, end_year, pasture_shares, 
   standardization, group_prices, biofuel_focus, aggregate_K, control_interact_sj, control_interact_Pk, reference_crop, control_direct, 
   rfs_reg, original_rfs_treatments, rfs_lead, rfs_lag, original_exposure_rfs, group_exposure_rfs, control_absolute_rfs, control_all_absolute_rfs, sjpos, fe, distribution, invhypsin, conley_cutoff, se, boot_cluster, coefs_to_aggregate, 
   output, glm_iter)




make_main_reg <- function(outcome_variable = "driven_loss_commodity", # one of "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
                          endo_vars = "pastures",
                          start_year = 2001, 
                          end_year = 2019, 
                          region = "brazil", # one of "brazil", "southam", or "seasia" 
                          
                          pasture_shares = FALSE, # if TRUE, and crop_j = "Fodder", then qj is proxied with the share of pasture area in 2000. 
                          standardization = "_std2", # one of "", "_std", or "_std2"
                          
                          rfs_rando = "", # either "between", "within", or any other string. If one of the former two, randomization inference of the type is performed
                          original_rfs_treatments = c("statute_conv"),
                          rfs_lead = 0,
                          rfs_lag = 0,
                          rfs_fya = 0, 
                          rfs_pya = 0,
                          aggr_dyn = TRUE, # whether to report aggregate coefficients of all leads and lags ("all") or leads and lags separately ("leadlag"), or no aggregate (any other string)
                          original_exposure_rfs = "eaear_Soy_compo",
                          group_exposure_rfs = FALSE, # if original_exposure_rfs is length 1, it does noe matter whether this option is TRUE or FALSE
                          
                          control_all_absolute_rfs = TRUE,
                          annual_rfs_controls = FALSE,
                          control_pasture = FALSE,
                          pasture_trend = FALSE,
                          # remaining = FALSE, # should remaining forest be controlled for STOP DOING THIS BECAUSE IT INTRODUCES NICKELL BIAS
                          s_trend = TRUE,
                          fc_trend = FALSE,
                          fc_s_trend = FALSE,
                          
                          sjpos = FALSE,
                          
                          fe = "grid_id + country_year", 
                          distribution = "quasipoisson",#  "quasipoisson", 
                          invhypsin = TRUE, # if distribution is gaussian, should the dep. var. be transformed to inverse hyperbolic sine?
                          se = "twoway",# # passed to vcov argument. Currently, one of "cluster", "twoway", or an object of the form: 
                          # - "exposure2ways" for the two-way clustering to be customed as original_exposure_rfs + country_year, and not along FE.   
                          # - vcov_conley(lat = "lat", lon = "lon", cutoff = 100, distance = "spherical")
                          # with cutoff the distance, in km, passed to fixest::vcov_conley, if se = "conley"  
                          glm_iter = 25,
                          output = "coef_table" # one of "data", est_object, or "coef_table" 
){
  
  # and this matrix is used to select exposure type, either eaear or aesi 
  if(grepl("eaear_", original_exposure_rfs) & rfs_reg){
    exp_matmap <- eaear_mapmat
  }else{
    exp_matmap <- mapmat_si
  }
  #### PREPARE NEEDED VARIABLE NAMES
  # this does not involve data, just arguments of the make_reg function
  # original_ names are used to get generic covariate names in coefficient tables (irrespective of modelling choices)

  
  # this is all possible exposures, necessary in every specificaton
  exposure_rfs <- paste0(original_exposure_rfs, standardization)
  
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
  
  # Keep only in data the useful variables 
  d <- dplyr::select(d, all_of(c("grid_id", "year", "lat", "lon", "country_name", "country_year", "continent_name", 
                                 "fc_2000", "fc_2009", "remaining_fc", "pasture_share_2000", 
                                 outcome_variable, endo_vars, exp_matmap[, "Crops"]
                                 ))) 
  
  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", unique(c(all_rfs_treatments, rfs_treatments)))], by = c("year"))#, all_treatments
  
  
  # If we want the exposure to deforestation for pasture to be proxied with the share of pasture area in 2000, rather than suitability index, 
  if(pasture_shares){
    d[,grepl("Fodder", names(d))] <- d$pasture_share_2000
  }
  
  if((distribution == "gaussian" | estimated_effect == "beta") & invhypsin){
    # transform dependent variable, if gaussian GLM 
    d <- dplyr::mutate(d, !!as.symbol(outcome_variable) := asinh(!!as.symbol(outcome_variable)))
  }
 
  
  #### INSTRUMENTS #### 

  if(rfs_reg){
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
  }
  
  ### CONTROLS ####
  # it's important that this is not conditioned on anything so these objects exist
  controls <- c()
  
  # IT IS NORMAL THAT ABS. EXPOSURE CONTROLS ARE INTERACTED WITH 2fya AND 1pya: BOTH ARE AVERAGES OF TWO YEARS. 
  # 1pya is the avg of current and lag1 values, which are grouped together because they reflect the same kind of mechanism. 
  if(control_all_absolute_rfs){
    all_abs <- exp_matmap[, "Crops"][!(exp_matmap[, "Crops"] %in% original_exposure_rfs)]
    
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
  
  
  # from here, if we are in rfs process, sj may refer to the group of exposures
  if((rfs_reg) & group_exposure_rfs){
    sj <- "grp_exp"
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

  used_vars <- unique(c("grid_id", "year", "lat", "lon", "country_name", "country_year", "continent_name", #"remaining_fc", "accu_defo_since2k", # "sj_year",
                       original_exposure_rfs, original_rfs_treatments, rfs_treatments, outcome_variable, endo_vars, instruments, controls))
  
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
  # Though we want identical samples to compare distributional assumptions)
  if(distribution == "gaussian"){
    # # - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove
    obstormv <- obs2remove(fml = as.formula(paste0(outcome_variable, " ~ ", fe)),
                           data = d,
                           family = "poisson")
    # this is necessary to handle cases when there is no obs. to remove
    if(length(obstormv)>0){
      d <- d[-obstormv,]
    } else {
      d <- d
    }
    rm(obstormv)
  }
  
  if(fe != "grid_id + country_year"){
    # # - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove
    obstormv <- obs2remove(fml = as.formula(paste0(outcome_variable, " ~ grid_id + country_year")),
                           data = d,
                           family = "poisson")
    # this is necessary to handle cases when there is no obs. to remove
    if(length(obstormv)>0){
      d <- d[-obstormv,]
    } else {
      d <- d
    }
    rm(obstormv)
  }
  
  # if we don't remove obs that have no remaining forest, and if we don't remove always 0 outcome obs. in previous step, d IS A BALANCED PANEL !!
  # this is important for the cluster boostrap approach that follows
  if(!nrow(d) == length(unique(d$grid_id))*length(unique(d$year))){
    warning("data is not balanced")
  }
  d_clean <- d
  rm(d)
  
  ### FORMULAE ####
  
  # endo var should be cropland extent in t+l (lead) 
  fml_1st <- as.formula(paste0(endo_vars,
                               " ~ ", 
                               paste0(instruments, collapse = "+"),
                               " + ",
                               paste0(controls, collapse = "+"),
                               " | ", 
                               fe))
  
  fml_2nd <- as.formula(paste0(outcome_variable,
                               " ~ ", endo_vars,
                               " + est_res_1st +", 
                               paste0(controls, collapse = "+"),
                               " | ",
                               fe))
  
  ### REGRESSIONS
  # myfun_data <- d_clean
  # fsf <- fml_1st
  # ssf <- fml_2nd
  
  ctrl_fun_endo <- function(myfun_data, fsf, ssf){
    
    # 1st stage 
    est_1st <- fixest::feols(fsf, 
                             data = myfun_data, 
                             fixef.rm = "perfect",
                             nthreads = 3, 
                             notes = TRUE)
    # how the vcov is computed does not change the residuals
    
    # save estimated residuals
    myfun_data$est_res_1st <- est_1st$residuals
    # 2nd stage
    est_2nd <- fixest::feglm(ssf, 
                             data = myfun_data, 
                             family = distribution, 
                             fixef.rm = "perfect",
                             nthreads = 3,
                             glm.iter = glm_iter,
                             notes = TRUE)
    
    # statistics we want to evaluate the variance of:
    return(est_2nd$coefficients[endo_vars])
  }
  
  ## bootstrap with 2-way clustering

  # names and numbers of clusters for cluster var 1 
  cluster_var1 <- "grid_id"
  par_list_dim1 <- list(cluster_variable = cluster_var1, 
                         cluster_names = unique(d_clean[,cluster_var1]),
                         number_clusters = length(unique(d_clean[,cluster_var1])))
  
  # names and numbers of clusters for cluster var 2 
  cluster_var2 <- "year"
  par_list_dim2 <- list(cluster_variable = cluster_var2, 
                        cluster_names = unique(d_clean[,cluster_var2]),
                        number_clusters = length(unique(d_clean[,cluster_var2])))
  
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
  
  # test_boot_d <- ran.gen_cluster_blc(original_data = d_clean,
  #                                    arg_list = par_list_dim1)
  # 
  # dim(test_boot_d)
  # dim(d_clean)
  
  # # test new clusters are not duplicated (correct if anyDuplicated returns 0)
  # base::anyDuplicated(test_boot_d[,c(cluster_var1,cluster_var2)])
  
  # Compute inference statistics by bootstrap
  bootstraped_1 <- boot(data = d_clean, 
                      statistic = ctrl_fun_endo, # 2 first arguments do not need to be called.
                      # the first one, arbitrarily called "myfun_data" is passed the previous "data" argument 
                      fsf = fml_1st,
                      ssf = fml_2nd,
                      ran.gen = ran.gen_cluster_blc,
                      mle = par_list_dim1,
                      sim = "parametric",
                      R = 2)
  
  bootstraped_2 <- boot(data = d_clean, 
                      statistic = ctrl_fun_endo, # 2 first arguments do not need to be called.
                      # the first one, arbitrarily called "myfun_data" is passed the previous "data" argument 
                      fsf = fml_1st,
                      ssf = fml_2nd,
                      ran.gen = ran.gen_cluster_blc,
                      mle = par_list_dim2,
                      sim = "parametric",
                      R = 2)
  
  bootstraped_12 <- boot(data = d_clean, 
                      statistic = ctrl_fun_endo, # 2 first arguments do not need to be called.
                      # the first one, arbitrarily called "myfun_data" is passed the previous "data" argument 
                      fsf = fml_1st,
                      ssf = fml_2nd,
                      ran.gen = ran.gen_cluster_blc,
                      mle = par_list_dim12,
                      sim = "parametric",
                      R = 2)
  norm.ci(bootstraped, index = 1) # equivalent to : 
  
  # Store only information necessary, in a dataframe. otherwise the output of fixest estimation is large and we can't collect too many at the same time (over loops)  
  # either there are several elemnts in instruments, and then we want to aggregate them, or there is only one. 
  # In both cases, we are interested in a one-line output
  
  # handle SE computation flexibly within feglm now, through argument vcov
  
  # in case of unit fe, adjust also the vcov argument, because two way wont work, but we might still want to cluster on country year
  if(fe == "grid_id" & se == "twoway"){
    se <- as.formula("~ grid_id + country_year")
  }
  if(fe == "country_year" & se == "twoway"){
    se <- as.formula("~ grid_id + country_year")
  }
  if(se == "exposure2ways"){
    se <- as.formula(paste0("~ ", original_exposure_rfs, " + country_year"))
  }
  
  if(estimated_effect == "alpha"){
    alpha_reg_res <- fixest::feglm(alpha_model,
                                   data = d_clean, 
                                   family = distribution,# "gaussian",#  # "poisson" ,
                                   vcov = se,
                                   # this is just to get the same p value by recomputing by hand below. 
                                   # see https://cran.r-project.org/web/packages/fixest/vignettes/standard_errors.html
                                   # ssc = ssc(cluster.df = "conventional", t.df = "conventional"),
                                   # glm.iter = 25,
                                   #fixef.iter = 100000,
                                   nthreads = 3,
                                   fixef.rm = "perfect",
                                   glm.iter = glm_iter,
                                   notes = TRUE, 
                                   verbose = 4) 
  }    
    
  df_res <- summary(alpha_reg_res)$coeftable
  
  fixest_df <- degrees_freedom(alpha_reg_res, type = "t")
    
  # store info on SE method
  if(se == "conley"){
    se_info <- paste0("Conley (",conley_cutoff,"km)")
  }
  if(se == "twoway"){
    # " If the two variables were used as fixed-effects in the estimation, 
    # you can leave it blank with vcov = "twoway""
    se_info <- "Two-way clustered (grid cell - year)"
  }
  if(se == "cluster"){
    se_info <- paste0("Clustered (",boot_cluster,")")
  }
  if(estimated_effect != "alpha"){
    se_info <- paste0("Clustered (",boot_cluster,")")
  }
  
  # output wanted
  if(output == "est_object"){
    toreturn <- alpha_reg_res
  }
  
  if(output == "data"){
    toreturn <- list(df_res, d_clean)
  }
  # if(output == "est_obj"){
  #   toreturn <- sum_res
  # }
  if(output == "coef_table"){
    toreturn <- df_res
  }
  
  if(rfs_rando=="between" | rfs_rando=="within"){
    toreturn <- list(pvals_tstats = raninf_pval_tstats, 
                     pvals_tstats_possided = raninf_possided_pval_tstats,
                     pvals_tstats_negsided = raninf_negsided_pval_tstats,
                     pvals_coeffs = raninf_pval_coeffs, 
                     all_coeffs = all_coeffs, 
                     all_tstats = all_tstats)
  }
  
  
  rm(d_clean, df_res)
  return(toreturn)
  rm(toreturn)
}

