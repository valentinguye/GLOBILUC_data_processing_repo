

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
                   "ggplot2", "dotwhisker", #"tmap",# "leaflet", "htmltools"
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


### This matrix is used for maping crops from GAEZ with commodities from price data sets
mapmat_data <- c(
  "Banana","Banana",
  "Barley", "Barley",
  "Beef", "Fodder", 
  "Orange", "Citrus", # Citrus sinensis in both GAEZ and FAO
  "Cocoa", "Cocoa",
  "Coconut_oil", "Coconut",
  "Coffee", "Coffee",
  "Cotton", "Cotton",
  "Groundnut", "Groundnut",
  "Maize", "Maizegrain",
  "Oat", "Oat",
  "Olive_oil", "Olive",
  "Palm_oil", "Oilpalm",
  "Rapeseed_oil", "Rapeseed",
  "Rice", "Rice",
  "Rubber", "Rubber",
  "Sorghum", "Sorghum2",
  "Soybean", "Soybean",
  # "Soybean_oil", "Soybean",
  # "Soybean_meal", "Soybean",
  "Sugar", "Sugar",
  # "Sugarbeet", "Sugarbeet",
  # "Sugarcane", "Sugarcane",
  "Sunflower_oil", "Sunflower",
  "Tea", "Tea",
  "Tobacco", "Tobacco",
  "Wheat", "Wheat")

mapmat <- matrix(data = mapmat_data, 
                 nrow = length(mapmat_data)/2,
                 ncol = 2, 
                 byrow = TRUE)

colnames(mapmat) <- c("Prices", "Crops")

# there are different ways of grouping crops in a relevant way, depending on the exposure crop, sj
meats = c("Beef", "Chicken", "Pork", "Sheep")
vegetable_oils = c("Palm_oil", "Rapeseed_oil", "Soybean_oil", "Sunflower_oil")
cereals_feeds = c("Barley", "Maize", "Oat", "Rice", "Soybean", "Soybean_meal", "Wheat")
biofuel_feedstocks = c("Maize", "Palm_oil", "Rapeseed_oil", "Sorghum", "Soybean_oil", "Sugar")

maplist <- list(Fodder = list(meats = meats, cereals_feeds = cereals_feeds, vegetable_oils = vegetable_oils, 
                              Crude_oil = "Crude_oil"), 
                Soybean = list(meats = meats, cereals_feeds = cereals_feeds, vegetable_oils = vegetable_oils, 
                               Crude_oil = "Crude_oil", Sorghum = "Sorghum", Sugar = "Sugar", Coconut_oil = "Coconut_oil"), 
                Oilpalm = list(meats = meats, cereals_feeds = cereals_feeds, vegetable_oils = vegetable_oils, 
                               Crude_oil = "Crude_oil", Sorghum = "Sorghum", Sugar = "Sugar", Coconut_oil = "Coconut_oil"), 
                Cocoa = list(Cocoa = "Cocoa", Coffee = "Coffee", vegetable_oils = vegetable_oils, Sugar = "Sugar", Crude_oil = "Crude_oil"), 
                Coffee = list(Coffee = "Coffee", Cocoa = "Cocoa", Tea = "Tea", Tobacco = "Tobacco", Crude_oil = "Crude_oil"))

### MAIN DATA SET ### 
main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "driverloss_aesi_long_final.Rdata"))

### PRICE DATA ###
prices <- readRDS(here("temp_data", "prepared_international_prices.Rdata"))

### HELPER FUNCTION TO BOOTSTRAP CLUSTER SE 
# will tell boot::boot how to sample data at each replicate of statistic
ran.gen_cluster_blc <- function(original_data, arg_list){
  # start <- Sys.time()
  # to store 
  cl_boot_dat <- list()
  
  cluster_var <- arg_list[["cluster_variable"]]
  
  # non-unique names of clusters (repeated when there is more than one obs. in a cluster) 
  nu_cl_names <- as.character(original_data[,cluster_var]) 
  
  # sample, in the vector of names of clusters as many draws as there are clusters, with replacement
  sample_cl <- sample(arg_list[["cluster_names"]], 
                      arg_list[["number_clusters"]], 
                      replace = TRUE) 
  
  # because of replacement, some names are sampled more than once
  # we need to give them a new cluster identifier, otherwise a cluster sampled more than once 
  # will be "incorrectly treated as one large cluster rather than two distinct clusters" (by the fixed effects) (Cameron and Miller, 2015)    
  sample_cl_tab <- table(sample_cl)
  
  for(n in 1:max(sample_cl_tab)){ # from 1 to the max number of times a cluster was sampled bc of replacement
    # vector to select obs. that are within the clusters sampled n times. 
    # seems slightly faster to construct the names_n object beforehand 
    # start <- Sys.time()
    names_n <- names(sample_cl_tab[sample_cl_tab == n])
    sel <- nu_cl_names %in% names_n
    # end <- Sys.time()
    
    # select data accordingly to the cluster sampling (duplicating n times observations from clusters sampled n times)
    clda <- original_data[sel,][rep(seq_len(sum(sel)), n), ]
    
    #identify row names without periods, and add ".0" 
    row.names(clda)[grep("\\.", row.names(clda), invert = TRUE)] <- paste0(grep("\\.", row.names(clda), invert = TRUE, value = TRUE),".0")
    
    # add the suffix due to the repetition after the existing cluster identifier. 
    clda[,cluster_var] <- paste0(clda[,cluster_var], sub(".*\\.","_",row.names(clda)))
    
    # stack the bootstrap samples iteratively 
    cl_boot_dat[[n]] <- clda
  }
  # end <- Sys.time()
  return(bind_rows(cl_boot_dat))
}




### TEMPORARY OBJECTS 
outcome_variable = "driven_loss" # "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
start_year = 2001
end_year = 2019
continent = "all"
further_lu_evidence = "none"
original_sj = "Fodder"# ,"Fodder",  "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber" in GAEZ spelling, one, part, or all of the 6 main drivers of deforestation: 
# for the k variables hypothesized in overleaf for palm oil, feglm quasipoisson converges within 25 iter.
# original_treatments <-  c("meats", "vegetable_oils", "cereals_feed", "biofuel_feedstocks") # in price spelling. One, part, or all of the full set of commodities having a price-AESI match.
# "Beef" , , , "Cocoa", "Coffee", "Rubber"
# "Rapeseed_oil", "Sunflower_oil","Rice", "Wheat", "Maize", "Sugar", "Sorghum")
# 
focal_j_extra_price = "Soybean"
extra_price_k = c() # in price spelling. One, part, or all of the full set of commodities NOT having a price-AESI match.
# ,"Chicken", "Pork", "Crude_oil

pasture_shares <- FALSE
j_soy = "Soybean"
fcr = 7.2
standardization = "_std2"
price_info = "_lag1"
estimated_effect = "alpha"
group_prices <- TRUE
control_interact_sj <- FALSE
control_interact_Pk <- FALSE
reference_crop = c("Oat", "Olive")#
control_direct = FALSE
# sjpj_lag = "_lag1" # either "" or "_lag1" or "_lag2"
# skpk_lag = "_lag1" # either "" or "_lag1" or "_lag2"SjPj = TRUE
# SjPj = TRUE
# SkPk = TRUE
remaining <- TRUE
sjpos <- TRUE # should the sample be restricted to cells where sj is positive? 
# open_path <- FALSE
# commoXcommo <- "Fodder_X_Beef"
fe = "grid_id + country_year" # + country_year
distribution = "gaussian"
invhypsin = TRUE
conley_cutoff <- 100
se = "twoway"
boot_cluster ="grid_id"
#coefstat = "confint"
coefs_to_aggregate <- c("Fodder",  "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber")
output = "coef_table"
glm_iter <- 25 



make_main_reg <- function(outcome_variable = "driven_loss", # one of "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
                          start_year = 2001, 
                          end_year = 2020, 
                          continent = "all", # one of "Africa", "America", "Asia", or "all"
                          further_lu_evidence = "none", # either "none", "sbqt_direct_lu", or "sbqt"mode_lu"
                          original_sj = c("Fodder"),  # in GAEZ spelling, one, part, or all of the 6 main drivers of deforestation: 
                          # , "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber"
                          # note that this is the vector that defines what is removed from full_control
                          original_treatments = c("Soybean"), # in price spelling. One, part, or all of the full set of commodities having a price-AESI match.
                          # , "Palm_oil", "Cocoa", "Coffee", "Rubber", 
                          # "Rapeseed_oil", "Sunflower_oil","Rice", "Wheat", "Maize", "Sugar", "Sorghum"
                          focal_j_extra_price = "Soybean", # in GAEZ spelling, the crop of interest if the treatment is price that has no match in GAEZ (an "extra price" in previous saying)
                          # extra_price_k = c(), # "Chicken", "Pork", "Sheep", "Crude_oil" in price spelling. One, part, or all of the full set of commodities NOT having a price-AESI match.
                          pasture_shares = FALSE, # if TRUE, and crop_j = "Fodder", then qj is proxied with the share of pasture area in 2000. 
                          standardization = "_std2", # one of "", "_std", or "_std2"
                          price_info = "_lag1", # one of "lag1", "_2pya", "_3pya", "_4pya", "_5pya",
                          estimated_effect = "alpha",# if "alpha", estimates the (aggregated or not) cross-elasticity effect. If different from "alpha" (e.g. "gamma") it estimates the ILUC effect from price covariation. 
                          group_prices = TRUE,
                          control_interact_sj = FALSE, 
                          control_interact_Pk = FALSE,
                          reference_crop = "Olive",
                          control_direct = FALSE,
                          # SjPj = FALSE,
                          # SkPk = TRUE,
                          # sjpj_lag = "_lag1", # either "" or "_lag1" or "_lag2"
                          # skpk_lag = "_lag1", # either "" or "_lag1" or "_lag2"
                          sjpos = TRUE, # should the sample be restricted to cells where sj is positive? 
                          remaining = TRUE, # should remaining forest be controlled for 
                          # open_path = FALSE,
                          # commoXcommo = "Fodder_X_Beef",
                          #commo_m = c(""), comment coder ça pour compatibilité avec loops over K_commo ? 
                          fe = "grid_id + country_year", 
                          distribution = "gaussian",#  "quasipoisson", 
                          invhypsin = TRUE, # if distribution is gaussian, should the dep. var. be transformed to inverse hyperbolic sine?
                          se = "twoway", # passed to vcov argument in fixest::summary. Currently, one of "cluster", "twoway", or an object of the form:    
                          # vcov_conley(lat = "lat", lon = "lon", cutoff = 100, distance = "spherical")
                          # with cutoff the distance, in km, passed to fixest::vcov_conley, if se = "conley"  
                          boot_cluster ="grid_id",
                          # old argument: cluster ="grid_id", # the cluster level if se = "cluster" (i.e. one way)
                          # coefstat = "confint", # one of "se", "tstat", "confint"
                          glm_iter = 25,
                          coefs_to_aggregate = c("Fodder",  "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber"),
                          output = "coef_table" # one of "data", or "coef_table" 
){
  
  
  #### PREPARE NEEDED VARIABLE NAMES
  # this does not involve data, just arguments of the make_reg function
  # original_ names are used to get generic covariate names in coefficient tables (irrespective of modelling choices)
  
  ## Standardized suitability index to find in the main data 
  # this just those for which focus is set in current specification
  sj <- paste0(original_sj, standardization)  
  # this is all possible exposures, necessary in every specificaton 
  original_all_exposures <- mapmat[,"Crops"]
  
  # (maybe for alpha only currently) identify the L - 1 exposures that will be interacted with the treatment.
  # this is necessary to avoid perfect colinearity with time FE
  # if(original_treatments %in% mapmat[,"Prices"]){
  #   original_all_exposures_but_k <- original_all_exposures[original_all_exposures != mapmat[mapmat[,"Prices"]==original_treatments, "Crops"]]
  # }else{
  #   # handles cases when treatments is an extra price with no match in GAEZ 
  #   original_all_exposures_but_k <- original_all_exposures[original_all_exposures != focal_j_extra_price]
  # }
  
  # standardize them
  all_exposures <- paste0(original_all_exposures, standardization)
  # all_exposures_but_k <- paste0(original_all_exposures_but_k, standardization)
  
  ## Price variable names to find in the price data
  original_treatments <- unlist(maplist[[original_sj]])
  
  # this just those for which focus is set in current specification
  treatments <- paste0(original_treatments, price_info)
  # this is all possible treatments, necessary in every specificaton 
  original_all_treatments <- mapmat[,"Prices"]
  all_treatments <- paste0(original_all_treatments, price_info)
  
  #### MAKE THE VARIABLES NEEDED IN THE DATA
  # manipulate a different data set so that original one can be provided to all functions and not read again every time. 
  d <- main_data
  
  # Keep only in data the useful variables 
  d <- dplyr::select(d, all_of(c("grid_id", "year", "lat", "lon", "country_name", "country_year", "continent_name", "remaining_fc", "accu_defo_since2k",
                                 outcome_variable, "pasture_share",
                                 unique(c(sj, all_exposures))))) #
  
  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", unique(c(treatments, all_treatments)))], by = c("year"))
  
  
  # If we want the exposure to deforestation for pasture to be proxied with the share of pasture area in 2000, rather than suitability index, 
  # if(pasture_shares){
  #   d[,grepl("Fodder", names(d))] <- d$pasture_share
  # }
  
  # transform dependent variable, if gaussian GLM
  if(distribution == "gaussian" & invhypsin){
    d <- dplyr::mutate(d, !!as.symbol(outcome_variable) := asinh(!!as.symbol(outcome_variable)))
  }
  
  
  if(group_prices){
    ## make groups specific to the present sj
    # remove Pj from crops to group
    original_Pj <- mapmat[mapmat[,"Crops"]==original_sj, "Prices"]
    Pj <- paste0(original_Pj, price_info)
    original_all_treatments_but_j <- original_treatments[original_treatments != original_Pj]
    all_treatments_but_j <- paste0(original_all_treatments_but_j, price_info)
    
    sj_group_list <- list()
    for(sj_grp in names(maplist[[original_sj]])){
      sj_group_list[[sj_grp]] <- all_treatments_but_j[original_all_treatments_but_j %in% maplist[[original_sj]][[sj_grp]]]
    }
    
    # sj_group_list <- list(meats = all_treatments_but_j[original_all_treatments_but_j %in% maplist[[original_sj]][["meats"]]], 
    #                       cereals_feeds = all_treatments_but_j[original_all_treatments_but_j %in% maplist[[original_sj]][["cereals_feeds"]]], 
    #                       vegetable_oils = all_treatments_but_j[original_all_treatments_but_j %in% maplist[[original_sj]][["vegetable_oils"]]], 
    #                       biofuel_feedstocks = all_treatments_but_j[original_all_treatments_but_j %in% maplist[[original_sj]][["biofuel_feedstocks"]]],
                            
    
    used_group_names <- names(sj_group_list[lengths(sj_group_list)>0]) # currently useless
    for(grp in used_group_names){
      # log the prices before grouping them
      for(commo in sj_group_list[[grp]]){
        d <- mutate(d, 
                    !!as.symbol(commo) := log(!!as.symbol(commo)))
      }
      # average the logarithms of prices (preserving the price info in the name, and thus replacing non-grouped crops)
      d <- mutate(d, 
                  !!as.symbol(paste0(grp,price_info)) := rowMeans(across(.cols = (any_of(sj_group_list[[grp]])))))
    }
    
    # log Pj 
    d <- mutate(d, !!as.symbol(Pj) := log(!!as.symbol(Pj)))
    
    original_grouped_treatments <- c(original_Pj, used_group_names)
    grouped_treatments <- paste0(original_grouped_treatments, price_info)
    
    regressors <- c()
    for(Pk in grouped_treatments){
      varname <- paste0(original_sj, "_X_", original_grouped_treatments[match(Pk, grouped_treatments)])
      regressors <- c(regressors, varname)
      # Log prices so their variations are comparable
      d <- mutate(d,
                  !!as.symbol(varname) := (!!as.symbol(sj)) * (!!as.symbol(Pk)) )# don't loge Pk again, it as been pre logged above. 
      
    }
    
  }else{
    
    ### Main regressors
    # code if we consider there are several regressors, not coded like this currently
    regressors <- c()
    for(Pk in treatments){
      varname <- paste0(original_sj, "_X_", original_treatments[match(Pk, treatments)])
      regressors <- c(regressors, varname)
      # Log prices so their variations are comparable
      d <- mutate(d,
                  !!as.symbol(varname) := (!!as.symbol(sj)) * log( (!!as.symbol(Pk)) ))#+1
  
    }
    rm(varname, Pk)
    
    # regressors <- paste0(original_sj, "_X_", original_treatments)
    # d <- mutate(d, !!as.symbol(regressors) := (!!as.symbol(sj)) * log( (!!as.symbol(treatments)) ))#+1
  }  
  
  ### Controls - mechanisms
  # it's important that this is not conditioned on anything so these objects exist
  alpha_controls <- c()
  gamma_controls <- c() 
  delta_controls <- c() 
  
  # add remainging forest cover as a control
  if(remaining){
    alpha_controls <- c(alpha_controls, "remaining_fc")
    gamma_controls <- c(gamma_controls, "remaining_fc")
    delta_controls <- c(delta_controls, "remaining_fc")
  }
  
    # if we are to group these controls in one single variable
  if(control_interact_sj){
    sj_interactions <- c()
    for(Pk in all_treatments[all_treatments!=treatments]){
      varname <- paste0(original_sj, "_X_", original_all_treatments[match(Pk, all_treatments)])
      sj_interactions <- c(sj_interactions, varname)
      # Log prices so their variations are comparable
      d <- mutate(d,
                  !!as.symbol(varname) := (!!as.symbol(sj)) * log( (!!as.symbol(Pk)) ))#+1
      }
    rm(varname, Pk)
    
    d <- mutate(d, one_ctrl = rowSums(across(.cols = (any_of(sj_interactions)))))
    alpha_controls <- c(alpha_controls, "one_ctrl")
  }
  
  
  # add indirect interactions with sj
  if(control_interact_Pk){

    # code like this allows to control which controls to drop
    s_ref <- paste0(reference_crop, standardization)
    # handle cases where s_ref has several elements
    d <- dplyr::mutate(d, s_ref = rowSums(across(.cols = (any_of(s_ref)))))
    
    # for(si in all_exposures[!(all_exposures %in% c(sj, s_ref))]){
    #   varname <- paste0(original_all_exposures[match(si, all_exposures)], "_X_", original_treatments)
    #   alpha_controls <- c(alpha_controls, varname)
    #   # Log prices so their variations are comparable
    #   d <- mutate(d,
    #               !!as.symbol(varname) := (!!as.symbol(si)) * log( (!!as.symbol(treatments)) ))#+1
    # }
    
    # for(si in all_exposures_but_k[all_exposures_but_k != sj]){ # do not put the regressor a second time
    #   varname <- paste0(original_all_exposures_but_k[match(si, all_exposures_but_k)], "_X_", original_treatments)
    #   alpha_controls <- c(alpha_controls, varname)
    #   # Log prices so their variations are comparable
    #   d <- mutate(d,
    #               !!as.symbol(varname) := (!!as.symbol(si)) * log( (!!as.symbol(treatments)) ))#+1
    # }
    # rm(varname, si)
    
    # Code below to control for other indirect interactions with Pk through controlling only for (1-sj-s_ref)Pk  
    # (which is supposed to be equivalent to main approach)

    varname <- "one_ctrl"
    alpha_controls <- c(alpha_controls, varname)
    # Log prices so their variations are comparable
    d <- mutate(d,
                !!as.symbol(varname) := (1 - !!as.symbol(sj) - s_ref ) * log( !!as.symbol(treatments) ) )#
    rm(varname)
    
    # Code below to control for other indirect interactions with Pk through controlling only for skPk 
    # (which is bad, because it is equivalent to (1-sk)Pk and thus throws sjPk once again)
    # original_sk <- mapmat[mapmat[,"Prices"]==original_treatments, "Crops"]
    # sk <- paste0(original_sk, standardization)
    # varname <- paste0(original_sk, "_X_", original_treatments)
    # alpha_controls <- c(alpha_controls, varname)
    # # Log prices so their variations are comparable
    # d <- mutate(d,
    #             !!as.symbol(varname) := (!!as.symbol(sk)) * log( (!!as.symbol(treatments)) ))#+1
    # rm(original_sk, sk, varname)
    
    
  }
  
  # add direct interactions
  if(control_direct){
    for(Pm in all_treatments[all_treatments != treatments]){ # do not put the direct interaction with the treatment of interest, otherwise perfect colinearity
      # find the suitability index that matches Pm
      original_Pm <- original_all_treatments[match(Pm, all_treatments)]
      original_sm <- mapmat[mapmat[,"Prices"]==original_Pm, "Crops"]
      sm <- paste0(original_sm, standardization)
      varname <- paste0(original_sm, "_X_", original_Pm)
      alpha_controls <- c(alpha_controls, varname)
      # Log prices so their variations are comparable
      d <- mutate(d,
                  !!as.symbol(varname) := (!!as.symbol(sm)) * log( (!!as.symbol(Pm)) ))#+1
    }
    rm(varname, original_Pm, original_sm, sm, Pm)
  }
  
  ## MODEL SPECIFICATION FORMULAE - alpha model
  if(length(alpha_controls) > 0){ 
    alpha_model <- as.formula(paste0(outcome_variable,
                                     " ~ ",
                                     paste0(regressors, collapse = "+"),
                                     " + ",
                                     paste0(alpha_controls, collapse = "+"),
                                     " | ",
                                     fe))
  }else{
    alpha_model <- as.formula(paste0(outcome_variable,
                                     " ~ ",
                                     paste0(regressors, collapse = "+"),
                                     " | ",
                                     fe))
  }
  
  # it is not needed to do all that if we only estimate alphas 
  if(estimated_effect != "alpha"){
    # construct all interaction variables
    indiv_controls <- c()
    d$full_control <- 0
    for(Pk in all_treatments){
      for(si in all_exposures){
        varname <- paste0(original_all_exposures[match(si, all_exposures)], "_X_", original_all_treatments[match(Pk, all_treatments)])
        # indiv_controls <- c(indiv_controls, varname)
        # d <- mutate(d, 
        #             !!as.symbol(varname) := log( (!!as.symbol(si)) * (!!as.symbol(Pk)) +1))
        # note that mutate will overwrite the already existing varnames created in the regressors step above. 
        
        # we need to add up the new terms together as we create them, otherwise there are 529 new columns added to the 4M+ rows dataframe and this is memory intensive
        if(!(varname%in%regressors)){ # add the interaction term if it's not a regressor already
          d <- mutate(d, full_control := ( full_control + log( (!!as.symbol(si)) * (!!as.symbol(Pk)) +1 ) ) )
        }
      }
    }
    
    # c(log(d[1,"Banana_lag1"]*d[1,"Barley_std2"] +1), log(d[1,"Barley_lag1"]*d[1,"Barley_std2"] +1))%>%sum()
    
    rm(varname, Pk, si)
    
    delta_controls <- c(delta_controls, "full_control")
    
    ## Partial control (to estimate gamma)
    # identify the variables to put in the partial control term, but this changes depending on the question we are answering.
    # yet, it is not necessary to condition: the loop handles it, 
    # bc in all cases, the removed controls capture the confounding covariation between the treatment and the prices of the main crops directly driving deforestation 
    d$part_control <- d$full_control
    for(si in sj){
      # from the inputed si, find the corresponding price to match it to, construct the SjPj interaction, and remove it from the control. 
      original_si <- original_sj[match(si,sj)]
      original_Pj <- mapmat[mapmat[,"Crops"]==original_Sj,"Prices"]
      Pj <- all_treatments[match(original_Pj,original_all_treatments)]
      # it's important that part_control be on both sides of the equal, and that it starts from full_control, to handle cases when original_sj has several elements. 
      d <- mutate(d, part_control := ( part_control - (log( (!!as.symbol(si)) * (!!as.symbol(Pj)) +1) ) ) ) # (note that the regressors are already taken out of full_control)
    }
    
    rm(original_Sj, original_Pj, si, Pj)
    gamma_controls <- c(gamma_controls, "part_control")
    
    
    ## MODEL SPECIFICATION FORMULAE - gamma and delta models
    gamma_model <- as.formula(paste0(outcome_variable,
                                     " ~ ",
                                     paste0(regressors, collapse = "+"),
                                     " + ",
                                     paste0(gamma_controls, collapse = "+"),
                                     " | ",
                                     fe)) 
    
    delta_model <- as.formula(paste0(outcome_variable,
                                     " ~ ",
                                     paste0(regressors, collapse = "+"),
                                     " + ",
                                     paste0(delta_controls, collapse = "+"),
                                     " | ",
                                     fe)) 
  }
  
  
  ### KEEP OBSERVATIONS THAT: 
  
  # - are suitable to crop j 
  if(sjpos){
    d <- dplyr::filter(d, !!as.symbol(sj) > 0)  
  }
  
  # - are in study period 
  d <- dplyr::filter(d, year >= start_year)
  d <- dplyr::filter(d, year <= end_year)
  
  # - are in study area
  if(continent != "all"){
    d <- dplyr::filter(d, continent_name == continent)
  }
  
  # have remaining forest
  # d <- dplyr::filter(d, remaining_fc > 0)
  
  if(estimated_effect == "alpha"){
    used_vars <- unique(c("grid_id", "year", "lat", "lon", "country_name", "country_year", "continent_name",
                          outcome_variable, regressors, alpha_controls))
  }else{
    used_vars <- unique(c("grid_id", "year", "lat", "lon", "country_name", "country_year", "continent_name",
                          outcome_variable, regressors, gamma_controls, delta_controls))
  }
  
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
  
  # stop using obstorm as it is deprecated in the present version of fixest, and did weird things 
  # removing all obs with outcome_variable == 0 (not only when it's always the case within the fe dimension)
  # # - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove
  # obstormv <- obs2remove(fml = as.formula(paste0(outcome_variable, " ~ ", fe)),
  #                        d, 
  #                        family = "poisson")
  # # this is necessary to handle cases when there is no obs. to remove
  # if(length(obstormv)>0){
  #   d_clean <- d[-obstormv,]
  # } else {
  #   d_clean <- d
  # }
  # rm(d, obstormv)
  
  # AND because we don't remove obs that have no remaining forest, d IS A BALANCED PANEL !!
  if(!nrow(d) == length(unique(d$grid_id))*length(unique(d$year))){
    stop("data is not balanced")
  }
  d_clean <- d
  rm(d)
  
  ### REGRESSIONS
  
  # Store only information necessary, in a dataframe. otherwise the output of fixest estimation is large and we can't collect too many at the same time (over loops)  
  # either there are several elemnts in regressors, and then we want to aggregate them, or there is only one. 
  # In both cases, we are interested in a one-line output
  df_res <- data.frame(estimate = NA, std.error = NA, t.statistic = NA, p.value = NA, observations = NA, inference = NA)  
  
  # handle SE computation flexibly within feglm now, through argument vcov
  if(estimated_effect == "alpha"){
    alpha_reg_res <- fixest::feglm(alpha_model,
                                   data = d_clean, 
                                   family = distribution,# "gaussian",#  # "poisson" ,
                                   vcov = se,
                                   # this is just to get the same p value by recomputing by hand below. 
                                   # see https://cran.r-project.org/web/packages/fixest/vignettes/standard_errors.html
                                   ssc = ssc(cluster.df = "conventional", t.df = "conventional"),
                                   # glm.iter = 25,
                                   #fixef.iter = 100000,
                                   nthreads = 3,
                                   glm.iter = glm_iter,
                                   notes = TRUE, 
                                   verbose = 4)  
    
    # keep only the coeff estimate for the sjPk term of interest
    df_res <- summary(alpha_reg_res)$coeftable[regressors, ]
    
    ## MAKE AGGREGATE RESULTS 
    #if(aggr_J){
    # select the coefficients to be addep up (code works in cases when there is only one regressor because we investigate specific crop-crop interactions)
    # (I have checked it)
    # df_res <- rbind(df_res, rep(NA, ncol(df_res)))
    # 
    # row.names(df_res) <- c(rownames(summary(alpha_reg_res)$coeftable), "aggr_J")
    # 
    # # select coefficients to aggregate 
    # select_coefs <- paste0(coefs_to_aggregate, "_X_", original_treatments)
    # 
    # # if the treatment is also one of the crops we want to aggregate over, remove the term, as it is not an indirect effect
    # select_coefs <- select_coefs[select_coefs != paste0(mapmat[mapmat[,"Prices"]==original_treatments, "Crops"], "_X_", original_treatments)]
    # 
    # df_res["aggr_J","Estimate"] <- alpha_reg_res$coefficients[select_coefs] %>% sum()
    # # select the part of the VCOV matrix that is to be used to compute the standard error of the sum 
    # # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
    # df_res["aggr_J","Std. Error"] <- alpha_reg_res$cov.scaled[select_coefs, select_coefs] %>% as.matrix() %>% sum() %>% sqrt()
    # 
    # df_res["aggr_J","t value"]  <- (df_res["aggr_J","Estimate"] - 0)/(df_res["aggr_J","Std. Error"])
    # 
    # df_res["aggr_J","Pr(>|t|)"]  <- (2*pnorm(abs(df_res["aggr_J","t value"]), lower.tail = FALSE)) 
    #}
  }else{
    
    
    ## START A BOOTSTRAP PROCEDURE FROM HERE
    
    make_estimate <- function(myfun_data){  
      gamma_reg_res <- fixest::feols(gamma_model,
                                     data = myfun_data, 
                                     # family = distribution,# "gaussian",#  # "poisson" ,
                                     # se are not computed in fixest anyways, so avoid any more computation than involved the default does
                                     # vcov = se,
                                     # ssc = ssc(cluster.df = "conventional", t.df = "conventional"),
                                     # glm.iter = 25,
                                     #fixef.iter = 100000,
                                     nthreads = 1,
                                     # glm.iter = glm_iter,
                                     notes = TRUE, 
                                     verbose = 4)  
      
      delta_reg_res <- fixest::feols(delta_model,
                                     data = myfun_data, 
                                     # family = distribution,# "gaussian",#  # "poisson" ,
                                     # se are not computed in fixest anyways, so avoid any more computation than involved the default does
                                     # vcov = se,
                                     # ssc = ssc(cluster.df = "conventional", t.df = "conventional"),
                                     # glm.iter = 25,
                                     #fixef.iter = 100000,
                                     nthreads = 1,
                                     # glm.iter = glm_iter,
                                     notes = TRUE, 
                                     verbose = 4)  
      
      gamma_aggr_coeff <- gamma_reg_res$coefficients[regressors] %>% sum()
      #gamma_aggr_se <- gamma_reg_res$cov.scaled[regressors,regressors] %>% as.matrix() %>% sum() %>% sqrt()
      
      delta_aggr_coeff <- delta_reg_res$coefficients[regressors] %>% sum()
      #delta_aggr_se <- delta_reg_res$cov.scaled[regressors,regressors] %>% as.matrix() %>% sum() %>% sqrt()
      
      # # this is our coefficient of interest
      beta <- gamma_aggr_coeff - delta_aggr_coeff
      
      return(beta)
      
      
      # # its SE is conservatively estimated as if delta and gamma coefficients had a null covariance.  
      # df_res[,"std.error"] <- gamma_aggr_se + delta_aggr_se 
      # 
      # df_res[,"observations"] <- gamma_reg_res$nobs
    }
    
    # methodology comes from Cameron and Miller (2015), and:
    # https://stats.stackexchange.com/questions/202916/cluster-boostrap-with-unequally-sized-clusters/202924#202924
    
    # list of parameters related to the dataset, and the clustering variable
    # names and numbers of clusters of size s
    par_list <- list(cluster_variable = boot_cluster, 
                     cluster_names = unique(d_clean[,boot_cluster]),
                     number_clusters = length(unique(d_clean[,boot_cluster])))
    
    start <- Sys.time()
    set.seed(145)
    bootstraped <- boot(data = d_clean, 
                        statistic = make_estimate, 
                        ran.gen = ran.gen_cluster_blc,
                        mle = par_list,
                        sim = "parametric",
                        # parallel = "multicore",
                        # ncpus = detectCores() - 1,
                        R = 400)
    
    end <- Sys.time()
    
    df_res[,"estimate"] <- bootstraped$t0 
    df_res[,"std.error"] <- sd(bootstraped$t)
    df_res[,"observations"] <- nrow(d_clean)
  }
  
  # compute t.stat and p-value of the aggregated estimate (yields similar quantities as in alpha_reg_res if only one coeff was aggregated)
  # df_res <- dplyr::mutate(df_res, t.statistic = (estimate - 0)/(std.error))
  
  # df_res <- dplyr::mutate(df_res, p.value =  (2*pnorm(abs(t.statistic), lower.tail = FALSE)) )
  
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
  #df_res$inference <- se_info  
  
  
  # output wanted
  if(output == "data"){
    toreturn <- list(df_res, d_clean)
  }
  # if(output == "est_obj"){
  #   toreturn <- sum_res
  # }
  if(output == "coef_table"){
    toreturn <- df_res
  }
  
  
  rm(d_clean, reg_res, df_res)
  return(toreturn)
  rm(toreturn)
}

## DEFINE HERE THE VALUES THAT WILL BE PASSED TO SPECIFICATION ARGUMENTS OF REGRESSION FUNCTION ### 
est_parameters <- list(outcome_variable  = "driven_loss",
                       standardization = "_std2",
                       price_info = "_lag1",
                       distribution = "gaussian", 
                       pasture_shares = FALSE)

#### ALPHA EFFECTS #### 


### DISAGGREGATED EFFECTS #### 
# Here, we run the regression over every j-k combinations of interest, and store and plot the estimates 

# infrastructure to store results
alpha_disaggr_res_list <- list()
elm <- 1

for(j_crop in c("Fodder", "Soybean", "Oilpalm")){#, "Cocoa", "Coffee", "Rubber"
  for(k_price in c(mapmat[,"Prices"])){#"Chicken", "Pork", "Sheep", "Crude_oil", 
    # uncomment to prevent estimation from direct effects 
    # if(mapmat[mapmat[,"Crops] == j_crop, "Prices"] != k_price){ 
    alpha_disaggr_res_list[[elm]] <- make_main_reg(original_sj = j_crop, 
                                                   original_treatments = k_price, 
                                                   estimated_effect = "alpha", 
                                                   fe = "grid_id + country_year", # 
                                                   se = "twoway",
                                                   control_interact_sj = TRUE,
                                                   control_interact_Pk = FALSE,
                                                   # reference_crop = c("Oat"),#, "Olive"
                                                   control_direct = FALSE,
                                                   sjpos = FALSE,
                                                   standardization = est_parameters[["standardization"]],
                                                   price_info = est_parameters[["price_info"]],
                                                   distribution = est_parameters[["distribution"]], 
                                                   pasture_shares = est_parameters[["pasture_shares"]]
    )
    
    names(alpha_disaggr_res_list)[elm] <- paste0(j_crop,"_X_",k_price)
    elm <- elm + 1
    # }
  }
}


est_filename <- paste0(est_parameters, collapse = "_") %>% paste0(".Rdata")
est_filename <- paste0("sjinteract_onectrl_sjpos", est_filename)
saveRDS(alpha_disaggr_res_list, here("temp_data","reg_results", "alpha", est_filename))


nodir <- readRDS(here("temp_data","reg_results", "alpha", "indctrl_driven_loss__std2__lag1_gaussian_FALSE.Rdata"))
dir <- readRDS(here("temp_data","reg_results", "alpha", "indctrl_dirctrl_driven_loss__std2__lag1_gaussian_FALSE.Rdata"))

df <- bind_rows(alpha_disaggr_res_list) %>% as.data.frame()
row.names(df) <- paste0(names(alpha_disaggr_res_list))


s <- d_clean[1:999,]
s$first <- c(rep(1, 333), rep(0, 666))
s$second <- c(rep(0, 333), rep(1, 333), rep(0, 333))
s$third <- c(rep(0, 666), rep(1, 333))

feols(fml = as.formula("driven_loss ~ first + third"), data = s)

std <- readRDS(here("temp_data","reg_results", "alpha", "disaggr_driven_loss__std__lag1_gaussian_FALSE.Rdata"))
std1 <- readRDS(here("temp_data","reg_results", "alpha", "disaggr_driven_loss__std1__lag1_gaussian_FALSE.Rdata"))
std2 <- readRDS(here("temp_data","reg_results", "alpha", "disaggr_driven_loss__std2__lag1_gaussian_FALSE.Rdata"))

# let's plot everything first (insignificant ones are a matter of absence of true effect or precision, but not bias)

# prepare regression outputs in a tidy data frame readable by dwplot
df <- bind_rows(alpha_disaggr_res_list)

df$model <- gsub(pattern = "_X_.*$", x = row.names(df), replacement = "") # replace everything after the first underscore with nothing
df$term <- sub(pattern = ".+?(_X_)", x = row.names(df), replacement = "") # replace everyting before the first underscore with nothing

# necessary for dotwhisker to recognize those columns
names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"

title <- paste0("Moderation effects of commodity prices on the influence of agro-ecological suitability on pan-tropical deforestation, NODIR_CTRL 2001-2019") # Cross-effects

# If we want to add brackets on y axis to group k commodities. But not necessarily relevant, as some crops as in several categories. 
# %>%  add_brackets(brackets)
# brackets <- list(c("Oil crops", "Soybean oil", "Palm oil", "Olilve oil", "Rapeseed oil", "Sunflower oil", "Coconut oil"), 
#                  c("Biofuel feedstock", "Sugar", "Maize"))

{dwplot(df,
        dot_args = list(size = 2),
        whisker_args = list(size = 0.5),
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(Beef = "Beef",
                         Chicken = "Chicken",
                         Cocoa = "Cocoa", 
                         Coconut_oil = "Coconut oil", 
                         Coffee = "Coffee", 
                         Crude_oil = "Crude oil",
                         Maize = "Maize",
                         Olive_oil = "Olive oil", 
                         Palm_oil = "Palm oil",
                         Pork = "Pork", 
                         Rapeseed_oil = "Rapeseed oil", 
                         Rice = "Rice",
                         Rubber = "Rubber",
                         Sheep = "Sheep", 
                         Sorghum = "Sorghum",
                         Soybean = "Soybean",
                         Soybean_meal = "Soybean meal", 
                         Sugar = "Sugar", 
                         Sunflower_oil = "Sunflower oil", 
                         Tea = "Tea", 
                         Tobacco = "Tobacco",
                         Wheat = "Wheat"
    )) +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    ggtitle(title) +
    theme(plot.title = element_text(face="bold", size=c(10)),
          legend.position = c(0.007, 0.001),
          legend.justification = c(0, 0), 
          legend.background = element_rect(colour="grey80"),
          legend.title = element_blank())}  


### EFFECTS AGGREGATED OVER J MAIN DRIVERS #### 
# Here, we run the regression over every J-k combinations of interest, and store and plot the estimates 

# infrastructure to store results
alpha_aggr_J_res_list <- list()
elm <- 1

for(k_price in c(mapmat[,"Prices"], "Chicken", "Pork", "Sheep", "Crude_oil")){
  
  alpha_aggr_J_res_list[[elm]] <- make_main_reg(original_sj = c("Fodder", "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber"), 
                                                original_treatments = k_price)
  
  names(alpha_aggr_J_res_list)[elm] <- paste0("6 main drivers_X_",k_price)
  elm <- elm + 1
}

est_filename <- paste0(est_parameters, collapse = "_") %>% paste0(".Rdata")
est_filename <- paste0("aggr_J", est_filename)
saveRDS(alpha_aggr_J_res_list, here("temp_data","reg_results", "alpha", est_filename))

# let's plot everything first (insignificant ones are a matter of absence of true effect or precision, but not bias)

# prepare regression outputs in a tidy data frame readable by dwplot
df <- bind_rows(alpha_aggr_J_res_list)

df$model <- gsub(pattern = "_X_.*$", x = row.names(df), replacement = "") # replace everything after the first underscore with nothing
df$term <- sub(pattern = ".+?(_X_)", x = row.names(df), replacement = "") # replace everyting before the first underscore with nothing

# necessary for dotwhisker to recognize those columns
names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"

title <- paste0("Moderation effects of commodity prices on the influence of agro-ecological suitability on pan-tropical deforestation, 2001-2019, aggregated over the 6 main drivers") # Cross-effects

# If we want to add brackets on y axis to group k commodities. But not necessarily relevant, as some crops as in several categories. 
# %>%  add_brackets(brackets)
# brackets <- list(c("Oil crops", "Soybean oil", "Palm oil", "Olilve oil", "Rapeseed oil", "Sunflower oil", "Coconut oil"), 
#                  c("Biofuel feedstock", "Sugar", "Maize"))

{dwplot(df,
        dot_args = list(size = 2),
        whisker_args = list(size = 0.5),
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(Beef = "Beef",
                         Chicken = "Chicken",
                         Cocoa = "Cocoa", 
                         Coconut_oil = "Coconut oil", 
                         Coffee = "Coffee", 
                         Crude_oil = "Crude oil",
                         Maize = "Maize",
                         Olive_oil = "Olive oil", 
                         Palm_oil = "Palm oil",
                         Pork = "Pork", 
                         Rapeseed_oil = "Rapeseed oil", 
                         Rice = "Rice",
                         Rubber = "Rubber",
                         Sheep = "Sheep", 
                         Sorghum = "Sorghum",
                         Soybean = "Soybean",
                         Soybean_meal = "Soybean meal", 
                         Sugar = "Sugar", 
                         Sunflower_oil = "Sunflower oil", 
                         Tea = "Tea", 
                         Tobacco = "Tobacco",
                         Wheat = "Wheat"
    )) +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    ggtitle(title) +
    theme(plot.title = element_text(face="bold", size=c(10)),
          legend.position = c(0.007, 0.001),
          legend.justification = c(0, 0), 
          legend.background = element_rect(colour="grey80"),
          legend.title = element_blank())}  


### EFFECTS AGGREGATED OVER K PRICES ####
# Here, we run the regression over every j-K combinations of interest, and store and plot the estimates 

# infrastructure to store results
alpha_aggr_K_res_list <- list()
elm <- 1

for(j_crop in c("Fodder", "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber")){#
  
  alpha_aggr_K_res_list[[elm]] <- make_main_reg(original_sj = j_crop, 
                                                original_treatments = c(mapmat[,"Prices"], "Chicken", "Pork", "Sheep", "Crude_oil"))
  
  names(alpha_aggr_K_res_list)[elm] <- paste0(j_crop, "_X_all prices")
  elm <- elm + 1
}

est_filename <- paste0(est_parameters, collapse = "_") %>% paste0(".Rdata")
est_filename <- paste0("aggr_K", est_filename)
saveRDS(alpha_aggr_K_res_list, here("temp_data","reg_results", "alpha", est_filename))

# let's plot everything first (insignificant ones are a matter of absence of true effect or precision, but not bias)

# prepare regression outputs in a tidy data frame readable by dwplot
df <- bind_rows(alpha_aggr_K_res_list)

df$model <- gsub(pattern = "_X_.*$", x = row.names(df), replacement = "") # replace everything after the first underscore with nothing
df$term <- sub(pattern = ".+?(_X_)", x = row.names(df), replacement = "") # replace everyting before the first underscore with nothing

# necessary for dotwhisker to recognize those columns
names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"

title <- paste0("Moderation effects of commodity prices on the influence of agro-ecological suitability on pan-tropical deforestation, 2001-2019, aggregated over 28 commodity prices") # Cross-effects

# If we want to add brackets on y axis to group k commodities. But not necessarily relevant, as some crops as in several categories. 
# %>%  add_brackets(brackets)
# brackets <- list(c("Oil crops", "Soybean oil", "Palm oil", "Olilve oil", "Rapeseed oil", "Sunflower oil", "Coconut oil"), 
#                  c("Biofuel feedstock", "Sugar", "Maize"))

{dwplot(df,
        dot_args = list(size = 2),
        whisker_args = list(size = 0.5),
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(Beef = "Beef",
                         Chicken = "Chicken",
                         Cocoa = "Cocoa", 
                         Coconut_oil = "Coconut oil", 
                         Coffee = "Coffee", 
                         Crude_oil = "Crude oil",
                         Maize = "Maize",
                         Olive_oil = "Olive oil", 
                         Palm_oil = "Palm oil",
                         Pork = "Pork", 
                         Rapeseed_oil = "Rapeseed oil", 
                         Rice = "Rice",
                         Rubber = "Rubber",
                         Sheep = "Sheep", 
                         Sorghum = "Sorghum",
                         Soybean = "Soybean",
                         Soybean_meal = "Soybean meal", 
                         Sugar = "Sugar", 
                         Sunflower_oil = "Sunflower oil", 
                         Tea = "Tea", 
                         Tobacco = "Tobacco",
                         Wheat = "Wheat"
    )) +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    ggtitle(title) +
    theme(plot.title = element_text(face="bold", size=c(10)),
          legend.position = c(0.007, 0.001),
          legend.justification = c(0, 0), 
          legend.background = element_rect(colour="grey80"),
          legend.title = element_blank())}  



#### BETA EFFECTS #### 

### DISAGGREGATED EFFECTS #### 
# Here, we run the regression over every j-k combinations of interest, and store and plot the estimates 

# infrastructure to store results
beta_disaggr_res_list <- list()
elm <- 1

for(j_crop in c("Soybean")){#"Fodder", "Oilpalm", , "Cocoa", "Coffee", "Rubber"
  for(k_price in c("Beef", "Pork", "Maize", "Wheat", "Orange")){#c(mapmat[,"Prices"], "Chicken", "Pork", "Sheep", "Crude_oil")
    # this makes sure that we estimate only models where cross-commodity terms are in the regressors
    if(mapmat[mapmat[,"Crops"] == j_crop, "Prices"] != k_price){  
      beta_disaggr_res_list[[elm]] <- make_main_reg(original_sj = j_crop, 
                                                    original_treatments = k_price, 
                                                    estimated_effect = "beta", 
                                                    fe = "grid_id",
                                                    boot_cluster = "grid_id",
                                                    standardization = est_parameters[["standardization"]],
                                                    price_info = est_parameters[["price_info"]],
                                                    distribution = est_parameters[["distribution"]], 
                                                    pasture_shares = est_parameters[["pasture_shares"]]
      )
      
      names(beta_disaggr_res_list)[elm] <- paste0(j_crop,"_X_",k_price)
      elm <- elm + 1
    }
  }
}

est_filename <- paste0(est_parameters, collapse = "_") %>% paste0(".Rdata")
est_filename <- paste0("disaggr_", est_filename)
saveRDS(beta_disaggr_res_list, )

old_res <- readRDS(here("temp_data","reg_results", "beta", paste0("unblc_",est_filename)))

std2 <- readRDS(here("temp_data","reg_results", "beta", "disaggr_driven_loss__std2__lag1_gaussian_FALSE.Rdata"))




