

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


### MAIN DATA SET ### 
main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "driverloss_aesi_long_final.Rdata"))

### PRICE DATA ###
prices <- readRDS(here("temp_data", "prepared_international_prices.Rdata"))

### HELPER FUNCTION TO BOOTSTRAP CLUSTER SE 
# will tell boot::boot how to sample data at each replicate of statistic
ran.gen_cluster <- function(original_data, arg_list){
  
  cl_boot_dat <- NULL
  
  cluster_var <- arg_list[["cluster_variable"]]
  
  # non-unique names of clusters (repeated when there is more than one obs. in a cluster) 
  nu_cl_names <- as.character(original_data[,cluster_var]) 
  
  for(s in arg_list[["unique_sizes"]]){
    # sample, in the vector of names of clusters of size s, as many draws as there are clusters of that size, with replacement
    sample_cl_s <- sample(arg_list[["cluster_names"]][[s]], 
                          arg_list[["number_clusters"]][[s]], 
                          replace = TRUE) 
    
    # because of replacement, some names are sampled more than once
    sample_cl_s_tab <- table(sample_cl_s)
    
    # we need to give them a new cluster identifier, otherwise a cluster sampled more than once 
    # will be "incorrectly treated as one large cluster rather than two distinct clusters" (by the fixed effects) (Cameron 2015)    
    
    for(n in 1:max(sample_cl_s_tab)){ # from 1 to the max number of times a name was sampled bc of replacement
      # vector to select obs. that are within the sampled clusters. 
      sel <- nu_cl_names %in% names(sample_cl_s_tab[sample_cl_s_tab == n])
      
      # replicate the selected data n times
      clda <- original_data[sel,][rep(seq_len(nrow(original_data[sel,])), n), ]
      
      #identify row names without periods, and add ".0" 
      row.names(clda)[grep("\\.", row.names(clda), invert = TRUE)] <- paste0(grep("\\.", row.names(clda), invert = TRUE, value = TRUE),".0")
      
      # add the suffix due to the repetition after the existing cluster identifier. 
      clda[,cluster_var] <- paste0(clda[,cluster_var], sub(".*\\.","_",row.names(clda)))
      
      # stack the bootstrap samples iteratively 
      cl_boot_dat <- rbind(cl_boot_dat, clda) 
    }
  }  
  return(cl_boot_dat)
}

### TEMPORARY OBJECTS 
outcome_variable = "driven_loss" # "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
start_year = 2001
end_year = 2019
continent = "all"
further_lu_evidence = "none"
original_exposures = c("Soybean") # ,"Fodder",  "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber" in GAEZ spelling, one, part, or all of the 6 main drivers of deforestation: 

# for the k variables hypothesized in overleaf for palm oil, feglm quasipoisson converges within 25 iter.
original_treatments <-  c("Beef") # in price spelling. One, part, or all of the full set of commodities having a price-AESI match.
# "Beef" , , , "Cocoa", "Coffee", "Rubber"
# "Rapeseed_oil", "Sunflower_oil","Rice", "Wheat", "Maize", "Sugar", "Sorghum")
# mapmat[,"Prices"], "Chicken", "Pork", "Sheep", "Crude_oil"
extra_price_k = c() # in price spelling. One, part, or all of the full set of commodities NOT having a price-AESI match.
# ,"Chicken", "Pork", "Crude_oil

pasture_shares <- FALSE
j_soy = "Soybean"
fcr = 7.2
standardization = "_std2"
price_info = "_lag1"
estimated_effect = "beta"
# sjpj_lag = "_lag1" # either "" or "_lag1" or "_lag2"
# skpk_lag = "_lag1" # either "" or "_lag1" or "_lag2"SjPj = TRUE
# SjPj = TRUE
# SkPk = TRUE
remaining <- TRUE
# open_path <- FALSE
# commoXcommo <- "Fodder_X_Beef"
fe = "grid_id" # + country_year
distribution = "gaussian"
invhypsin = TRUE
conley_cutoff <- 100
se = "twoway"
boot_cluster ="grid_id"
#coefstat = "confint"
output = "coef_table"
glm_iter <- 25 



make_main_reg <- function(outcome_variable = "driven_loss", # one of "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
                          start_year = 2001, 
                          end_year = 2020, 
                          continent = "all", # one of "Africa", "America", "Asia", or "all"
                          further_lu_evidence = "none", # either "none", "sbqt_direct_lu", or "sbqt"mode_lu"
                          original_exposures = c("Fodder"),  # in GAEZ spelling, one, part, or all of the 6 main drivers of deforestation: 
                          # , "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber"
                          # note that this is the vector that defines what is removed from full_control
                          original_treatments = c("Soybean"), # in price spelling. One, part, or all of the full set of commodities having a price-AESI match.
                          # , "Palm_oil", "Cocoa", "Coffee", "Rubber", 
                          # "Rapeseed_oil", "Sunflower_oil","Rice", "Wheat", "Maize", "Sugar", "Sorghum"
                          # extra_price_k = c(), # "Chicken", "Pork", "Sheep", "Crude_oil" in price spelling. One, part, or all of the full set of commodities NOT having a price-AESI match.
                          pasture_shares = FALSE, # if TRUE, and crop_j = "Fodder", then qj is proxied with the share of pasture area in 2000. 
                          standardization = "_std2", # one of "", "_std", or "_std2"
                          price_info = "_lag1", # one of "lag1", "_2pya", "_3pya", "_4pya", "_5pya",
                          estimated_effect = "alpha",# if "alpha", estimates the (aggregated or not) cross-elasticity effect. If different from "alpha" (e.g. "gamma") it estimates the ILUC effect from price covariation. 
                          # SjPj = FALSE,
                          # SkPk = TRUE,
                          # sjpj_lag = "_lag1", # either "" or "_lag1" or "_lag2"
                          # skpk_lag = "_lag1", # either "" or "_lag1" or "_lag2"
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
                          output = "coef_table" # one of "data", or "coef_table" 
){
  
  
  #### PREPARE NEEDED VARIABLE NAMES
  # this does not involve data, just arguments of the make_reg function
  # original_ names are used to get generic covariate names in coefficient tables (irrespective of modelling choices)
  
  ## Standardized suitability index to find in the main data 
  # this just those for which focus is set in current specification
  exposures <- paste0(original_exposures, standardization)  
  # this is all possible exposures, necessary in every specificaton 
  original_all_exposures <- mapmat[,"Crops"]
  all_exposures <- paste0(original_all_exposures, standardization)
  
  ## Price variable names to find in the price data
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
                                 unique(c(exposures, all_exposures))))) #
  
  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", unique(c(treatments, all_treatments)))], by = c("year"))
  
  
  # If we want the exposure to deforestation for pasture to be proxied with the share of pasture area in 2000, rather than suitability index, 
  if(pasture_shares){
    d[,grepl("Fodder", names(d))] <- d$pasture_share
  }
  
  # transform dependent variable, if gaussian GLM
  if(distribution == "gaussian" & invhypsin){
    d <- dplyr::mutate(d, !!as.symbol(outcome_variable) := asinh(!!as.symbol(outcome_variable)))
  }
  
  
  ### Main regressors
  
  regressors <- c()
  for(Pk in treatments){
    for(Sj in exposures){
      varname <- paste0(original_exposures[match(Sj, exposures)], "_X_", original_treatments[match(Pk, treatments)])
      regressors <- c(regressors, varname)
      # Log them so their variations are comparable
      d <- mutate(d, 
                  !!as.symbol(varname) := log( (!!as.symbol(Sj)) * (!!as.symbol(Pk)) +1))
      
    }
  }
  rm(varname, Sj, Pk)
  
  
  
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
      for(Sj in all_exposures){
        varname <- paste0(original_all_exposures[match(Sj, all_exposures)], "_X_", original_all_treatments[match(Pk, all_treatments)])
        # indiv_controls <- c(indiv_controls, varname)
        # d <- mutate(d, 
        #             !!as.symbol(varname) := log( (!!as.symbol(Sj)) * (!!as.symbol(Pk)) +1))
        # note that mutate will overwrite the already existing varnames created in the regressors step above. 
        
        # we need to add up the new terms together as we create them, otherwise there are 529 new columns added to the 4M+ rows dataframe and this is memory intensive
        if(!(varname%in%regressors)){ # add the interaction term if it's not a regressor already
          d <- mutate(d, full_control := ( full_control + log( (!!as.symbol(Sj)) * (!!as.symbol(Pk)) +1 ) ) )
        }
      }
    }
    
    # c(log(d[1,"Banana_lag1"]*d[1,"Barley_std2"] +1), log(d[1,"Barley_lag1"]*d[1,"Barley_std2"] +1))%>%sum()
    
    rm(varname, Pk, Sj)
    
    delta_controls <- c(delta_controls, "full_control")
    
    ## Partial control (to estimate gamma)
    # identify the variables to put in the partial control term, but this changes depending on the question we are answering.
    # yet, it is not necessary to condition: the loop handles it, 
    # bc in all cases, the removed controls capture the confounding covariation between the treatment and the prices of the main crops directly driving deforestation 
    d$part_control <- d$full_control
    for(Sj in exposures){
      # from the inputed Sj, find the corresponding price to match it to, construct the SjPj interaction, and remove it from the control. 
      original_Sj <- original_exposures[match(Sj,exposures)]
      original_Pj <- mapmat[mapmat[,"Crops"]==original_Sj,"Prices"]
      Pj <- all_treatments[match(original_Pj,original_all_treatments)]
      # it's important that part_control be on both sides of the equal, and that it starts from full_control, to handle cases when original_exposures has several elements. 
      d <- mutate(d, part_control := ( part_control - (log( (!!as.symbol(Sj)) * (!!as.symbol(Pj)) +1) ) ) ) # (note that the regressors are already taken out of full_control)
    }
    
    rm(original_Sj, original_Pj, Sj, Pj)
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
  
  # - are in study period 
  d <- dplyr::filter(d, year >= start_year)
  d <- dplyr::filter(d, year <= end_year)
  
  # - are in study area
  if(continent != "all"){
    d <- dplyr::filter(d, continent_name == continent)
  }
  
  # have remaining forest
  d <- dplyr::filter(d, remaining_fc > 0)
  
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
  
  # - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove
  obstormv <- obs2remove(fml = as.formula(paste0(outcome_variable, " ~ ", fe)),
                         d, 
                         family = "poisson")
  # this is necessary to handle cases when there is no obs. to remove
  if(length(obstormv)>0){
    d_clean <- d[-obstormv,]
  } else {
    d_clean <- d
  }
  
  rm(d, obstormv)
  
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
    
    # select the coefficients to be addep up (code works in cases when there is only one regressor because we investigate specific crop-crop interactions)
    # (I have checked it)
    df_res[,"estimate"] <- alpha_reg_res$coefficients[regressors] %>% sum()
    # select the part of the VCOV matrix that is to be used to compute the standard error of the sum 
    # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
    df_res[,"std.error"] <- alpha_reg_res$cov.scaled[regressors,regressors] %>% as.matrix() %>% sum() %>% sqrt()
    
    df_res[,"observations"] <- alpha_reg_res$nobs
    
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
    
    # get the different cluster sizeS. This is necessary to cluster bootstraping with clusters of different sizes. 
    # methodology comes from:
    # https://stats.stackexchange.com/questions/202916/cluster-boostrap-with-unequally-sized-clusters/202924#202924
    
    # this is a vector of same length as cluster_names. Values are the number of observations (the size) per cluster. 
    sizes <- table(d_clean[,boot_cluster]) # in original code this was not grid_id but boot_cluster, which was "parcel_id". 
    u_sizes <- sort(unique(sizes))
    cl_names <- list()
    n_clusters <- list()
    for(s in u_sizes){
      # names and numbers of clusters of size s
      cl_names[[s]] <- names(sizes[sizes == s])
      n_clusters[[s]] <- length(cl_names[[s]])
    }
    #sum(unlist(n_clusters)) == length(unique(d_clean[,boot_cluster]))
    
    # start <- Sys.time()
    set.seed(145)
    bootstraped <- boot(data = d_clean, 
                        statistic = make_estimate, 
                        ran.gen = ran.gen_cluster,
                        mle = list(cluster_variable = boot_cluster, 
                                   unique_sizes = u_sizes, 
                                   cluster_names = cl_names, 
                                   number_clusters = n_clusters),
                        sim = "parametric",
                        # parallel = "multicore",
                        # ncpus = detectCores() - 1,
                        R = 200)
    
    # end <- Sys.time()
    
    df_res[,"estimate"] <- bootstraped$t0 
    df_res[,"std.error"] <- sd(bootstraped$t)
    df_res[,"observations"] <- nrow(d_clean)
  }
  
  # compute t.stat and p-value of the aggregated estimate (yields similar quantities as in alpha_reg_res if only one coeff was aggregated)
  df_res <- dplyr::mutate(df_res, t.statistic = (estimate - 0)/(std.error))
  
  df_res <- dplyr::mutate(df_res, p.value =  (2*pnorm(abs(t.statistic), lower.tail = FALSE)) )
  
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
  df_res$inference <- se_info  
  
  
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



### DEFINE HERE THE VALUES THAT WILL BE PASSED TO SPECIFICATION ARGUMENTS OF REGRESSION FUNCTION ### 
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

for(j_crop in c("Fodder", "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber")){#
  for(k_price in c(mapmat[,"Prices"], "Chicken", "Pork", "Sheep", "Crude_oil")){
    # uncomment to prevent estimation from direct effects 
    # if(mapmat[mapmat$Crops == j_crop, "Prices"] == k_price){ 
    alpha_disaggr_res_list[[elm]] <- make_main_reg(original_exposures = j_crop, 
                                                   original_treatments = k_price, 
                                                   estimated_effect = "alpha", 
                                                   fe = "grid_id + country_year", 
                                                   se = "twoway",
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
est_filename <- paste0("disaggr_", est_filename)
saveRDS(alpha_disaggr_res_list, here("temp_data","reg_results", "alpha", est_filename))

# let's plot everything first (insignificant ones are a matter of absence of true effect or precision, but not bias)

# prepare regression outputs in a tidy data frame readable by dwplot
df <- bind_rows(alpha_disaggr_res_list)

df$model <- gsub(pattern = "_X_.*$", x = row.names(df), replacement = "") # replace everything after the first underscore with nothing
df$term <- sub(pattern = ".+?(_X_)", x = row.names(df), replacement = "") # replace everyting before the first underscore with nothing

# necessary for dotwhisker to recognize those columns
names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"

title <- paste0("Moderation effects of commodity prices on the influence of agro-ecological suitability on pan-tropical deforestation, 2001-2019") # Cross-effects

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
  
  alpha_aggr_J_res_list[[elm]] <- make_main_reg(original_exposures = c("Fodder", "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber"), 
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
  
  alpha_aggr_K_res_list[[elm]] <- make_main_reg(original_exposures = j_crop, 
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
  for(k_price in c("Maize", "Wheat")){#c(mapmat[,"Prices"], "Chicken", "Pork", "Sheep", "Crude_oil")
    # this makes sure that we estimate only models where cross-commodity terms are in the regressors
    if(mapmat[mapmat$Crops == j_crop, "Prices"] == k_price){  
      beta_disaggr_res_list[[elm]] <- make_main_reg(original_exposures = j_crop, 
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
saveRDS(beta_disaggr_res_list, here("temp_data","reg_results", "beta", est_filename))






