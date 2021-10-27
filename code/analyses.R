
##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
neededPackages = c("plyr", "dplyr", "here", #"tibble", "data.table",
                   "foreign", # "readxl",
                   "raster", "rgdal",  "sp", "sf", # "spdep",
                   "DataCombine",
                   "knitr", "kableExtra",
                   "fixest", #,"msm", "car",  "sandwich", "lmtest", "boot", "multcomp",
                   "ggplot2", "dotwhisker", #"tmap",# "leaflet", "htmltools"
                   "foreach"
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
  "Orange", "Citrus",
  "Cocoa", "Cocoa",
  "Coconut_oil", "Coconut",
  "Coffee", "Coffee",
  "Cotton", "Cotton",
  "Groundnut", "Groundnut",
  "Maize", "Maizegrain",
  "Oat", "Oat",
  "Orange", "Citrus", # Citrus sinensis in both GAEZ and FAO
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

# crop_prices <- list("Fodder" = c("Soybean", "Palm_oil", "Cocoa", "Coffee", "Rubber", "Maize", "Sugar", "Wheat", "Barley", "Oat"), 
#                     "Soybean" = c("Beef", "Palm_oil", "Cocoa", "Coffee", "Rubber", "Rapeseed_oil", "Sunflower_oil",  "Maize", "Sugar", "Wheat"),
#                     "Oilpalm" = c("Beef", "Soybean", "Cocoa", "Coffee", "Rubber", "Rapeseed_oil", "Sunflower_oil", "Maize", "Sugar", "Wheat"),
#                     "Cocoa" = c("Beef", "Soybean", "Palm_oil", "Coffee", "Rubber", "Rapeseed_oil", "Sunflower_oil", "Sugar"),
#                     "Coffee" = c("Beef", "Soybean", "Palm_oil", "Cocoa", "Rubber", "Tea", "Sugar", "Tobacco"),
#                     "Rubber" = c("Beef", "Soybean", "Palm_oil", "Coffee", "Cocoa", "Rapeseed_oil", "Sunflower_oil",  "Maize", "Sugar", "Wheat", "Barley", "Oat")
# )

# We do not run robustness checks for Beef and Coffee, as they have nothing significant in the main estimation. 
# crops_to_runover <- c("Soybean", "Oilpalm", "Cocoa", "Rubber")


### READ ALL POSSIBLE DATASETS HERE 
# but not all together because of memory issues. 
# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aesi_long_final.Rdata"))

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_aesi_long_final.Rdata"))

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aesi_long_final.Rdata"))

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aeay_long_final.Rdata"))

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_aeay_long_final.Rdata"))

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aeay_long_final.Rdata"))

prices <- readRDS(here("temp_data", "prepared_international_prices.Rdata"))
# 
# rm(outcome_variable, start_year, end_year, crop_j, j_soy, price_k, extra_price_k, SjPj, SkPk, fe, distribution, output, se, cluster, 
#    controls, regressors, outcome_variable)

outcome_variable = "driven_loss" # "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
start_year = 2001
end_year = 2019
continent = "all"
further_lu_evidence = "none"
crop_j = c("Fodder", "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber") # in GAEZ spelling
j_soy = "Soybean"
standardized = FALSE
fcr = 7.2
# for the k variables hypothesized in overleaf for palm oil, feglm quasipoisson converges within 25 iter.
# but maybe not with skPk controls.
price_k <- c("Beef", "Soybean", "Palm_oil", "Cocoa", "Coffee", "Rubber")#, 
             # "Rapeseed_oil", "Sunflower_oil","Rice", "Wheat", "Maize", "Sugar", "Sorghum")

extra_price_k = c() # , "Sheep","Chicken", "Pork"
price_info = "lag1"
sjpj_lag = "_lag2" # either "" or "_lag1" or "_lag2"
skpk_lag = "_lag1" # either "" or "_lag1" or "_lag2"SjPj = TRUE
# SjPj = TRUE
# SkPk = TRUE
remaining <- TRUE
pasture_shares <- FALSE
open_path <- FALSE
commoXcommo <- "Fodder_X_Beef"
fe = "grid_id + country_year"
distribution = "quasipoisson"
conley_cutoff <- 100
se = "twoway"
cluster ="grid_id"
#coefstat = "confint"
output = "coef_table"
glm_iter <- 25



make_main_reg <- function(outcome_variable = "driven_loss", # one of "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
                          start_year = 2001, 
                          end_year = 2020, 
                          continent = "all", # one of "Africa", "America", "Asia", or "all"
                          further_lu_evidence = "none", # either "none", "sbqt_direct_lu", or "sbqt"mode_lu"
                          crop_j = c("Fodder", "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber"), # in GAEZ spelling
                          price_k = c("Beef", "Soybean", "Palm_oil", "Cocoa", "Coffee", "Rubber", 
                                      "Rapeseed_oil", "Sunflower_oil","Rice", "Wheat", "Maize", "Sugar", "Sorghum"), # in prices spelling
                          extra_price_k = c(), # One of "Crude_oil", "Chicken", "Pork", "Sheep" 
                          SjPj = TRUE,
                          SkPk = TRUE,
                          price_info = "lag1", # one of "lag1", "2pya", "3pya", "4pya", "5pya",
                          sjpj_lag = "_lag1", # either "" or "_lag1" or "_lag2"
                          skpk_lag = "_lag1", # either "" or "_lag1" or "_lag2"
                          remaining = TRUE, # should remaining forest be controlled for 
                          pasture_shares = FALSE, # if TRUE, and crop_j = "Fodder", then qj is proxied with the share of pasture area in 2000. 
                          open_path = FALSE,
                          commoXcommo = "Fodder_X_Beef",
                          #commo_m = c(""), comment coder ça pour compatibilité avec loops over K_commo ? 
                          fe = "grid_id + country_year", 
                          distribution = "quasipoisson",#  "quasipoisson", 
                          se = "twoway", # passed to vcov argument in fixest::summary. Currently, one of "cluster", "twoway", or "conley".  
                          conley_cutoff = 100, # the distance cutoff, in km, passed to fixest::vcov_conley, if se = "conley"  
                          cluster ="grid_id", # the cluster level if se = "cluster" (i.e. one way)
                          # coefstat = "confint", # one of "se", "tstat", "confint"
                          glm_iter = 25,
                          output = "est_obj" # one of "data", "est_obj", "coef_table" 
){
  
  
  #### PREPARE NEEDED VARIABLE NAMES
  # this does not involve data, just arguments of the make_reg function 
  
  ## Suitability index variable names
  # Explicit the names of the sk variables corresponding to commodities in price_k only (i.e. not in price_j nor extra_price_k)
  # match function is important because suitability_k must have commodities in the same order as price_k 
  crop_k <- mapmat[,"Crops"][match(price_k, mapmat[,"Prices"])]

  # Standardized suitability index or not
  if(standardized){
    suitability_j <- paste0(crop_j, "_std")  
    suitability_k <- paste0(crop_k, "_std")  
  } else{
    suitability_j <- crop_j
    suitability_k <- crop_k
  }
  
  ## Price variable names
  # Explicit the name of the price of j based on crop_j (GAEZ spelling)
  price_j <- mapmat[match(crop_j, mapmat[,"Crops"]),"Prices"]

  # Group the names of the different prices 
  original_price_j <- price_j
  original_price_k <- price_k # /!\ this one is used in making rk 
  original_price_reg <- c(price_k, extra_price_k)
  
  # lag the different price sets
  price_reg <- paste0(original_price_reg, "_", price_info)
  price_j <- paste0(original_price_j, sjpj_lag)
  price_k <- paste0(original_price_k, skpk_lag)
  
  
  #### MAKE THE VARIABLES NEEDED IN THE DATA
  #d <- main_data
  if(outcome_variable == "first_loss" | outcome_variable == "nd_first_loss"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aesi_long_final.Rdata"))}
  if(outcome_variable == "firstloss_glassgfc"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_aesi_long_final.Rdata"))}
  if(outcome_variable == "phtf_loss"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aesi_long_final.Rdata"))}
  if(outcome_variable == "driven_loss"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "driverloss_aesi_long_final.Rdata"))}
  
  # Keep only in data the useful variables 
  d <- dplyr::select(d, all_of(c("grid_id", "year", "lat", "lon", "country_name", "country_year", "continent_name", "remaining_fc", "accu_defo_since2k",
                                 outcome_variable, "pasture_share",
                                 suitability_j, suitability_k))) 
  
  # If we want the exposure to deforestation for pasture to be proxied with the share of pasture area in 2000, rather than suitability index, 
  if(pasture_shares){
    d[,grepl("Fodder", names(d))] <- d$pasture_share
  }
  
  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", unique(c(price_reg, price_j, price_k)))], by = c("year"))
  
  # Main regressors
  regressors <- c()
  for(Pk in price_reg){
    for(Sj in suitability_j){
      varname <- paste0(crop_j[match(Sj, suitability_j)], "_X_", original_price_reg[match(Pk, price_reg)])
      regressors <- c(regressors, varname)
      d <- mutate(d, 
                  !!as.symbol(varname) := (!!as.symbol(Sj)*(!!as.symbol(Pk))))
    }
  }
  rm(varname, Sj, Pk)
  
  
  ## Mechanisms
  # direct_effects_var <- paste0(mapmat[,"Crops"], "_X_", mapmat[,"Prices"]) 
  if(open_path){
    regressors <- regressors[!(regressors %in% commoXcommo)]
  }
  
  
  ## Log them so their variations are comparable
  for(reg in regressors){
    d <- mutate(d, !!as.symbol(reg) := log(!!as.symbol(reg)+1))
  }
    
  ## Controls
  controls <- c() # it's important that this is not conditioned on SkPk nor sjPj
  
  # add remainging forest cover as a control
  if(remaining){
    controls <- c(controls, "remaining_fc")
  }
  
  
  # MODEL SPECIFICATION FORMULAE
  if(length(controls) > 0){
    fe_model <- as.formula(paste0(outcome_variable,
                                  " ~ ",
                                  paste0(regressors, collapse = "+"),
                                  " + ",
                                  paste0(controls, collapse = "+"),
                                  " | ",
                                  fe))
  }else{
    fe_model <- as.formula(paste0(outcome_variable,
                                  " ~ ",
                                  paste0(regressors, collapse = "+"),
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
  
  used_vars <- c("grid_id", "year", "lat", "lon", "country_name", "country_year", "continent_name",
                 outcome_variable, regressors, controls)
  
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
  
  # is.na(d$Oilpalm) %>% sum()
  
  
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
  
  if(distribution != "negbin"){ # i.e. if it's poisson or quasipoisson or gaussian
    reg_res <- fixest::feglm(fe_model,
                             data = d_clean, 
                             family = distribution,# "gaussian",#  # "poisson" ,
                             # glm.iter = 25,
                             #fixef.iter = 100000,
                             nthreads = 3,
                             glm.iter = glm_iter,
                             notes = TRUE, 
                             verbose = 4)  
    
    
  }
  
  # Now keep only information necessary, otherwise the output of fixest estimation is large and we can't collect too many at the same time (over loops)  
  # this is necessary to compute SE as we want to.  
  if(se == "conley"){
    sum_res <- summary(reg_res,
                       vcov = vcov_conley(lat = "lat", lon = "lon", cutoff = conley_cutoff, distance = "spherical"))
    
    # store info on SE method
    se_info <- paste0("Conley (",conley_cutoff,"km)")
  }
  if(se == "twoway"){
    # " If the two variables were used as fixed-effects in the estimation, 
    # you can leave it blank with vcov = "twoway""
    sum_res <- summary(reg_res,
                       vcov = se)
    
    se_info <- "Clustered (grid cell - country*year"
  }
  if(se == "cluster"){
    sum_res <- summary(reg_res,
                       #vcov = se, 
                       cluster = cluster)
    
    se_info <- paste0("Clustered (",cluster,")")
  }
  
  df_res <- sum_res$coeftable %>% as.data.frame()
  # add a column with the number of observations
  df_res$Observations <- sum_res$nobs
  df_res$`Standard errors` <- se_info
  
  # output wanted
  if(output == "data"){
    toreturn <- list(reg_res, d_clean)
  }
  if(output == "est_obj"){
    toreturn <- sum_res
  }
  if(output == "coef_table"){
    toreturn <- df_res
  }
  
  
  rm(d_clean, reg_res, df_res)
  return(toreturn)
  rm(toreturn)
}


#### MAIN FIGURE #### 
df <- df_res


df <- make_main_reg(crop_j = c("Fodder", "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber"),
                         price_k =  c("Beef", "Soybean", "Palm_oil", "Cocoa", "Coffee", "Rubber"),
                         output = "coef_table"
)


# prepare regression outputs in a tidy data frame readable by dwplot

# select estimates to display (significant ones)
df <- df[df$`Pr(>|t|)`<0.1,]
df <- df[rownames(df)!="remaining_fc",]
# df <- df[c("Soybean_X_Beef", 
#            "Soybean_X_Palm_oil",
#            "Soybean_X_Sugar", 
#            "Soybean_X_Maize",
#            "Oilpalm_X_Sugar",  
#            "Oilpalm_X_Sunflower_oil",
#            "Cocoa_X_Coffee",
#            "Cocoa_X_Rapeseed_oil", 
#            "Rubber_X_Sugar"),]

df$model <- gsub(pattern = "_X_.*$", x = row.names(df), replacement = "") # replace everything after the first underscore with nothing
df$term <- sub(pattern = ".+?(_X_)", x = row.names(df), replacement = "") # replace everyting before the first underscore with nothing

names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"

title <- paste0("Indirect effects of commodity prices on agriculture-driven tropical deforestation")

# If we want to add brackets on y axis to group k commodities. But not necessarily relevant, as some crops as in several categories. 
# %>%  add_brackets(brackets)
# brackets <- list(c("Oil crops", "Soybean oil", "Palm oil", "Olilve oil", "Rapeseed oil", "Sunflower oil", "Coconut oil"), 
#                  c("Biofuel feedstock", "Sugar", "Maize"))
{dwplot(df,
        dot_args = list(size = 2),
        whisker_args = list(size = 1),
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(Sugar = "Sugar", 
                         Maize = "Maize",
                         Crude_oil = "Crude oil",
                         Palm_oil = "Palm oil",
                         Soybean_oil = "Soybean oil",                       
                         Olive_oil = "Olive oil", 
                         Rapeseed_oil = "Rapeseed oil", 
                         Sunflower_oil = "Sunflower oil", 
                         Coconut_oil = "Coconut oil", 
                         Soybean = "Soybean", 
                         Soybean_meal = "Soybean meal", 
                         Cocoa = "Cocoa", 
                         Coffee = "Coffee", 
                         Tea = "Tea", 
                         Tobacco = "Tobacco"
    )) +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    ggtitle(title) +
    theme(plot.title = element_text(face="bold", size=c(10)),
          legend.position = c(0.007, 0.001),
          legend.justification = c(0, 0), 
          legend.background = element_rect(colour="grey80"),
          legend.title = element_blank())}  




#### COMPARE COEFFICIENTS #### 
direct_effects_var <- paste0(mapmat[,"Crops"], "_X_", mapmat[,"Prices"]) 

# Net effects
df_res_total <- make_main_reg(output = "coeff_table")

## Effects including effects from supply substitution
# With Soybean
df_res_wo_Fodder <- make_main_reg(open_path = TRUE,
                                   commoXcommo = "Fodder_X_Beef",
                                   output = "coeff_table")

# With Soybean
df_res_wo_Soybean <- make_main_reg(open_path = TRUE,
                                   commoXcommo = "Soybean_X_Soybean",
                                    output = "coeff_table")

# With Oilpalm
df_res_wo_Oilpalm <- make_main_reg(open_path = TRUE,
                                   commoXcommo = "Oilpalm_X_Palm_oil",
                                   output = "coeff_table")


compare_APEs_across_groups <- function(full_coeffs, net_coeffs, m0 = 0, alternative = "two.sided", 
                                       # this argument identifies the commodity k for which we want to compute the indirect
                                       # effect through price shocks 
                                       # Give in GAEZ spelling
                                       j_removed = "Fodder"
                                       ) { 
  # important that it's selecting the "1" (resp. 2) row, and not the "Estimate" (resp. "SE") row in cases where there are several 
  # rows with the same names, i.e. if APEs of interaction effects have been computed as well. 
  
  # full_coeffs always has less regressors, since net_coeffs is the output from specification with all possible interactions.
  
  # define coefficients for which we make the comparison 
  coefficients <- rownames(full_coeffs)[grepl(j_removed, rownames(full_coeffs))]
  for(reg in coefficients){
  
    coeff1 <- full_coeffs[reg,"estimate"]
    coeff2 <- net_coeffs[reg,"estimate"]
    sigma1 <- full_coeffs[reg,"std.error"]
    sigma2 <- net_coeffs[reg,"std.error"]
    
    statistic <- (ape1 - ape2 - m0) / sqrt(sigma1 + sigma2)
    
    pval <- if (alternative == "two.sided") { 
      2 * pnorm(abs(statistic), lower.tail = FALSE) 
    } else if (alternative == "less") { 
      pnorm(statistic, lower.tail = TRUE) 
    } else { 
      pnorm(statistic, lower.tail = FALSE) 
    } 
    # LCL <- (M1 - M2 - S * qnorm(1 - alpha / 2)) UCL <- (M1 - M2 + S * qnorm(1 - alpha / 2)) value <- list(mean1 = M1, mean2 = M2, m0 = m0, sigma1 = sigma1, sigma2 = sigma2, S = S, statistic = statistic, p.value = p, LCL = LCL, UCL = UCL, alternative = alternative) 
    # print(sprintf("P-value = %g",p)) # print(sprintf("Lower %.2f%% Confidence Limit = %g", 
    # alpha, LCL)) # print(sprintf("Upper %.2f%% Confidence Limit = %g", # alpha, UCL)) return(value) } test <- t.test_knownvar(dat1$sample1, dat1$sample2, V1 = 1, V2 = 1 )
    return(pval)
  }
}

