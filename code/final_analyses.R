
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
                   "ggplot2", "dotwhisker"# "leaflet", "htmltools"
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
                 first_loss = "First forest loss",
                 Soybean_meal = "Soybean meal",
                 Soybean_oil = "Soybean oil",
                 Olive_oil = "Olive oil",
                 Rapeseed_oil = "Rapeseed oil",
                 Sunflower_oil = "Sunflower oil",
                 Coconut_oil = "Coconut oil"))

### This matrix is used for maping crops from GAEZ with commodities from price data sets
mapmat_data <- c(
"Banana","Banana",
"Barley", "Barley",
"Beef", "fodder_crops", 
"Orange", "Citrus",
"Cocoa", "Cocoa",
"Coconut_oil", "Coconut",
"Coffee", "Coffee",
"Cotton", "fibre_crops",
"Groundnuts", "Groundnut",
"Maize", "Maize",
"Oat", "Oat",
"Olive_oil", "Olive",
"Palm_oil", "Oilpalm",
"Rapeseed_oil", "Rapeseed",
"Rice", "rice_crops",
"Sorghum", "Sorghum",
"Soybeans", "Soybean",
"Soybean_oil", "Soybean",
"Soybean_meal", "Soybean",
"Sugar", "sugar_crops",
"Sunflower_oil", "Sunflower",
"Tea", "Tea",
"Tobacco", "Tobacco",
"Wheat", "Wheat",
"cereal_crops", "cereal_crops")

mapmat <- matrix(data = mapmat_data, 
                 nrow = length(mapmat_data)/2,
                 ncol = 2, 
                 byrow = TRUE)

colnames(mapmat) <- c("Prices", "Crops")


### READ ALL POSSIBLE DATASETS HERE 
# but not all together because of memory issues. 
# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aesi_long_final.Rdata"))

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_aesi_long_final.Rdata"))

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aesi_long_final.Rdata"))

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_acay_long_final.Rdata"))

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_acay_long_final.Rdata"))

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_acay_long_final.Rdata"))

prices <- readRDS(here("temp_data", "prepared_prices.Rdata"))

rm(outcome_variable, start_year, end_year, crop_j, j_soy, price_k, extra_price_k, SjPj, SkPk, fe, distribution, output, se, cluster, 
   controls, regressors, outcome_variable)

outcome_variable = "first_loss" # "first_loss", "firstloss_glassgfc", "phtf_loss"
start_year = 2002
end_year = 2020
price_info = "lag1"
further_lu_evidence = "none"
crop_j = "Oilpalm"
j_soy = "Soybean_oil"
# for the k variables hypothesized in overleaf for palm oil, feglm quasipoisson converges within 25 iter. 
# but maybe not with skPk controls. 
price_k <- c("Soybean_oil", "Rapeseed_oil", "Sunflower_oil", 
             "Coconut_oil", "Soybeans","Sugar", "Maize")
# c( "Sugar", "Maize") 
# "Barley",  "Chicken", "Sheep", "Banana", "Beef", "Olive_oil",
#               "Orange",  "Cotton",  "Groundnuts",  "Rubber", "Sorghum","Cocoa",  "Coffee",
#                "Rice",   "Wheat",  "Palm_oil", ), 
#                "Tea", "Tobacco",  "Oat",  
#               , "Pork")
extra_price_k = c() # , "Sheep", "Pork", "Chicken"
SjPj = TRUE
SkPk = TRUE
fe = "grid_id + country_year"
distribution = "quasipoisson"
se = "cluster"
cluster ="grid_id"
coefstat = "confint"
output = "coef_table"

# "Banana", "Barley", "Beef", 
# "Orange", "Cocoa", "Coconut_oil", "Coffee", "Cotton", "Rice", "Groundnuts", 
# "Maize", "Palm_oil", "Rubber", "Sorghum", "Soybean_oil", 
# "Sugar", "Tea", "Tobacco", "Wheat", "Oat", "Olive_oil", "Rapeseed_oil", 
# "Sunflower_oil"
# rm(dataset, start_year, end_year, crop_j, j_soy, price_k, SjPj, SkPk, fe, distribution, output, se, cluster,     controls, regressors, outcome_variable)

make_reg_acay <- function(outcome_variable = "first_loss", # one of "first_loss", "firstloss_glassgfc", "phtf_loss"
                          start_year = 2002, 
                          end_year = 2020, 
                          further_lu_evidence = "none", # either "none", "sbqt_direct_lu", or "sbqt"mode_lu"
                          crop_j = "Oilpalm", # in GAEZ spelling
                          j_soy = "Soybean_oil", # in case crop_j is Soybean, which price should be used? One of "Soybeans", "Soybean_oil", "Soybean_meal".
                          price_k = c("Sugar", "Maize"), # in prices spelling
                          extra_price_k = c(), # One of "Crude_oil", "Chicken", "Pork", "Sheep" 
                          price_info = "lag1", # one of "lag1", "2pya", "3pya", "4pya", "5pya",
                          SjPj = TRUE,
                          SkPk = FALSE,
                          #commo_m = c(""), comment coder ça pour compatibilité avec loops over K_commo ? 
                          fe = "grid_id + country_year", 
                          distribution = "gaussian",#  "quasipoisson", 
                          se = "cluster", 
                          cluster ="grid_id",
                          # coefstat = "confint", # one of "se", "tstat", "confint"
                          output = "coef_table" # one of "data", "est_obj", "coef_table" 
){
  
  
  #### PREPARE NEEDED VARIABLE NAMES
  # this does not involve data, just arguments of the make_reg function 
  
  ## Price variable names
  # Explicit the name of the price of j based on crop_j (GAEZ spelling)
  price_j <- mapmat[mapmat[,"Crops"]==crop_j,"Prices"]
  names(price_j) <- NULL
  
  # Group the names of the different prices 
  original_price_all <- c(price_j, price_k, extra_price_k)
  original_price_k <- price_k
  # lag and log
  price_all <- paste0(original_price_all, "_", price_info)
  price_all <- paste0("ln_", price_all)
  price_k <- paste0(original_price_k, "_", price_info)
  price_k <- paste0("ln_", price_k)
  
  ## Revenue variable names
  # To determine land use, only the potential revenue of soybeans is considered (for simplicity) 
  all_crop_prices <- mapmat[mapmat[,"Prices"]!="Soybean_oil" & mapmat[,"Prices"]!="Soybean_meal" ,"Prices"]
  all_crop_acay <- mapmat[mapmat[,"Prices"]!="Soybean_oil" & mapmat[,"Prices"]!="Soybean_meal" ,"Crops"]
  
  # We will construct expected revenues for all crops - not just those in price_k or crop_j.
  # Thus we take the first lag or the X past year average of the prices, but not the log.
  prices_4_revenue <- paste0(all_crop_prices, "_", price_info)
  
  # Explicit the name of the variable rj, the standardized potential revenue of j. 
  revenue_j <- paste0("R_", crop_j)
  revenue_j_std <- paste0(revenue_j, "_std")
  
  # Explicit the names of the rk variables corresponding to commodities in price_k only (i.e. not in price_j nor extra_price_k)
  # match function is important because revenue_k_std must have commodities in the same order as price_k 
  revenue_k_std <- mapmat[,"Crops"][match(original_price_k, mapmat[,"Prices"])]
  revenue_k_std <- paste0("R_", revenue_k_std, "_std")
  

  #### MAKE THE VARIABLES NEEDED IN THE DATA
  #d <- main_data
  if(outcome_variable == "first_loss"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_acay_long_final.Rdata"))}
  if(outcome_variable == "firstloss_glassgfc"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_acay_long_final.Rdata"))}
  if(outcome_variable == "phtf_loss"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_acay_long_final.Rdata"))}
  
  # # restrict the outcome variable to evidence of further LU
  # d <- dplyr::mutate(d, lu_evidence = TRUE)
  # if(outcome_variable == "first_loss" & further_lu_evidence != "none"){
  #    # , cropland, or tree plantation
  #   # if it is actually forest afterwards, this means that 
  #   
  #   if(crop_j == "fodder_crops"){
  #     d <- dplyr::mutate(d, lu_evidence = (!!as.symbol(further_lu_evidence) == 30)) # either direct or mode subsequent lu are grassland
  #   }
  #   if(crop_j == "Soybean"){
  #     d <- dplyr::mutate(d, lu_evidence = (!!as.symbol(further_lu_evidence) == 10)) # either direct or mode subsequent lu are cropland
  #   }
  #   if(crop_j %in% c("Oilpalm", "Cocoa", "Coffee")){
  #     d <- dplyr::mutate(d, lu_evidence = (sbqt_direct_lu != 20 & sbqt_mode_lu == 20)) # is not forest in the year directly after, but is forest again in the mode subsequent lu. 
  #   }
  # }
  
  
  ### PREPARE rj, the standardized achievable revenues
  
  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", prices_4_revenue)], by = "year")
  
  # Agro-climatic achievable yield data from GAEZ are in t/ha (as shown by the Map tool on the platform)
  # Prices have been converted to $/t in prepare_prices.R
  for(i in 1:length(prices_4_revenue)){
    price_i <- prices_4_revenue[i]
    revenue_i <- paste0("R_", all_crop_acay[i])
    acay_i <- all_crop_acay[i]
    d <- dplyr::mutate(d, 
                       !!as.symbol(revenue_i) := !!as.symbol(acay_i)*!!as.symbol(price_i))
  }
  rm(price_i, revenue_i, acay_i)                   

  # scale fodder crop revenues by a feed conversion ratio. 
  # Following Alexander et al. 2016, 25 tons of feed transform into 1 ton of beef meat (in edible weight). 
  # Thus, every ton of feed counted in R_fodder_crops (coming from agro-climatically achievable yields in ton/ha) is scaled to 1/25 ton of beef meat  
  d <- dplyr::mutate(d, R_fodder_crops = R_fodder_crops/25)
  
  
  ## Standardize the revenue variable(s) needed here -> stream dj = Rj/sumRi 
  
  # To understand this line, see https://dplyr.tidyverse.org/articles/rowwise.html#row-wise-summary-functions
  d <- dplyr::mutate(d, Ri_sum = rowSums(across(.cols = starts_with("R_", ignore.case = FALSE))))
  
  # if we control for SkPk, then we need the standardized revenues of many crops
  if(SkPk){
    # divide each Revenue variable by the sum of them, to standardize.
    d <- dplyr::mutate(d, across(.cols = starts_with("R_", ignore.case = FALSE),
                                 .fns = ~./Ri_sum, 
                                 .names = paste0("{.col}", "_std"))) 
  }else{
    # divide only Revenue of j
    d <- dplyr::mutate(d, !!as.symbol(revenue_j_std) := !!as.symbol(revenue_j)/Ri_sum)
  }
  
  ## Construct the max indicatrice -> stream dj = 1[Rj = maxRi]
  # d <- dplyr::mutate(d, across(.cols = starts_with("R_", ignore.case = FALSE),
  #                              .names = paste0("{.col}", "_max")))

  # Keep only usefull variables
  vars <- c("grid_id", "year", "lon", "lat", "country_name", "country_year", outcome_variable, 
              names(d)[grepl(pattern ="_std", x= names(d))])
  d <- dplyr::select(d, all_of(vars))


  ### MAKE FINAL REGRESSION VARIABLES 
  
  # keep only the cells with positive rj (since we need to divide by rj)
  d <- dplyr::filter(d, !!as.symbol(revenue_j_std) > 0)
  
  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", price_all)], by = "year")
  
  # Main regressors
  regressors <- c()
  for(Pk in price_all){
    varname <- paste0(crop_j, "_X_", original_price_all[match(Pk, price_all)])
    regressors <- c(regressors, varname)
    d <- mutate(d, 
                !!as.symbol(varname) := !!as.symbol(Pk)/(!!as.symbol(revenue_j_std)))
  }
  rm(varname)
  
  # Remove the first of the regressors if we do not control for SjPj. The first is necessart SjPj because of the construction of price_all that has price_j first. 
  if(!SjPj){
    regressors <- regressors[-1]
  }
  
  # Controls  
  controls <- c() # it's important that this is not conditioned on SkPk
  if(SkPk){
    for(Pk in price_k){
      rk <- revenue_k_std[match(Pk, price_k)]
      varname <- paste0("ctrl_", original_price_k[match(Pk, price_k)])
      controls <- c(controls, varname)
      d <- mutate(d, 
                  !!as.symbol(varname) := !!as.symbol(rk)*!!as.symbol(Pk))
    }
  }
  
  
  ### MODEL SPECIFICATION FORMULAE
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
  
  # # - have a lu_evidence 
  # d <- dplyr::filter(d, lu_evidence)
  
  used_vars <- c("grid_id", "year", "country_name", "country_year", 
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
  
  # is.na(d$Oilpalm_std) %>% sum()
  
  
  # - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove
  d_clean <- d[-obs2remove(fml = as.formula(paste0(outcome_variable, " ~ ", fe)),
                           d, 
                           family = "poisson"),]
  rm(d)
  
  ### REGRESSIONS
  
  if(distribution != "negbin"){ # i.e. if it's poisson or quasipoisson or gaussian
    reg_res <- fixest::feglm(fe_model,
                             data = d_clean, 
                             family = "gaussian",#distribution,#   # "poisson" ,
                             # glm.iter = 25,
                             #fixef.iter = 100000,
                             nthreads = 3,
                             notes = TRUE, 
                             verbose = 4)  
    
    
  }
  
  # Now keep only information necessary, otherwise the output of fixest estimation is large and we can't collect too many at the same time (over loops)  
  # this is necessary to compute SE as we want to.  
  df_res <- summary(reg_res, se = se,
                    cluster = cluster)$coeftable
  
  # ci <- confint(reg_res, se =se, cluster = cluster, level = 0.95)   
  # 
  # df_res <- cbind(df_res, ci)
  
  # add a column with the number of observations
  df_res$Observations <- reg_res$nobs
  # mat[row.names(mat)=="Observations",] <- mat[row.names(mat)=="Observations",] %>% formatC(digits = 0, format = "f")
  
  # add a row with the number of clusters
  df_res$Clusters <- length(unique(d_clean[,cluster]))
  # mat[row.names(mat)=="Clusters",] <- mat[row.names(mat)=="Clusters",] %>% formatC(digits = 0, format = "f")
  
  
  if(output == "data"){
    toreturn <- list(reg_res, d_clean)
  }
  if(output == "est_obj"){
    toreturn <- reg_res
  }
  if(output == "coef_table"){
    toreturn <- df_res
  }
  
  
  rm(d_clean, reg_res, df_res)
  return(toreturn)
  rm(toreturn)

}


make_reg_aesi <- function(outcome_variable = "first_loss", # one of "first_loss", "firstloss_glassgfc", "phtf_loss"
                          start_year = 2002, 
                          end_year = 2020, 
                          further_lu_evidence = "none", # either "none", "sbqt_direct_lu", or "sbqt"mode_lu"
                          crop_j = "Oilpalm", # in GAEZ spelling
                          j_soy = "Soybean_oil", # in case crop_j is Soybean, which price should be used? One of "Soybeans", "Soybean_oil", "Soybean_meal".
                          price_k = c("Sugar", "Maize"), # in prices spelling
                          extra_price_k = c(), # One of "Crude_oil", "Chicken", "Pork", "Sheep" 
                          price_info = "lag1", # one of "lag1", "2pya", "3pya", "4pya", "5pya",
                          SjPj = TRUE,
                          SkPk = FALSE,
                          #commo_m = c(""), comment coder ça pour compatibilité avec loops over K_commo ? 
                          fe = "grid_id + country_year", 
                          distribution = "gaussian",#  "quasipoisson", 
                          se = "cluster", 
                          cluster ="grid_id",
                          # coefstat = "confint", # one of "se", "tstat", "confint"
                          output = "coef_table" # one of "data", "est_obj", "coef_table" 
                          ){

  #### PREPARE NEEDED VARIABLE NAMES
  # this does not involve data, just arguments of the make_reg function 

  
  ## Price variable names
  # Explicit the name of the price of j based on crop_j (GAEZ spelling)
  price_j <- mapmat[mapmat[,"Crops"]==crop_j,"Prices"]
  names(price_j) <- NULL
  
  # Group the names of the different prices 
  original_price_all <- c(price_j, price_k, extra_price_k)
  original_price_k <- price_k
  # lag and log
  price_all <- paste0(original_price_all, "_", price_info)
  price_all <- paste0("ln_", price_all)
  price_k <- paste0(original_price_k, "_", price_info)
  price_k <- paste0("ln_", price_k)

  ## Suitability index variable names

  # Explicit the name of the variable sj, the standardized suitability index of j. 
  suitability_j_std <- paste0(crop_j, "_std")  
  
  # Explicit the names of the sk variables corresponding to commodities in price_k only (i.e. not in price_j nor extra_price_k)
  # match function is important because suitability_k_std must have commodities in the same order as price_k 
  crop_k <- mapmat[,"Crops"][match(original_price_k, mapmat[,"Prices"])]
  suitability_k_std <- paste0(crop_k, "_std")
  
  
  
  
  # ### SELECT PRICE AND SUITABILITY VARIABLES 
  # # first, save price_k for later purpose 
  # original_price_j <- price_j 
  # original_price_k <- price_k 
  # original_extra_price_k <- extra_price_k 
  # original_price_all <- c(price_j, price_k, extra_price_k)
  # 
  # # For j
  # if(price_lag != 0){price_j <- paste0(price_j,"_lag",price_lag)}
  # price_j <- paste0("ln_", price_j)
  # 
  # # For k 
  # if(price_lag != 0){price_k <- paste0(price_k,"_lag",price_lag)}
  # # don't condition that, as we won't use non logged prices a priori. 
  # price_k <- paste0("ln_", price_k)
  # 
  # # For extra commodities k 
  # if(length(extra_price_k) > 0){
  #   if(price_lag != 0){extra_price_k <- paste0(extra_price_k,"_lag",price_lag)}
  #   extra_price_k <- paste0("ln_", extra_price_k)
  # }
  
  # # group all prices transformed
  # price_all <- c(price_j, price_k, extra_price_k)
  

  
  
  #### MAKE THE VARIABLES NEEDED IN THE DATA
  #d <- main_data
  if(outcome_variable == "first_loss"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aesi_long_final.Rdata"))}
  if(outcome_variable == "firstloss_glassgfc"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_aesi_long_final.Rdata"))}
  if(outcome_variable == "phtf_loss"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aesi_long_final.Rdata"))}
  
  # Keep only in data the useful variables 
  d <- dplyr::select(d, all_of(c("grid_id", "year", "country_name", "country_year", outcome_variable, 
                               suitability_j_std, suitability_k_std))) 
    
  # and keep only the cells with positive suitability for crop j 
  d <- dplyr::filter(d, !!as.symbol(suitability_j_std) > 0)

  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", price_all)], by = "year")
  
  
  # Main regressors
  regressors <- c()
  for(Pk in price_all){
    varname <- paste0(crop_j, "_X_", original_price_all[match(Pk, price_all)])
    regressors <- c(regressors, varname)
    d <- mutate(d, 
                !!as.symbol(varname) := !!as.symbol(Pk)/(!!as.symbol(suitability_j_std)))
  }
  rm(varname)
  
  # Remove the first of the regressors if we do not control for SjPj. The first is necessart SjPj because of the construction of price_all that has price_j first. 
  if(!SjPj){
    regressors <- regressors[-1]
  }
  
  # Controls  
  controls <- c() # it's important that this is not conditioned on SkPk
  if(SkPk){
    for(Pk in price_k){
      sk <- suitability_k_std[match(Pk, price_k)]
      varname <- paste0("ctrl_", original_price_k[match(Pk, price_k)])
      controls <- c(controls, varname)
      d <- mutate(d, 
                  !!as.symbol(varname) := !!as.symbol(sk)*!!as.symbol(Pk))
    }
  }
  rm(varname)

  
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
  
  
  used_vars <- c("grid_id", "year", "country_name", "country_year", 
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

  # is.na(d$Oilpalm_std) %>% sum()
  
  
  # - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove
  d_clean <- d[-obs2remove(fml = as.formula(paste0(outcome_variable, " ~ ", fe)),
                                d, 
                                family = "poisson"),]
  
  ### REGRESSIONS
  
  if(distribution != "negbin"){ # i.e. if it's poisson or quasipoisson or gaussian
    reg_res <- fixest::feglm(fe_model,
                             data = d_clean, 
                             family = distribution,# "gaussian",#  # "poisson" ,
                             # glm.iter = 25,
                             #fixef.iter = 100000,
                             nthreads = 3,
                             notes = TRUE, 
                             verbose = 4)  

   
  }
  
  # Now keep only information necessary, otherwise the output of fixest estimation is large and we can't collect too many at the same time (over loops)  
  # this is necessary to compute SE as we want to.  
  df_res <- summary(reg_res, se = se,
                    cluster = cluster)$coeftable
  
  # ci <- confint(reg_res, se =se, cluster = cluster, level = 0.95)   
  # 
  # df_res <- cbind(df_res, ci)
  
  # add a column with the number of observations
  df_res$Observations <- reg_res$nobs
  # mat[row.names(mat)=="Observations",] <- mat[row.names(mat)=="Observations",] %>% formatC(digits = 0, format = "f")
  
  # add a row with the number of clusters
  df_res$Clusters <- length(unique(d_clean[,cluster]))
  # mat[row.names(mat)=="Clusters",] <- mat[row.names(mat)=="Clusters",] %>% formatC(digits = 0, format = "f")
  
  
  if(output == "data"){
    toreturn <- list(reg_res, d_clean)
  }
  if(output == "est_obj"){
    toreturn <- reg_res
  }
  if(output == "coef_table"){
    toreturn <- df_res
  }

  
  rm(d, d_clean, reg_res, df_res)
  return(toreturn)
  rm(toreturn)
}


# This function takes as main input the table of results outputted from one regression. It's purpose is to be applied over a list of such results. 
make_table_mat <- function(df_res, 
                           rounding=4
                           ){
  
  # Note if there are SkPk controls
  SkPk <- any(grepl(pattern = "ctrl_", x = row.names(df_res)))
  # Keep only variable of interest
  df_res <- df_res[!grepl(pattern = "ctrl_", x = row.names(df_res)),]  
  
  df_res <- round(df_res, digits = rounding)
  
  
  # make a one column matrix with all coeficient point estimates and std error next to each others
  # for each coefficient we want the point estimate and the std error (hence *2). In addition we want # obs. and # clusters (hence +2)
  
  mat_res <- NULL
  for(i in 1:nrow(df_res)){
    mat_res <- rbind(mat_res, as.matrix(df_res[i,"Estimate"]), as.matrix(paste0("(",df_res[i,"Std. Error"],")")))
  }
  
  # make row names
  k_crops <- row.names(df_res)
  k_crops <- sub(pattern = ".+?(_X_)", x = k_crops, replacement = "") # replace everything before the first underscore with nothing
  names <- c()
  for(crop in k_crops){
    names <- c(names, crop, "")
  }
  row.names(mat_res) <- names 
  
  # Mark whether there were controls or not
  if(SkPk){
    mat_res <- rbind(mat_res, "X")
  } else {
    mat_res <- rbind(mat_res, "")
  }
  row.names(mat_res)[nrow(mat_res)] <- "Controls"
  
  # Indicate the number of observations and clusters
  mat_res <- rbind(mat_res, unique(df_res$Observations))
  row.names(mat_res)[nrow(mat_res)] <- "Observations"

  mat_res <- rbind(mat_res, unique(df_res$Clusters))
  row.names(mat_res)[nrow(mat_res)] <- "Clusters"
  return(mat_res)
  
}


#### DESCRIPTIVE STATISTICS #### 
# make_des_stats <- function(outcome_variable, 
#                            crop_j,
#                            start_year = 2002, 
#                            end_year = 2020, 
#                            qj = "rj", 
#                            price_info = "lag1")


# For each crop j, we produce two quantities, for aesi and acay. 
# Across crops, these quantities differ in the forest loss measure, and (the time and) space it is aggregated to. 

# For acay we need to compute the revenue here. This helper function does it
make_rj <- function(data, crop_j, price_info){
  ## Revenue variable names
  # To determine land use, only the potential revenue of soybeans is considered (for simplicity) 
  all_crop_prices <- mapmat[mapmat[,"Prices"]!="Soybean_oil" & mapmat[,"Prices"]!="Soybean_meal" ,"Prices"]
  all_crop_acay <- mapmat[mapmat[,"Prices"]!="Soybean_oil" & mapmat[,"Prices"]!="Soybean_meal" ,"Crops"]
  
  # We will construct expected revenues for all crops - not just those in price_k or crop_j.
  # Thus we take the first lag or the X past year average of the prices, but not the log.
  prices_4_revenue <- paste0(all_crop_prices, "_", price_info)
  
  # Explicit the name of the variable rj, the standardized potential revenue of j. 
  revenue_j <- paste0("R_", crop_j)
  revenue_j_std <- paste0(revenue_j, "_std")
  

  ### MAKE THE VARIABLES NEEDED IN THE DATA
  
  ### PREPARE rj, the standardized achievable revenues
  
  # Merge only the prices needed, not the whole price dataframe
  data <- left_join(data, prices[,c("year", prices_4_revenue)], by = "year")
  
  # Agro-climatic achievable yield data from GAEZ are in t/ha (as shown by the Map tool on the platform)
  # Prices have been converted to $/t in prepare_prices.R
  for(i in 1:length(prices_4_revenue)){
    price_i <- prices_4_revenue[i]
    revenue_i <- paste0("R_", all_crop_acay[i])
    acay_i <- all_crop_acay[i]
    data <- dplyr::mutate(data, 
                       !!as.symbol(revenue_i) := !!as.symbol(acay_i)*!!as.symbol(price_i))
  }
  rm(price_i, revenue_i, acay_i)                   
  
  # scale fodder crop revenues by a feed conversion ratio. 
  # Following Alexander et al. 2016, 25 tons of feed transform into 1 ton of beef meat (in edible weight). 
  # Thus, every ton of feed counted in R_fodder_crops (coming from agro-climatically achievable yields in ton/ha) is scaled to 1/25 ton of beef meat  
  data <- dplyr::mutate(data, R_fodder_crops = R_fodder_crops/25)
  
  
  ## Standardize the revenue variable(s) needed here -> stream dj = Rj/sumRi 
  
  # To understand this line, see https://dplyr.tidyverse.org/articles/rowwise.html#row-wise-summary-functions
  data <- dplyr::mutate(data, Ri_sum = rowSums(across(.cols = starts_with("R_", ignore.case = FALSE))))
  
  data <- dplyr::mutate(data, !!as.symbol(revenue_j_std) := !!as.symbol(revenue_j)/Ri_sum)
  
  ## Construct the max indicatrice -> stream dj = 1[Rj = maxRi]
  # data <- dplyr::mutate(data, across(.cols = starts_with("R_", ignore.case = FALSE),
  #                              .names = paste0("{.col}", "_max")))
  
  return(data)
  
}

# Infrastructure to store results
defo_table <- matrix(nrow = 5, ncol = 5, data = "")
colnames(defo_table) <- c("Crop", "Region", "Goldman et al. 2020 estimate", "AESI estimate", "ACAY estimate")
row.names(defo_table) <- c("Pasture", "Oil palm", "Soybean", "Cocoa", "Coffee")

# First we need to prepare some shapes for Soy and Cattle
countries <- st_read(here("input_data", "Global_LSIB_Polygons_Detailed"))

# We simplify the continent but not Brazil, as the border effect on deforestation is significant (Burgess) 
brazil <- countries[countries$COUNTRY_NA=="Brazil",] %>% st_geometry() %>% 
  st_transform("EPSG:31970","+proj=utm +zone=16 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") 

southam <- countries[countries$COUNTRY_NA=="Argentina" |
                       countries$COUNTRY_NA=="Bolivia" |
                       countries$COUNTRY_NA=="Brazil" | countries$COUNTRY_NA=="Isla Brasilera (disp)" | 
                       countries$COUNTRY_NA=="Chile" |   
                       countries$COUNTRY_NA=="Colombia" | 
                       countries$COUNTRY_NA=="Ecuador" | 
                       countries$COUNTRY_NA=="French Guiana (Fr)" |
                       countries$COUNTRY_NA=="Guyana" |
                       countries$COUNTRY_NA=="Panama" |
                       countries$COUNTRY_NA=="Paraguay" | 
                       countries$COUNTRY_NA=="Peru" | 
                       countries$COUNTRY_NA=="Suriname" | 
                       countries$COUNTRY_NA=="Trinidad & Tobago" | 
                       countries$COUNTRY_NA=="Uruguay" | 
                       countries$COUNTRY_NA=="Venezuela"  ,] %>% st_union() %>% st_geometry() %>% 
  st_transform("EPSG:31970","+proj=utm +zone=16 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") %>% 
  st_simplify(dTolerance = 10000) # simplifies by 10km


### OIL PALM 
## AESI
# Summarize at global (tropical) scale, first loss,in 2001-2015. 
d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aesi_long_final.Rdata"))
d <- dplyr::select(d, grid_id, year, country_name, first_loss, Oilpalm_std)
d <- dplyr::mutate(d, Y_oilpalm = first_loss*Oilpalm_std)
d <- dplyr::filter(d, year >= 2001 & year <= 2015)
summary(d$first_loss)
summary(d$Y_oilpalm)

### SOY
# For soy, we input different data, prepared in south_america_track.R
# This is necessary because we want to compare to a measure of deforestation to soy in Goldman et al. 2020 that aggregates the whole continent.

## AESI
d <- readRDS(here("temp_data", "merged_datasets", "southam_aoi", "glass_aesi_long_final.Rdata"))
d <- dplyr::select(d, grid_id, year, lon, lat, first_loss, Soybean_std)
d <- dplyr::mutate(d, Y_soybean = first_loss*Soybean_std)

# add up annual deforestation over the period
d <- dplyr::filter(d, year >= 2001 & year <= 2015)

accu <- ddply(d, "grid_id", summarise, accu_defo = sum(Y_soybean, na.rm = TRUE))
d_cs <- d[!duplicated(d$grid_id), c("grid_id", "year", "lon", "lat")]
accu <- left_join(accu, d_cs, by = "grid_id")


# Restrict to cells in South America 
accu <- st_as_sf(accu, coords=c("lon", "lat"), crs = 4326)
accu <- st_transform(accu, "EPSG:31970","+proj=utm +zone=16 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") 

plot(st_geometry(accu))

sgbp <- st_within(accu, southam)

accu_southam <- accu[lengths(sgbp)>0,]

defo_table["Soybean", "AESI estimate"] <- round(sum(accu_southam$accu_defo)/1e6,digits=2) 

## ACAY
d <- readRDS(here("temp_data", "merged_datasets", "southam_aoi", "glass_acay_long_final.Rdata"))
# here we need to construct rj 
d <- make_rj(data = d, crop_j = "Soybean", price_info = "lag1")
d <- dplyr::select(d, grid_id, year, lon, lat, first_loss, R_Soybean_std)
d <- dplyr::mutate(d, Y_soybean = first_loss*R_Soybean_std)

# add up annual deforestation over the period
d <- dplyr::filter(d, year >= 2001 & year <= 2015)
accu <- ddply(d, "grid_id", summarise, accu_defo = sum(Y_soybean, na.rm = TRUE))
d_cs <- d[!duplicated(d$grid_id), c("grid_id", "lon", "lat")]
accu <- left_join(accu, d_cs, by = "grid_id")


# Restrict to cells in South America 
# not necessary to compute sgbp again (which is long), just reuse that of aesi
# (Even if there are missing in rj that were not in sj, there is no filtering based on this.)

# accu <- st_as_sf(accu, coords=c("lon", "lat"), crs = 4326)
# accu <- st_transform(accu, "EPSG:31970","+proj=utm +zone=16 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") 

# sgbp <- st_within(accu, southam)

accu_southam <- accu[lengths(sgbp)>0,]

defo_table["Soybean", "ACAY estimate"] <- round(sum(accu_southam$accu_defo)/1e6,digits=2) 



### COCOA
d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aesi_long_final.Rdata"))
d <- dplyr::select(d, grid_id, year, country_name, firstloss_glassgfc, Cocoa_std)
d <- dplyr::mutate(d, Y_Cocoa = firstloss_glassgfc*Cocoa_std)

### COFFEE
d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aesi_long_final.Rdata"))
d <- dplyr::select(d, grid_id, year, country_name, firstloss_glassgfc, Coffee_std)
d <- dplyr::mutate(d, Y_coffee = firstloss_glassgfc*Coffee_std)
















### Pasture and Soy are compared in out-of-tropical-boundary areas.  

## Limit to 2001-2015 and sum over these years
first_loss <- readRDS(here("temp_data", "merged_datasets", "southam_aoi", "glass_aesi_long_final.Rdata"))

# select 2001-2015 period 
first_loss <- dplyr::filter(first_loss, year >= 2001 & year <= 2015)

# Add up annual aggregated LUCFP (result is a single layer with cell values = the sum of annual cell values over the selected time period)
accu_firstloss <- calc(first_loss, fun = sum, na.rm = TRUE)

## Extract in Brazil and in all South America

countries <- st_read(here("input_data", "Global_LSIB_Polygons_Detailed"))

brazil <- countries[countries$COUNTRY_NA=="Brazil",] %>% st_geometry() %>% as("Spatial")

southam <- countries[countries$COUNTRY_NA=="Argentina" |
                       countries$COUNTRY_NA=="Bolivia" |
                       countries$COUNTRY_NA=="Brazil" | countries$COUNTRY_NA=="Isla Brasilera (disp)" | 
                       countries$COUNTRY_NA=="Chile" |   
                       countries$COUNTRY_NA=="Colombia" | 
                       countries$COUNTRY_NA=="Ecuador" | 
                       countries$COUNTRY_NA=="French Guiana (Fr)" |
                       countries$COUNTRY_NA=="Guyana" |
                       countries$COUNTRY_NA=="Panama" |
                       countries$COUNTRY_NA=="Paraguay" | 
                       countries$COUNTRY_NA=="Peru" | 
                       countries$COUNTRY_NA=="Suriname" | 
                       countries$COUNTRY_NA=="Trinidad & Tobago" | 
                       countries$COUNTRY_NA=="Uruguay" | 
                       countries$COUNTRY_NA=="Venezuela"  ,] %>% st_union() %>% st_geometry() %>% as("Spatial")

firstloss_brazil <- raster::extract(x = accu_firstloss, 
                                    y =  brazil,
                                    fun = sum, 
                                    na.rm = TRUE) 

firstloss_southam <- raster::extract(x = accu_firstloss, 
                                     y =  southam,
                                     fun = sum, 
                                     na.rm = TRUE) 



#### CONSTRUCT SPECIFICATION COMPARATIVE TABLES ####

# We build 5 tables, one for each j crop. 
# These tables compare first_loss and phtf_loss, for the same period, 2002, 2020, for acay (rj) and aesi (sj) models, with and without controlling for SkPk. 

### BEEF
K_beef <- c("Soybeans", "Soybean_meal", "Palm_oil", 
            "cereal_crops") # in prices spelling "Barley", "Oat", "Sorghum", "Maize", "Rice"

res_list_beef <- list()
elm <- 1

# Forest definition
outcome_variableS <- c("firstloss_glassgfc", "phtf_loss")

# with controls or not 
control_ornot <- c(FALSE, TRUE)


for(OV in outcome_variableS){
  for(CTRL in control_ornot){
    res_list_beef[[elm]] <- make_reg_acay(outcome_variable = OV,
                                               start_year = 2002, end_year = 2020,
                                               crop_j = "fodder_crops",
                                               price_k = K_beef,
                                               SkPk= CTRL,
                                               extra_price_k = c())#"Sheep", "Pork", "Chicken"
    names(res_list_beef)[elm] <- paste0(OV,"_",CTRL, "_acay")
    elm <- elm + 1
    
    res_list_beef[[elm]] <- make_reg_aesi(outcome_variable = OV,
                                               start_year = 2002, end_year = 2020,
                                               crop_j = "fodder_crops",
                                               price_k = K_beef,
                                               SkPk= CTRL,
                                               extra_price_k = c())#"Sheep", "Pork", "Chicken"
    names(res_list_beef)[elm] <- paste0(OV,"_",CTRL, "_aesi")
    elm <- elm + 1
  }
}

# Prepare matrix for table
rm(ape_mat)
ape_mat <- bind_cols(lapply(res_list_beef, FUN = make_table_mat)) %>% as.matrix()
row.names(ape_mat) <- make_table_mat(res_list_beef[[1]]) %>% row.names()
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = "Indirect effects of global commodity markets on deforestation for cattle, 2002-2020") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Estimates" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1, "First loss" = 4,"Primary forest loss" = 4),
                   align = "c",
                   strikeout = F) %>%
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "7em",
              latex_valign = "b") 



### PALM OIL
K_oilpalm <- c("Soybean_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil",
               "Sugar", "Maize") # in prices spelling "Olive_oil", "Soybeans", "Soybean_meal",


res_list_oilpalm <- list()
elm <- 1

# Forest definition
outcome_variableS <- c("firstloss_glassgfc", "phtf_loss")

# with controls or not 
control_ornot <- c(FALSE, TRUE)

for(OV in outcome_variableS){
  for(CTRL in control_ornot){
    res_list_oilpalm[[elm]] <- make_reg_acay(outcome_variable = OV,
                                               start_year = 2002, end_year = 2020,
                                               crop_j = "Oilpalm",
                                               price_k = K_oilpalm,
                                               SkPk= CTRL,
                                               extra_price_k = c())#
    names(res_list_oilpalm)[elm] <- paste0(OV,"_",CTRL, "_acay")
    elm <- elm + 1
    
    res_list_oilpalm[[elm]] <- make_reg_aesi(outcome_variable = OV,
                                               start_year = 2002, end_year = 2020,
                                               crop_j = "Oilpalm",
                                               price_k = K_oilpalm,
                                               SkPk= CTRL,
                                               extra_price_k = c())#
    names(res_list_oilpalm)[elm] <- paste0(OV,"_",CTRL, "_aesi")
    elm <- elm + 1
  }
}

# Prepare matrix for table
rm(ape_mat)
ape_mat <- bind_cols(lapply(res_list_oilpalm, FUN = make_table_mat)) %>% as.matrix()
row.names(ape_mat) <- make_table_mat(res_list_oilpalm[[1]]) %>% row.names()
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = "Indirect effects of global commodity markets on deforestation for oil palm, 2002-2020") %>% #
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Estimates" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1, "First loss" = 4,"Primary forest loss" = 4),
                   align = "c",
                   strikeout = F) %>%
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "7em",
              latex_valign = "b") 



### SOY 
K_soy <- c("Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", 
           "Sugar", "Maize") # in prices spelling "Olive_oil",


res_list_soy <- list()
elm <- 1

# Forest definition
outcome_variableS <- c("firstloss_glassgfc", "phtf_loss")

# with controls or not 
control_ornot <- c(FALSE, TRUE)

for(OV in outcome_variableS){
  for(CTRL in control_ornot){
    res_list_soy[[elm]] <- make_reg_acay(outcome_variable = OV,
                                             start_year = 2002, end_year = 2020,
                                             crop_j = "Soybean",
                                             price_k = K_soy,
                                             SkPk= CTRL,
                                             extra_price_k = c())#
    names(res_list_soy)[elm] <- paste0(OV,"_",CTRL, "_acay")
    elm <- elm + 1
    
    res_list_soy[[elm]] <- make_reg_aesi(outcome_variable = OV,
                                             start_year = 2002, end_year = 2020,
                                             crop_j = "Soybean",
                                             price_k = K_soy,
                                             SkPk= CTRL,
                                             extra_price_k = c())#
    names(res_list_soy)[elm] <- paste0(OV,"_",CTRL, "_aesi")
    elm <- elm + 1
  }
}

# Prepare matrix for table
rm(ape_mat)
ape_mat <- bind_cols(lapply(res_list_soy, FUN = make_table_mat)) %>% as.matrix()
row.names(ape_mat) <- make_table_mat(res_list_soy[[1]]) %>% row.names()
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = "Indirect effects of global commodity markets on deforestation for soy, 2002-2020") %>% #
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Estimates" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1, "First loss" = 4,"Primary forest loss" = 4),
                   align = "c",
                   strikeout = F) %>%
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "7em",
              latex_valign = "b") 

### COCOA 
K_cocoa <- c("Coffee", "Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", "Sugar") # in prices spelling


res_list_cocoa <- list()
elm <- 1

# Forest definition
outcome_variableS <- c("firstloss_glassgfc", "phtf_loss")

# with controls or not 
control_ornot <- c(FALSE, TRUE)

for(OV in outcome_variableS){
  for(CTRL in control_ornot){
    res_list_cocoa[[elm]] <- make_reg_acay(outcome_variable = OV,
                                             start_year = 2002, end_year = 2020,
                                             crop_j = "Cocoa",
                                             price_k = K_cocoa,
                                             SkPk= CTRL,
                                             extra_price_k = c())#
    names(res_list_cocoa)[elm] <- paste0(OV,"_",CTRL, "_acay")
    elm <- elm + 1
    
    res_list_cocoa[[elm]] <- make_reg_aesi(outcome_variable = OV,
                                             start_year = 2002, end_year = 2020,
                                             crop_j = "Cocoa",
                                             price_k = K_cocoa,
                                             SkPk= CTRL,
                                             extra_price_k = c())#
    names(res_list_cocoa)[elm] <- paste0(OV,"_",CTRL, "_aesi")
    elm <- elm + 1
  }
}

# Prepare matrix for table
rm(ape_mat)
ape_mat <- bind_cols(lapply(res_list_cocoa, FUN = make_table_mat)) %>% as.matrix()
row.names(ape_mat) <- make_table_mat(res_list_cocoa[[1]]) %>% row.names()
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = "Indirect effects of global commodity markets on deforestation for cocoa, 2002-2020") %>% #
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Estimates" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1, "First loss" = 4,"Primary forest loss" = 4),
                   align = "c",
                   strikeout = F) %>%
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "7em",
              latex_valign = "b") 

### COFFEE 
K_coffee <- c("Tea", "Cocoa", "Sugar", "Tobacco") # in prices spelling


res_list_coffee <- list()
elm <- 1

# Forest definition
outcome_variableS <- c("firstloss_glassgfc", "phtf_loss")

# with controls or not 
control_ornot <- c(FALSE, TRUE)

for(OV in outcome_variableS){
  for(CTRL in control_ornot){
    res_list_coffee[[elm]] <- make_reg_acay(outcome_variable = OV,
                                             start_year = 2002, end_year = 2020,
                                             crop_j = "Coffee",
                                             price_k = K_coffee,
                                             SkPk= CTRL,
                                             extra_price_k = c())#
    names(res_list_coffee)[elm] <- paste0(OV,"_",CTRL, "_acay")
    elm <- elm + 1
    
    res_list_coffee[[elm]] <- make_reg_aesi(outcome_variable = OV,
                                             start_year = 2002, end_year = 2020,
                                             crop_j = "Coffee",
                                             price_k = K_coffee,
                                             SkPk= CTRL,
                                             extra_price_k = c())#
    names(res_list_coffee)[elm] <- paste0(OV,"_",CTRL, "_aesi")
    elm <- elm + 1
  }
}

# Prepare matrix for table
rm(ape_mat)
ape_mat <- bind_cols(lapply(res_list_coffee, FUN = make_table_mat)) %>% as.matrix()
row.names(ape_mat) <- make_table_mat(res_list_coffee[[1]]) %>% row.names()
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = "Indirect effects of global commodity markets on deforestation for coffee, 2002-2020") %>% #
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Estimates" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1, "First loss" = 4,"Primary forest loss" = 4),
                   align = "c",
                   strikeout = F) %>%
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "7em",
              latex_valign = "b") 




#### RUN AND PLOT ACAY ####
OV <- "phtf_loss"
SY <- 1983
EY <- 2020
if(OV == "first_loss"){main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_acay_long_final.Rdata"))}
if(OV == "firstloss_glassgfc"){main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_acay_long_final.Rdata"))}
if(OV == "phtf_loss"){main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_acay_long_final.Rdata"))}

### BEEF
K_beef <- c("Soybeans", "Soybean_meal", "Palm_oil", 
             "cereal_crops") # in prices spelling "Barley", "Oat", "Sorghum", "Maize", "Rice"

beef_res <- make_reg_acay(outcome_variable = OV,
                             start_year = SY, end_year = EY,
                             crop_j = "fodder_crops",
                             price_k = K_beef,
                          SkPk= F,
                             extra_price_k = c())#"Sheep", "Pork", "Chicken"

### PALM OIL
K_palmoil <- c("Soybean_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil",
               "Sugar", "Maize") # in prices spelling "Olive_oil", "Soybeans", "Soybean_meal",

oilpalm_res <- make_reg_acay(outcome_variable = OV,
                             start_year = SY, end_year = EY,
                             crop_j = "Oilpalm",
                             price_k = K_palmoil,
                             extra_price_k = "Crude_oil")


### SOY 
K_soy <- c("Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", 
           "Sugar", "Maize") # in prices spelling "Olive_oil",


soy_res <- make_reg_acay(outcome_variable = OV,
                         start_year = SY, end_year = EY, 
                         crop_j = "Soybean", 
                         price_k = K_soy, 
                         extra_price_k = "Crude_oil")


### COCOA 
K_cocoa <- c("Coffee", "Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", "Sugar") # in prices spelling


cocoa_res <- make_reg_acay(outcome_variable = OV,
                           start_year = SY, end_year = EY, 
                           crop_j = "Cocoa", 
                           price_k = K_cocoa)


### COFFEE 
K_coffee <- c("Tea", "Cocoa", "Sugar", "Tobacco") # in prices spelling


coffee_res <- make_reg_acay(outcome_variable = OV,
                            start_year = SY, end_year = EY, 
                            crop_j = "Coffee", 
                            price_k = K_coffee)


### PLOT COEFFICIENTS ### 
# prepare regression outputs in a tidy data frame readable by dwplot
rm(df)
df <- rbind(beef_res, oilpalm_res, soy_res, cocoa_res, coffee_res)#
df$model <- gsub(pattern = "_X_.*$", x = row.names(df), replacement = "") # replace everything after the first underscore with nothing
df$term <- sub(pattern = ".+?(_X_)", x = row.names(df), replacement = "") # replace everyting before the first underscore with nothing

names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"

if(DS == "phtfl"){
  title <- paste0("Indirect effects of commodity prices on the main agricultural drivers of global primary deforestation from ",SY," to ",EY)
}else{
  title <- paste0("Indirect effects of commodity prices on the main agricultural drivers of global deforestation from ",SY," to ",EY)
}
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
                         Soybeans = "Soybeans", 
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



#### RUN AND PLOT AESI ####
DS <- "glass"
SY <- 1983
EY <- 2020
if(DS == "glass"){main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aesi_long_final.Rdata"))}
if(DS == "fl8320"){main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_aesi_long_final.Rdata"))}
if(DS == "phtfl"){main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aesi_long_final.Rdata"))}


### PALM OIL
K_palmoil <- c("Soybean_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil",
               "Sugar", "Maize") # in prices spelling "Olive_oil", "Soybeans", "Soybean_meal",

oilpalm_res <- make_reg_aesi(dataset = DS,
                              start_year = SY, end_year = EY,
                              crop_j = "Oilpalm",
                              price_k = K_palmoil,
                              extra_price_k = "Crude_oil")


### SOY 
K_soy <- c("Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", 
           "Sugar", "Maize") # in prices spelling "Olive_oil",


soy_res <- make_reg_aesi(dataset = DS,
                        start_year = SY, end_year = EY, 
                        crop_j = "Soybean", 
                        price_k = K_soy, 
                        extra_price_k = "Crude_oil")


### COCOA 
K_cocoa <- c("Coffee", "Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", "Sugar") # in prices spelling


cocoa_res <- make_reg_aesi(dataset = DS,
                          start_year = SY, end_year = EY, 
                          crop_j = "Cocoa", 
                          price_k = K_cocoa)


### COFFEE 
K_coffee <- c("Tea", "Cocoa", "Sugar", "Tobacco") # in prices spelling


coffee_res <- make_reg_aesi(dataset = DS,
                            start_year = SY, end_year = EY, 
                            crop_j = "Coffee", 
                            price_k = K_coffee)


### PLOT COEFFICIENTS ### 
# prepare regression outputs in a tidy data frame readable by dwplot
rm(df)
df <- rbind(oilpalm_res, soy_res, cocoa_res, coffee_res)#
df$model <- gsub(pattern = "_.*$", x = row.names(df), replacement = "") # replace everything after the first underscore with nothing
df$term <- sub(pattern = ".+?(_)", x = row.names(df), replacement = "") # replace everyting before the first underscore with nothing

names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"

if(DS == "phtfl"){
  title <- paste0("Indirect effects of commodity prices on the main agricultural drivers of global primary deforestation from ",SY," to ",EY)
}else{
  title <- paste0("Indirect effects of commodity prices on the main agricultural drivers of global deforestation from ",SY," to ",EY)
}
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
                         Soybeans = "Soybeans", 
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



#### First loss (1983-2015) ####
main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aesi_long_final.Rdata"))

### PALM OIL
K_palmoil <- c("Soybean_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil",
               "Sugar", "Maize") # in prices spelling "Olive_oil", "Soybeans", "Soybean_meal",

oilpalm_res <- make_reg_aesi(#start_year = 1983, end_year = 2015,
  crop_j = "Oilpalm",
  price_k = K_palmoil,
  extra_price_k = "Crude_oil")


### SOY 
K_soy <- c("Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", 
           "Sugar", "Maize") # in prices spelling "Olive_oil",


soy_res <- make_reg_aesi(#start_year = 1983, end_year = 2015, 
  crop_j = "Soybean", 
  price_k = K_soy, 
  extra_price_k = "Crude_oil")


### COCOA 
K_cocoa <- c("Coffee", "Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", "Sugar") # in prices spelling


cocoa_res <- make_reg_aesi(#start_year = 1983, end_year = 2015, 
  crop_j = "Cocoa", 
  price_k = K_cocoa)


### COFFEE 
K_coffee <- c("Tea", "Cocoa", "Sugar", "Tobacco") # in prices spelling


coffee_res <- make_reg_aesi( #start_year = 1983, end_year = 2015, 
  crop_j = "Coffee", 
  price_k = K_coffee)


### PLOT COEFFICIENTS ### 
# prepare regression outputs in a tidy data frame readable by dwplot
df <- rbind(oilpalm_res, soy_res, cocoa_res, coffee_res)#
df$model <- gsub(pattern = "_.*$", x = row.names(df), replacement = "") # replace everything after the first underscore with nothing
df$term <- sub(pattern = ".+?(_)", x = row.names(df), replacement = "") # replace everyting before the first underscore with nothing

names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"

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
                         Soybeans = "Soybeans", 
                         Soybean_meal = "Soybean meal", 
                         Cocoa = "Cocoa", 
                         Coffee = "Coffee", 
                         Tea = "Tea", 
                         Tobacco = "Tobacco"
    )) +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    ggtitle("Indirect effects of commodity prices on the main agricultural drivers of deforestation") +
    theme(plot.title = element_text(face="bold", size=c(10)),
          legend.position = c(0.007, 0.001),
          legend.justification = c(0, 0), 
          legend.background = element_rect(colour="grey80"),
          legend.title = element_blank())}   #element_text("Direct drivers of tropical deforestation")


# table_title <- paste0("Indirect effects of commodity k prices on tropical deforestation for oil palm plantations, 1983-2015") 
# etable(oilpalm_res_list, 
#        digits = 1, 
#        tex = TRUE, 
#        title = table_title,
#        depvar = TRUE,
#        subtitles = paste0("k = ", K_commodities),
#        #drop = c("SkPk", "SjPj"),
#        #coefstat = "confint",
#        sdBelow = TRUE,
#        yesNo = "X",
#        fitstat = c("sq.cor"),
#        dict = TRUE, 
#        powerBelow = -7)




#### Primary forest loss (2002-2020) ####
main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aesi_long_final.Rdata"))

### PALM OIL
K_palmoil <- c("Soybean_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil",
               "Sugar", "Maize") # in prices spelling "Olive_oil", "Soybeans", "Soybean_meal",

oilpalm_res <- make_reg_aesi(dataset = "phtfl",
                             start_year = 2002, end_year = 2020,
                             crop_j = "Oilpalm",
                             price_k = K_palmoil,
                             extra_price_k = "Crude_oil")


### SOY 
K_soy <- c("Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", 
           "Sugar", "Maize") # in prices spelling "Olive_oil",


soy_res <- make_reg_aesi(dataset = "phtfl",
                         start_year = 2002, end_year = 2020, 
                         crop_j = "Soybean", 
                         price_k = K_soy, 
                         extra_price_k = "Crude_oil")


### COCOA 
K_cocoa <- c("Coffee", "Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", "Sugar") # in prices spelling


cocoa_res <- make_reg_aesi(dataset = "phtfl",
                           start_year = 2002, end_year = 2020, 
                           crop_j = "Cocoa", 
                           price_k = K_cocoa)


### COFFEE 
K_coffee <- c("Tea", "Cocoa", "Sugar", "Tobacco") # in prices spelling


coffee_res <- make_reg_aesi(dataset = "phtfl",
                            start_year = 2002, end_year = 2020, 
                            crop_j = "Coffee", 
                            price_k = K_coffee)


### PLOT COEFFICIENTS ### 
# prepare regression outputs in a tidy data frame readable by dwplot
df <- rbind(oilpalm_res, soy_res, cocoa_res, coffee_res)#
df$model <- gsub(pattern = "_.*$", x = row.names(df), replacement = "") # replace everything after the first underscore with nothing
df$term <- sub(pattern = ".+?(_)", x = row.names(df), replacement = "") # replace everyting before the first underscore with nothing

names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"

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
                         Soybeans = "Soybeans", 
                         Soybean_meal = "Soybean meal", 
                         Cocoa = "Cocoa", 
                         Coffee = "Coffee", 
                         Tea = "Tea", 
                         Tobacco = "Tobacco"
    )) +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    ggtitle("Indirect effects of commodity prices on the main agricultural drivers of deforestation") +
    theme(plot.title = element_text(face="bold", size=c(10)),
          legend.position = c(0.007, 0.001),
          legend.justification = c(0, 0), 
          legend.background = element_rect(colour="grey80"),
          legend.title = element_blank())}  



