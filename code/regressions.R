
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

rm(outcome_variable, start_year, end_year, crop_j, j_soy, price_k, extra_price_k, standardized_si, price_lag, SjPj, SkPk, fe, distribution, output, se, cluster, 
   controls, regressors, outcome_variable)

outcome_variable = "first_loss" # "first_loss", "firstloss_glassgfc", "phtf_loss"
start_year = 1983
end_year = 2020
price_info = "lag1"
further_lu_evidence = "none"
crop_j = "fodder_crops"
j_soy = "Soybean_oil"
# for the k variables hypothesized in overleaf for palm oil, feglm quasipoisson converges within 25 iter. 
# but maybe not with skPk controls. 
price_k <- K_beef
# c( "Soybean_oil", "Rapeseed_oil", "Sunflower_oil", 
#               "Coconut_oil", "Soybeans","Sugar", "Maize") 
# "Barley",  "Chicken", "Sheep", "Banana", "Beef", "Olive_oil",
#               "Orange",  "Cotton",  "Groundnuts",  "Rubber", "Sorghum","Cocoa",  "Coffee",
#                "Rice",   "Wheat",  "Palm_oil", ), 
#                "Tea", "Tobacco",  "Oat",  
#               , "Pork")
extra_price_k = c("Sheep", "Pork", "Chicken") # , 
standardized_si = TRUE
price_lag = 1
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
# rm(dataset, start_year, end_year, crop_j, j_soy, price_k, standardized_si, price_lag, SjPj, SkPk, fe, distribution, output, se, cluster,     controls, regressors, outcome_variable)

make_reg_acay <- function(outcome_variable = "first_loss", # one of "first_loss", "firstloss_glassgfc", "phtf_loss"
                          price_info = "lag1", # one of "lag1", "2pya", "3pya", "4pya", "5pya",
                          start_year = 2002, 
                          end_year = 2015, 
                          further_lu_evidence = "none", # either "none", "sbqt_direct_lu", or "sbqt"mode_lu"
                          crop_j = "Oilpalm", # in GAEZ spelling
                          j_soy = "Soybean_oil", # in case crop_j is Soybean, which price should be used? One of "Soybeans", "Soybean_oil", "Soybean_meal".
                          price_k = c("Sugar", "Maize"), # in prices spelling
                          extra_price_k = c(), # One of "Crude_oil", "Chicken", "Pork", "Sheep" 
                          SjPj = TRUE,
                          SkPk = FALSE,
                          #commo_m = c(""), comment coder ça pour compatibilité avec loops over K_commo ? 
                          fe = "grid_id + country_year", 
                          distribution = "gaussian",#  "quasipoisson", 
                          se = "cluster", 
                          cluster ="grid_id",
                          coefstat = "confint", # one of "se", "tstat", "confint"
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
  if(outcome_variable=="first_loss"){
    vars <- c("grid_id", "year", "lon", "lat", "country_name", "country_year", outcome_variable, "sbqt_direct_lu", "sbqt_mode_lu",
              names(d)[grepl(pattern ="_std", x= names(d))])
  }else {
    vars <- c("grid_id", "year", "lon", "lat", "country_name", "country_year", outcome_variable, 
              names(d)[grepl(pattern ="_std", x= names(d))])
  } 
  d <- dplyr::select(d, all_of(vars))


  ### MAKE FINAL REGRESSION VARIABLES 
  
  # keep only the cells with positive rj (since we need to divide by rj)
  d <- dplyr::filter(d, !!as.symbol(revenue_j_std) > 0)
  
  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", price_all)], by = "year")
  
  # Main regressors
  regressors <- c()
  for(Pk in price_all){
    varname <- paste0(crop_j, "_r_", original_price_all[match(Pk, price_all)])
    regressors <- c(regressors, varname)
    d <- mutate(d, 
                !!as.symbol(varname) := !!as.symbol(Pk)/(!!as.symbol(revenue_j_std)))
  }
  rm(varname)
  
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
  
  # this is necessary to compute SE as we want to.  
  reg_res <- summary(reg_res, se = se,
                     cluster = cluster)
  
  # Now keep only information necessary, otherwise the output of fixest estimation is large and we can't collect too many at the same time (over loops)  
  df_res <- reg_res$coeftable
  
  # Keep only variable of interest, and rename it
  df_res <- df_res[!grepl(pattern = "ctrl_", x = row.names(df_res)),]
  
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
                          end_year = 2015, 
                          further_lu_evidence = "none", # either "none", "sbqt_direct_lu", or "sbqt"mode_lu"
                          crop_j = "Oilpalm", # in GAEZ spelling
                          j_soy = "Soybean_oil", # in case crop_j is Soybean, which price should be used? One of "Soybeans", "Soybean_oil", "Soybean_meal".
                          price_k = c("Sugar", "Maize"), # in prices spelling
                          extra_price_k = c(), # One of "Crude_oil", "Chicken", "Pork", "Sheep" 
                          standardized_si = TRUE,
                          price_lag = 1, 
                          SjPj = TRUE,
                          SkPk = FALSE,
                          #commo_m = c(""), comment coder ça pour compatibilité avec loops over K_commo ? 
                          fe = "grid_id + country_year", 
                          distribution = "gaussian",#  "quasipoisson", 
                          se = "cluster", 
                          cluster ="grid_id",
                          coefstat = "confint", # one of "se", "tstat", "confint"
                          output = "coef_table" # one of "data", "est_obj", "coef_table" 
                          ){

   
  ### DATA 
  d <- main_data

  
  # ## Outcome variable
  # if(dataset=="glass"){
  #   outcome_variable <- "first_loss"} #   "sbqt_direct_lu"      "sbqt_mode_lu" 
  # if(dataset=="fl8320"){
  #   outcome_variable <- "firstloss_glassgfc"}
  # if(dataset=="phtfl"){
  #   outcome_variable <- "phtf_loss"}
  
  # restrict the outcome variable to evidence of further LU
  if(outcome_variable == "first_loss" & further_lu_evidence != "none"){
    if(crop_j == "fodder_crops"){
      d <- dplyr::mutate(d, lu_evidence = (!!as.symbol(further_lu_evidence) == 30)) # either direct or mode subsequent lu are grassland
    }
    if(crop_j == "Soybean"){
      d <- dplyr::mutate(d, lu_evidence = (!!as.symbol(further_lu_evidence) == 10)) # either direct or mode subsequent lu are cropland
    }
    if(crop_j %in% c("Oilpalm", "Cocoa", "Coffee")){
      d <- dplyr::mutate(d, lu_evidence = (sbqt_direct_lu != 20 & sbqt_mode_lu == 20)) # is not forest in the year directly after, but is forest again in the mode subsequent lu. 
    }
  }
  
#### SPECIFICATIONS    
  
  # We need to create the variables we will need, based on:
  # crop_j, price_k, whether the price is lagged, whether the suitability is standardized, and the controls we want 
  
  
  ### PREPARE EQUIVALENCES BETWEEN SUITABILITY AND PRICE
  
  ## Identify the price of j based on crop_j (GAEZ spelling)
  if(crop_j == "Pasture"){price_j <- "Beef"}
  if(crop_j == "Oilpalm"){price_j <- "Palm_oil"}
  if(crop_j == "Soybean"){price_j <- j_soy}
  if(crop_j == "Cocoa"){price_j <- "Cocoa"}
  if(crop_j == "Coffee"){price_j <- "Coffee"}
  
  ## Do the revert for k: from price of k to the GAEZ crop 
  
  # code like that allows to order corresponding crops in crop_k in the same order as in price_k, irrespectively of the order prices are passed to price_k 
  crop_k <- c()
  if("Banana"%in% price_k){crop_k[match("Banana", price_k)] <- "Banana"}
  if("Barley"%in% price_k){crop_k[match("Barley", price_k)] <- "Barley"}
  # if("Beef"%in% price_k){crop_k[match("Beef", price_k)] <- "Pasture"} # currently does not exist
  if("Orange"%in% price_k){crop_k[match("Orange", price_k)] <- "Citrus"}
  if("Cocoa"%in% price_k){crop_k[match("Cocoa", price_k)] <- "Cocoa"}
  if("Coconut_oil"%in% price_k){crop_k[match("Coconut_oil", price_k)] <- "Coconut"}
  if("Coffee"%in% price_k){crop_k[match("Coffee", price_k)] <- "Coffee"}
  if("Cotton"%in% price_k){crop_k[match("Cotton", price_k)] <- "fibre_crops"}
  if("Groundnuts"%in% price_k){crop_k[match("Groundnuts", price_k)] <- "Groundnut"}
  if("Maize"%in% price_k){crop_k[match("Maize", price_k)] <- "Maize"}
  if("Oat"%in% price_k){crop_k[match("Oat", price_k)] <- "Oat"} # aggregate different crops there ? like rye...
  if("Olive_oil"%in% price_k){crop_k[match("Olive_oil", price_k)] <- "Olive"}
  if("Palm_oil"%in% price_k){crop_k[match("Palm_oil", price_k)] <- "Oilpalm"}
  if("Rapeseed_oil"%in% price_k){crop_k[match("Rapeseed_oil", price_k)] <- "Rapeseed"}
  if("Rice"%in% price_k){crop_k[match("Rice", price_k)] <- "rice_crops"}
  if("Sorghum"%in% price_k){crop_k[match("Sorghum", price_k)] <- "Sorghum"}
  if("Soybeans"%in% price_k){crop_k[match("Soybeans", price_k)] <- "Soybean"}
  if("Soybean_oil"%in% price_k){crop_k[match("Soybean_oil", price_k)] <- "Soybean"}
  if("Soybean_meal"%in% price_k){crop_k[match("Soybean_meal", price_k)] <- "Soybean"}
  if("Sugar"%in% price_k){crop_k[match("Sugar", price_k)] <- "sugar_crops"}
  if("Sunflower_oil"%in% price_k){crop_k[match("Sunflower_oil", price_k)] <- "Sunflower"}
  if("Tea"%in% price_k){crop_k[match("Tea", price_k)] <- "Tea"}
  if("Tobacco"%in% price_k){crop_k[match("Tobacco", price_k)] <- "Tobacco"}
  if("Wheat"%in% price_k){crop_k[match("Wheat", price_k)] <- "Wheat"}
  if("cereal_crops"%in% price_k){crop_k[match("cereal_crops", price_k)] <- "cereal_crops"}
  if("oil_crops"%in% price_k){crop_k[match("oil_crops", price_k)] <- "oil_crops"}

  
  # handle commodities that do not have SI
  # if("Crude_oil" %in% price_k){crop_k <- c(crop_k, NULL)}
  # if("Rubber" %in% price_k){crop_k <- c(crop_k, NULL)}
  # if("Chicken" %in% price_k){crop_k <- c(crop_k, NULL)}
  # if("Pork" %in% price_k){crop_k <- c(crop_k, NULL)}
  # if("Sheep" %in% price_k){crop_k <- c(crop_k, NULL)}
  
  # this is mostly useless
  if(length(crop_k) == 0){SkPk <- FALSE}

  
  ### SELECT PRICE AND SUITABILITY VARIABLES 
  # first, save price_k for later purpose 
  original_price_j <- price_j 
  original_price_k <- price_k 
  original_extra_price_k <- extra_price_k 
  original_price_all <- c(price_j, price_k, extra_price_k)
  
  # For j
  if(price_lag != 0){price_j <- paste0(price_j,"_lag",price_lag)}
  price_j <- paste0("ln_", price_j)
  
  # For k 
  if(price_lag != 0){price_k <- paste0(price_k,"_lag",price_lag)}
  # don't condition that, as we won't use non logged prices a priori. 
  price_k <- paste0("ln_", price_k)

  # For extra commodities k 
  if(length(extra_price_k) > 0){
    if(price_lag != 0){extra_price_k <- paste0(extra_price_k,"_lag",price_lag)}
    extra_price_k <- paste0("ln_", extra_price_k)
  }
  
  # group all prices transformed
  price_all <- c(price_j, price_k, extra_price_k)
  
  # Identify the suitability index needed - don't condition that either for now, as the non stded indexes are not in the main dataset now (for memory issues) we'll see later how to conduct the rob check
  suitability_j <- paste0(crop_j, "_std")
  suitability_k <- paste0(crop_k, "_std")
  
  # Keep only in data the useful variables 
  d <- dplyr::select(d, all_of(c("grid_id", "year", "country_name", "country_year",
                               outcome_variable, suitability_j, suitability_k))) 
    
  # and keep only the cells with positive suitability for crop j 
  d <- dplyr::filter(d, !!as.symbol(suitability_j) > 0)

  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", price_all)], by = "year")
  
  
  ### MAKE NECESSARY VARIABLES 

  # Main regressors
  regressors <- c()
  for(k in price_all){
    varname <- paste0(crop_j, "_", original_price_all[match(k, price_all)])
    regressors <- c(regressors, varname)
    d <- mutate(d, 
                !!as.symbol(varname) := !!as.symbol(k)/(!!as.symbol(suitability_j)))
  }
  
  # Controls
  controls <- c()
  if(SkPk){
    for(i in 1:length(suitability_k)){
      varname <- paste0("ctrl_", original_price_k[i])
      controls <- c(controls, varname)
      d <- mutate(d, 
                  !!as.symbol(varname) := !!as.symbol(suitability_k[i])*!!as.symbol(price_k[i]))
    }
  }
  

  # # if(SkPk){controls <- c(controls, "SkPk")}
  # if(SjPj){controls <- c(controls, "SjPj")} 
  
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
  
  # this is necessary to compute SE as we want to.  
  reg_res <- summary(reg_res, se = se,
                    cluster = cluster)
  
  # Now keep only information necessary, otherwise the output of fixest estimation is large and we can't collect too many at the same time (over loops)  
  df_res <- reg_res$coeftable
  
  # Keep only variable of interest, and rename it
  df_res <- df_res[!grepl(pattern = "ctrl_", x = row.names(df_res)),]
  
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



