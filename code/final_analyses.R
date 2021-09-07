
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
                   "ggplot2", "dotwhisker", "tmap",# "leaflet", "htmltools"
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
                 first_loss = "First forest loss",
                 Soybean_meal = "Soybean meal",
                 Soybean_oil = "Soybean oil",
                 Olive_oil = "Olive oil",
                 Rapeseed_oil = "Rapeseed oil",
                 Sunflower_oil = "Sunflower oil",
                 Coconut_oil = "Coconut oil"))


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


### READ ALL POSSIBLE DATASETS HERE 
# but not all together because of memory issues. 
# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aesi_long_final.Rdata"))

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_aesi_long_final.Rdata"))

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aesi_long_final.Rdata"))

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aeay_long_final.Rdata"))

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_aeay_long_final.Rdata"))

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aeay_long_final.Rdata"))

prices <- readRDS(here("temp_data", "prepared_international_prices.Rdata"))



make_reg_aeay <- function(outcome_variable = "driven_loss", # one of "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
                          start_year = 2001, 
                          end_year = 2020, 
                          continent = "all", # one of "Africa", "America", "Asia", or "all"
                          further_lu_evidence = "none", # either "none", "sbqt_direct_lu", or "sbqt"mode_lu"
                          crop_j = "Oilpalm", # in GAEZ spelling
                          price_k = c("Sugar", "Maize"), # in prices spelling
                          extra_price_k = c(), # One of "Crude_oil", "Chicken", "Pork", "Sheep" 
                          price_info = "lag1", # one of "lag1", "2pya", "3pya", "4pya", "5pya",
                          fcr = 7.8, # feed conversion ratio according to Wilkinson, 2010. 
                          SjPj = TRUE,
                          SkPk = FALSE,
                          #commo_m = c(""), comment coder ça pour compatibilité avec loops over K_commo ? 
                          fe = "grid_id + country_name", 
                          distribution = "quasipoisson",#  "quasipoisson", 
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
  # Identify only Prices and Crops that can be matched 
  # See below, in conversion part, why we exclude coconut and cotton. 
  # To determine land use, only the potential revenue of soybeans is considered (for simplicity) 
  
  all_crop_prices <- mapmat[mapmat[,"Prices"]!="Soybean_oil" & 
                              mapmat[,"Prices"]!="Soybean_meal" & 
                              mapmat[,"Prices"]!="Coconut_oil" &
                              mapmat[,"Prices"]!="Cotton", "Prices"]
  all_crop_aeay <- mapmat[mapmat[,"Prices"]!="Soybean_oil" & 
                            mapmat[,"Prices"]!="Soybean_meal" & 
                            mapmat[,"Prices"]!="Coconut_oil" &
                            mapmat[,"Prices"]!="Cotton", "Crops"]
  
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
  if(outcome_variable == "first_loss" | outcome_variable == "nd_first_loss"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aeay_long_final.Rdata"))}
  if(outcome_variable == "firstloss_glassgfc"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_aeay_long_final.Rdata"))}
  if(outcome_variable == "phtf_loss"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aeay_long_final.Rdata"))}
  if(outcome_variable == "driven_loss"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "driverloss_aeay_long_final.Rdata"))}
  
  # # restrict the outcome variable to evidence of further LU
  # d <- dplyr::mutate(d, lu_evidence = TRUE)
  # if(outcome_variable == "first_loss" & further_lu_evidence != "none"){
  #    # , cropland, or tree plantation
  #   # if it is actually forest afterwards, this means that 
  #   
  #   if(crop_j == "Fodder"){
  #     d <- dplyr::mutate(d, lu_evidence = (!!as.symbol(further_lu_evidence) == 30)) # either direct or mode subsequent lu are grassland
  #   }
  #   if(crop_j == "Soybean"){
  #     d <- dplyr::mutate(d, lu_evidence = (!!as.symbol(further_lu_evidence) == 10)) # either direct or mode subsequent lu are cropland
  #   }
  #   if(crop_j %in% c("Oilpalm", "Cocoa", "Coffee")){
  #     d <- dplyr::mutate(d, lu_evidence = (sbqt_direct_lu != 20 & sbqt_mode_lu == 20)) # is not forest in the year directly after, but is forest again in the mode subsequent lu. 
  #   }
  # }
  
  ### Keep only variables needed
  # typically removes only aeay of crops that we do not use anyway because they match no price. 
  # and in glass, the sbqt_lu variables. 
  var2keep <- c("grid_id", "year", "lon", "lat", "country_name", "country_year", "continent_name", outcome_variable, all_crop_aeay)
  d <- dplyr::select(d, all_of(var2keep))
  rm(var2keep)
  ### PREPARE rj, the standardized achievable revenues
  
  # First, convert aeay quantities into TON/ha. Note that:
  # "For most crops the agro-climatic potential yield is given as kg dry weight per hectare. 
  # For alfalfa, miscanthus, switchgrass, reed canary grass, napier grass, pasture legumes and grasses the yield is given in 10 kg dry weight per hectare. 
  # For sugar beet and sugarcane (and hence Sugar, the max of them) yields are in kg sugar per hectare, 
  # and for oil palm and olives in kg oil per hectare. Cotton yield is given as kg lint per hectare." 
  # https://gaez.fao.org/pages/theme-details-theme-3
  names(d)
  # Convert those in kg/ha into ton/ha
  d <- dplyr::mutate(d, across(.cols = all_of(all_crop_aeay),
                               .fns = ~./1000)) 
  # Convert Nappier grass and alfalfa (components of Fodder) from now 10 tons to tons. 
  d <- dplyr::mutate(d, across(.cols = all_of("Fodder"),
                               .fns = ~.*10)) 
  
  # Those in different commodity format in GAEZ AEAY and in prices. 
  
  # We do not convert cotton from lint to kg, because we do not use cotton, as it is a crop appart.
  
  # For rapeseed and sunflower, the extraction rate seem very close to 1 
  # as per the crude fat to DM ratios at:
  # https://feedtables.com/content/rapeseed-oil
  # https://feedtables.com/content/sunflower-oil
  
  # Coconut cannot be converted to coconut oil, as it is made of only a by product of coconut; 
  # We would not correctly estimate the value of coconut if we counted the whole coconut yield as the copra byproduct;   
  
  # Olive and Oilpalm are expressed in oil already, as in the price data. Hence nothing to do. 
  
  # Convert fodder crop yield into beef by a feed conversion ratio. 
  # Following Alexander et al. 2016, 25 tons of feed transform into 1 ton of beef meat (in edible weight). 
  # Thus, every ton of feed counted in R_Fodder (coming from agro-climatically achievable yields in ton/ha) is scaled to 1/25 ton of beef meat  
  d <- dplyr::mutate(d, Fodder = Fodder/fcr)
  
  ## Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", prices_4_revenue)], by = c("year"))
  
  # Prices have been converted to $/t in prepare_prices.R
  for(i in 1:length(prices_4_revenue)){
    price_i <- prices_4_revenue[i]
    revenue_i <- paste0("R_", all_crop_aeay[i])
    aeay_i <- all_crop_aeay[i]
    d <- dplyr::mutate(d, 
                       !!as.symbol(revenue_i) := !!as.symbol(aeay_i)*!!as.symbol(price_i))
  }
  rm(price_i, revenue_i, aeay_i)                  
  
  
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
  vars <- c("grid_id", "year", "lon", "lat", "country_name", "country_year", "continent_name", outcome_variable, 
            names(d)[grepl(pattern ="_std", x= names(d))])
  d <- dplyr::select(d, all_of(vars))
  rm(vars)
  
  ### MAKE FINAL REGRESSION VARIABLES 
  
  # keep only the cells with positive rj (since we need to divide by rj)
  d <- dplyr::filter(d, !!as.symbol(revenue_j_std) > 0)
  
  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", price_all)], by = c("year"))
  
  # Main regressors
  regressors <- c()
  for(Pk in price_all){
    varname <- paste0(crop_j, "_X_", original_price_all[match(Pk, price_all)])
    regressors <- c(regressors, varname)
    d <- mutate(d, 
                !!as.symbol(varname) := !!as.symbol(Pk)/(!!as.symbol(revenue_j_std)))
  }
  rm(varname)
  
  # Remove the first of the regressors if we do not control for SjPj. The first is always SjPj because of the construction of price_all that has price_j first. 
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
                  !!as.symbol(varname) := !!as.symbol(Pk)/(!!as.symbol(rk)+0.00000001))
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
  
  # - are in study area
  if(continent != "all"){
    d <- dplyr::filter(d, continent_name == continent)
  }
  
  # # - have a lu_evidence 
  # d <- dplyr::filter(d, lu_evidence)
  
  used_vars <- c("grid_id", "year", "country_name", "country_year", "continent_name",
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
  obstormv <- obs2remove(fml = as.formula(paste0(outcome_variable, " ~ ", fe)),
                         d, 
                         family = "poisson")
  # this is necessary to handle cases when there is no obs. to remove
  if(length(obstormv)>0){
    d_clean <- d[-obstormv,]
  } else {
    d_clean <- d
  }
  
  rm(d)
  
  ### REGRESSIONS
  
  if(distribution != "negbin"){ # i.e. if it's poisson or quasipoisson or gaussian
    reg_res <- fixest::feglm(fe_model,
                             data = d_clean, 
                             family = distribution,#distribution,#   # "poisson" ,
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


make_reg_aesi <- function(outcome_variable = "driven_loss", # one of "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
                          start_year = 2001, 
                          end_year = 2020, 
                          continent = "all", # one of "Africa", "America", "Asia", or "all"
                          further_lu_evidence = "none", # either "none", "sbqt_direct_lu", or "sbqt"mode_lu"
                          crop_j = "Oilpalm", # in GAEZ spelling
                          price_k = c("Sugar", "Maize"), # in prices spelling
                          extra_price_k = c(), # One of "Crude_oil", "Chicken", "Pork", "Sheep" 
                          price_info = "lag1", # one of "lag1", "2pya", "3pya", "4pya", "5pya",
                          SjPj = TRUE,
                          SkPk = FALSE,
                          #commo_m = c(""), comment coder ça pour compatibilité avec loops over K_commo ? 
                          fe = "grid_id + country_name", 
                          distribution = "quasipoisson",#  "quasipoisson", 
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
  if(outcome_variable == "first_loss" | outcome_variable == "nd_first_loss"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aesi_long_final.Rdata"))}
  if(outcome_variable == "firstloss_glassgfc"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_aesi_long_final.Rdata"))}
  if(outcome_variable == "phtf_loss"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aesi_long_final.Rdata"))}
  if(outcome_variable == "driven_loss"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "driverloss_aesi_long_final.Rdata"))}
  
  # Keep only in data the useful variables 
  d <- dplyr::select(d, all_of(c("grid_id", "year", "country_name", "country_year", "continent_name", outcome_variable, 
                                 suitability_j_std, suitability_k_std))) 
  
  # and keep only the cells with positive suitability for crop j 
  # Important to do this here, in aesi case, as we divide by sj in main regressor construction
  d <- dplyr::filter(d, !!as.symbol(suitability_j_std) > 0)
  
  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", price_all)], by = c("year"))
  
  
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
                  !!as.symbol(varname) := !!as.symbol(Pk)/(!!as.symbol(sk)+0.00000001))
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
  
  # - are in study area
  if(continent != "all"){
    d <- dplyr::filter(d, continent_name == continent)
  }
  
  used_vars <- c("grid_id", "year", "country_name", "country_year",  "continent_name",
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
  obstormv <- obs2remove(fml = as.formula(paste0(outcome_variable, " ~ ", fe)),
                         d, 
                         family = "poisson")
  # this is necessary to handle cases when there is no obs. to remove
  if(length(obstormv)>0){
    d_clean <- d[-obstormv,]
  } else {
    d_clean <- d
  }
  
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
                           rounding=5
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
  
  ## Mark whether there were controls or not
  # if(SkPk){
  #   mat_res <- rbind(mat_res, "X")
  # } else {
  #   mat_res <- rbind(mat_res, "")
  # }
  # row.names(mat_res)[nrow(mat_res)] <- "Controls"
  
  # Indicate the number of observations and clusters
  mat_res <- rbind(mat_res, unique(df_res$Observations))
  row.names(mat_res)[nrow(mat_res)] <- "Observations"
  
  mat_res <- rbind(mat_res, unique(df_res$Clusters))
  row.names(mat_res)[nrow(mat_res)] <- "Clusters"
  return(mat_res)
  
}


#### CONSTRUCT SPECIFICATION COMPARATIVE TABLES ####

# We build 6 tables, one for each j crop. 
# These tables compare aeay (rj) and aesi (sj) models, with and without controlling for SkPk, with different price info  

# Final parameters 
SY <- 2001
EY <- 2019
PI <- "4pya"
CTRL <- TRUE
K_extra <- c("Chicken", "Pork")


### BEEF
K_beef <- c("Soybean", 
            "Maize", "Wheat", "Barley", "Oat") # in prices spelling "cereal_crops" "Barley", "Oat", "Sorghum", "Maize", "Rice"

res_list_beef <- list()
elm <- 1

res_list_beef[[elm]] <- make_reg_aeay(price_info = PI,
                                      start_year = SY, end_year = EY,
                                      crop_j = "Fodder",
                                      price_k = K_beef,
                                      SkPk= CTRL,
                                      extra_price_k = K_extra)#"Sheep", "Pork", "Chicken"
names(res_list_beef)[elm] <- paste0(PI,"_",CTRL, "_aeay")
elm <- elm + 1

res_list_beef[[elm]] <- make_reg_aesi(price_info = PI,
                                      start_year = SY, end_year = EY,
                                      crop_j = "Fodder",
                                      price_k = K_beef,
                                      SkPk= CTRL,
                                      extra_price_k = K_extra)#"Sheep", "Pork", "Chicken"
names(res_list_beef)[elm] <- paste0(PI,"_",CTRL, "_aesi")
elm <- elm + 1



# Prepare matrix for table
rm(ape_mat)
ape_mat <- bind_cols(lapply(res_list_beef, FUN = make_table_mat)) %>% as.matrix()
row.names(ape_mat) <- make_table_mat(res_list_beef[[1]]) %>% row.names()
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = paste0("Indirect effects of global commodity markets on deforestation for cattle, ",SY,"-",EY)) %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Commodities \ models" = 1,"PR" = 1,"SI" = 1),
                   bold = F,
                   align = "c") %>%
  # add_header_above(c(" " = 1, "lag1" = 4,"4pya" = 4),
  #                  align = "c",
  #                  strikeout = F) %>%
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "7em",
              latex_valign = "b") 



### PALM OIL
K_oilpalm <- c("Soybean", "Rapeseed_oil", "Sunflower_oil",
               "Sugar", "Maize") # in prices spelling. We cannot have "Coconut_oil", because GAEZ yield is only expressed as dry matter and conversion is not possible (see above)
# "Olive_oil", "Soybean", "Soybean_meal",


res_list_oilpalm <- list()
elm <- 1

res_list_oilpalm[[elm]] <- make_reg_aeay(price_info = PI,
                                         start_year = SY, end_year = EY,
                                         crop_j = "Oilpalm",
                                         price_k = K_oilpalm,
                                         SkPk= CTRL,
                                         extra_price_k = c())#
names(res_list_oilpalm)[elm] <- paste0(PI,"_",CTRL, "_aeay")
elm <- elm + 1

res_list_oilpalm[[elm]] <- make_reg_aesi(price_info = PI,
                                         start_year = SY, end_year = EY,
                                         crop_j = "Oilpalm",
                                         price_k = K_oilpalm,
                                         SkPk= CTRL,
                                         extra_price_k = c())#
names(res_list_oilpalm)[elm] <- paste0(PI,"_",CTRL, "_aesi")
elm <- elm + 1



# Prepare matrix for table
rm(ape_mat)
ape_mat <- bind_cols(lapply(res_list_oilpalm, FUN = make_table_mat)) %>% as.matrix()
row.names(ape_mat) <- make_table_mat(res_list_oilpalm[[1]]) %>% row.names()
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = paste0("Indirect effects of global commodity markets on deforestation for oil palm, ",SY,"-",EY)) %>% #
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Commodities \ models" = 1,"PR" = 1,"SI" = 1),
                   bold = F,
                   align = "c") %>%
  # add_header_above(c(" " = 1, "lag1" = 4,"4pya" = 4),
  #                  align = "c",
  #                  strikeout = F) %>%
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "7em",
              latex_valign = "b") 



### SOY 
K_soy <- c("Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Beef",
           "Sugar", "Maize") # in prices spelling "Olive_oil",


res_list_soy <- list()
elm <- 1

res_list_soy[[elm]] <- make_reg_aeay(price_info = PI,
                                     start_year = SY, end_year = EY,
                                     crop_j = "Soybean",
                                     price_k = K_soy,
                                     extra_price_k = extra_K_price,
                                     SkPk= CTRL)#
names(res_list_soy)[elm] <- paste0(PI,"_",CTRL, "_aeay")
elm <- elm + 1

res_list_soy[[elm]] <- make_reg_aesi(price_info = PI,
                                     start_year = SY, end_year = EY,
                                     crop_j = "Soybean",
                                     price_k = K_soy,
                                     extra_price_k = extra_K_price,
                                     SkPk= CTRL)#
names(res_list_soy)[elm] <- paste0(PI,"_",CTRL, "_aesi")
elm <- elm + 1



# Prepare matrix for table
rm(ape_mat)
ape_mat <- bind_cols(lapply(res_list_soy, FUN = make_table_mat)) %>% as.matrix()
row.names(ape_mat) <- make_table_mat(res_list_soy[[1]]) %>% row.names()
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = paste0("Indirect effects of global commodity markets on deforestation for soy, ",SY,"-",EY)) %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Commodities \ models" = 1,"PR" = 1,"SI" = 1),
                   bold = F,
                   align = "c") %>%
  # add_header_above(c(" " = 1, "lag1" = 4,"4pya" = 4),
  #                  align = "c",
  #                  strikeout = F) %>%
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "7em",
              latex_valign = "b") 

### COCOA 
K_cocoa <- c("Coffee", "Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Sugar") # in prices spelling


res_list_cocoa <- list()
elm <- 1

res_list_cocoa[[elm]] <- make_reg_aeay(price_info = PI,
                                       start_year = SY, end_year = EY,
                                       crop_j = "Cocoa",
                                       price_k = K_cocoa,
                                       SkPk= CTRL,
                                       extra_price_k = c())#
names(res_list_cocoa)[elm] <- paste0(PI,"_",CTRL, "_aeay")
elm <- elm + 1

res_list_cocoa[[elm]] <- make_reg_aesi(price_info = PI,
                                       start_year = SY, end_year = EY,
                                       crop_j = "Cocoa",
                                       price_k = K_cocoa,
                                       SkPk= CTRL,
                                       extra_price_k = c())#
names(res_list_cocoa)[elm] <- paste0(PI,"_",CTRL, "_aesi")
elm <- elm + 1



# Prepare matrix for table
rm(ape_mat)
ape_mat <- bind_cols(lapply(res_list_cocoa, FUN = make_table_mat)) %>% as.matrix()
row.names(ape_mat) <- make_table_mat(res_list_cocoa[[1]]) %>% row.names()
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = paste0("Indirect effects of global commodity markets on deforestation for cocoa, ",SY,"-",EY)) %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Commodities \ models" = 1,"PR" = 1,"SI" = 1),
                   bold = F,
                   align = "c") %>%
  # add_header_above(c(" " = 1, "lag1" = 4,"4pya" = 4),
  #                  align = "c",
  #                  strikeout = F) %>%
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

res_list_coffee[[elm]] <- make_reg_aeay(price_info = PI,
                                        start_year = SY, end_year = EY,
                                        crop_j = "Coffee",
                                        price_k = K_coffee,
                                        SkPk= CTRL,
                                        extra_price_k = c())#
names(res_list_coffee)[elm] <- paste0(PI,"_",CTRL, "_aeay")
elm <- elm + 1

res_list_coffee[[elm]] <- make_reg_aesi(price_info = PI,
                                        start_year = SY, end_year = EY,
                                        crop_j = "Coffee",
                                        price_k = K_coffee,
                                        SkPk= CTRL,
                                        extra_price_k = c())#
names(res_list_coffee)[elm] <- paste0(PI,"_",CTRL, "_aesi")
elm <- elm + 1



# Prepare matrix for table
rm(ape_mat)
ape_mat <- bind_cols(lapply(res_list_coffee, FUN = make_table_mat)) %>% as.matrix()
row.names(ape_mat) <- make_table_mat(res_list_coffee[[1]]) %>% row.names()
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = paste0("Indirect effects of global commodity markets on deforestation for coffee, ",SY,"-",EY)) %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Commodities \ models" = 1,"PR" = 1,"SI" = 1),
                   bold = F,
                   align = "c") %>%
  # add_header_above(c(" " = 1, "lag1" = 4,"4pya" = 4),
  #                  align = "c",
  #                  strikeout = F) %>%
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "7em",
              latex_valign = "b") 


### RUBBER 
K_rubber <- c("Palm_oil", "Soybean", "Rapeseed_oil", "Sunflower_oil", 
              "Sugar", "Maize") # in prices spelling

res_list_rubber <- list()
elm <- 1

res_list_rubber[[elm]] <- make_reg_aeay(price_info = PI,
                                        start_year = SY, end_year = EY,
                                        crop_j = "Rubber",
                                        price_k = K_rubber,
                                        SkPk= CTRL,
                                        extra_price_k = c())#
names(res_list_rubber)[elm] <- paste0(PI,"_",CTRL, "_aeay")
elm <- elm + 1

res_list_rubber[[elm]] <- make_reg_aesi(price_info = PI,
                                        start_year = SY, end_year = EY,
                                        crop_j = "Rubber",
                                        price_k = K_rubber,
                                        SkPk= CTRL,
                                        extra_price_k = c())#
names(res_list_rubber)[elm] <- paste0(PI,"_",CTRL, "_aesi")
elm <- elm + 1



# Prepare matrix for table
rm(ape_mat)
ape_mat <- bind_cols(lapply(res_list_rubber, FUN = make_table_mat)) %>% as.matrix()
row.names(ape_mat) <- make_table_mat(res_list_rubber[[1]]) %>% row.names()
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = paste0("Indirect effects of global commodity markets on deforestation for rubber, ",SY,"-",EY)) %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Commodities \ models" = 1,"PR" = 1,"SI" = 1),
                   bold = F,
                   align = "c") %>%
  # add_header_above(c(" " = 1, "lag1" = 4,"4pya" = 4),
  #                  align = "c",
  #                  strikeout = F) %>%
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "7em",
              latex_valign = "b") 

