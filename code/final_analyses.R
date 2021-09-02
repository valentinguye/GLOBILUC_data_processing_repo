
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
# 
# rm(outcome_variable, start_year, end_year, crop_j, j_soy, price_k, extra_price_k, SjPj, SkPk, fe, distribution, output, se, cluster, 
#    controls, regressors, outcome_variable)
# 
outcome_variable = "driven_loss" # "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
start_year = 2001
end_year = 2019
price_info = "lag1"
further_lu_evidence = "none"
crop_j = "Oilpalm"
j_soy = "Soybean"
fcr = 7.2
# for the k variables hypothesized in overleaf for palm oil, feglm quasipoisson converges within 25 iter.
# but maybe not with skPk controls.
price_k <- c("Soybean", "Rapeseed_oil", "Sunflower_oil",
             "Sugar", "Maize")
# c( "Sugar", "Maize")
# "Barley",  "Chicken", "Sheep", "Banana", "Beef", "Olive_oil",
#               "Orange",  "Cotton",  "Groundnut",  "Rubber", "Sorghum","Cocoa",  "Coffee",
#                "Rice",   "Wheat",  "Palm_oil", ),
#                "Tea", "Tobacco",  "Oat",
#               , "Pork")
extra_price_k = c() # , "Sheep", "Pork", "Chicken"
SjPj = TRUE
SkPk = FALSE
fe = "grid_id + country_year"
distribution = "gaussian"
se = "cluster"
cluster ="grid_id"
coefstat = "confint"
output = "coef_table"

# "Banana", "Barley", "Beef",
# "Orange", "Cocoa", "Coconut_oil", "Coffee", "Cotton", "Rice", "Groundnut",
# "Maize", "Palm_oil", "Rubber", "Sorghum", "Soybean_oil",
# "Sugar", "Tea", "Tobacco", "Wheat", "Oat", "Olive_oil", "Rapeseed_oil",
# "Sunflower_oil"
rm(outcome_variable, start_year, end_year, crop_j, j_soy, fcr, price_k, SjPj, SkPk, fe, distribution, output, se, cluster,     controls, regressors, outcome_variable)

make_reg_aeay <- function(outcome_variable = "driven_loss", # one of "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
                          start_year = 2001, 
                          end_year = 2020, 
                          further_lu_evidence = "none", # either "none", "sbqt_direct_lu", or "sbqt"mode_lu"
                          crop_j = "Oilpalm", # in GAEZ spelling
                          price_k = c("Sugar", "Maize"), # in prices spelling
                          extra_price_k = c(), # One of "Crude_oil", "Chicken", "Pork", "Sheep" 
                          price_info = "lag1", # one of "lag1", "2pya", "3pya", "4pya", "5pya",
                          fcr = 7.8, # feed conversion ratio according to Wilkinson, 2010. 
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
  var2keep <- c("grid_id", "year", "lon", "lat", "country_name", "country_year", outcome_variable, all_crop_aeay)
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
  vars <- c("grid_id", "year", "lon", "lat", "country_name", "country_year", outcome_variable, 
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
                          further_lu_evidence = "none", # either "none", "sbqt_direct_lu", or "sbqt"mode_lu"
                          crop_j = "Oilpalm", # in GAEZ spelling
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
  if(outcome_variable == "first_loss" | outcome_variable == "nd_first_loss"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aesi_long_final.Rdata"))}
  if(outcome_variable == "firstloss_glassgfc"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_aesi_long_final.Rdata"))}
  if(outcome_variable == "phtf_loss"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aesi_long_final.Rdata"))}
  if(outcome_variable == "driven_loss"){d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "driverloss_aesi_long_final.Rdata"))}
  
  # Keep only in data the useful variables 
  d <- dplyr::select(d, all_of(c("grid_id", "year", "country_name", "country_year", outcome_variable, 
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


# make_des_stats <- function(outcome_variable, 
#                            crop_j,
#                            start_year = 2002, 
#                            end_year = 2020, 
#                            qj = "rj", 
#                            price_info = "lag1")


# For each crop j, we produce two quantities, for aesi and acay. 
# Across crops, these quantities differ in the forest loss measure, and (the time and) space it is aggregated to. 

# data = d
# crops = c("Soybean", "Fodder", "Oilpalm", "Cocoa", "Coffee")
# qj = "continuous"
# price_info = "4pya"

# For acay we need to compute the revenue here. This helper function does it
make_rj <- function(data, crops, qj = "continuous", price_info = "lag1", fcr = 7.2){
  ## Revenue variable names
  # To determine land use, only the potential revenue of soybeans is considered (for simplicity) 
  all_crop_prices <- mapmat[mapmat[,"Prices"]!="Soybean_oil" & 
                              mapmat[,"Prices"]!="Soybean_meal" & 
                              mapmat[,"Prices"]!="Coconut_oil" &
                              mapmat[,"Prices"]!="Cotton", "Prices"]
  all_crop_acay <- mapmat[mapmat[,"Prices"]!="Soybean_oil" & 
                            mapmat[,"Prices"]!="Soybean_meal" & 
                            mapmat[,"Prices"]!="Coconut_oil" &
                            mapmat[,"Prices"]!="Cotton", "Crops"]
  

  ### PREPARE rj, the standardized achievable revenues
  
  # First, convert aeay quantities into TON/ha. Note that:
  # "For most crops the agro-climatic potential yield is given as kg dry weight per hectare. 
  # For alfalfa, miscanthus, switchgrass, reed canary grass, napier grass, pasture legumes and grasses the yield is given in 10 kg dry weight per hectare. 
  # For sugar beet and sugarcane (and hence Sugar, the max of them) yields are in kg sugar per hectare, 
  # and for oil palm and olives in kg oil per hectare. Cotton yield is given as kg lint per hectare." 
  # https://gaez.fao.org/pages/theme-details-theme-3
  names(data)
  # Convert those in kg/ha into ton/ha
  data <- dplyr::mutate(data, across(.cols = all_of(all_crop_acay),
                               .fns = ~./1000)) 
  # Convert Nappier grass and alfalfa (components of Fodder) from now 10 tons to tons. 
  data <- dplyr::mutate(data, across(.cols = all_of("Fodder"),
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
  # Every ton of feed counted in R_Fodder (coming from agro-climatically achievable yields currently expressed in ton/ha) is scaled to 1/fcr ton of beef meat  
  # Following Alexander et al. 2016, 25 tons of feed transform into 1 ton of beef meat (in edible weight). 
  data <- dplyr::mutate(data, Fodder = Fodder/fcr)

  # We will construct expected revenues for all crops - not just those in price_k or crop_j.
  # Thus we take the first lag or the X past year average of the prices, but not the log.
  prices_4_revenue <- paste0(all_crop_prices, "_", price_info)

  ### PREPARE rj, the standardized achievable revenues
  
  # Merge only the prices needed, not the whole price dataframe
  data <- left_join(data, prices[,c("year", prices_4_revenue)], by = c("year"))
  
  # Prices have been converted to $/t in prepare_prices.R
  for(i in 1:length(prices_4_revenue)){
    price_i <- prices_4_revenue[i]
    revenue_i <- paste0("R_", all_crop_acay[i])
    acay_i <- all_crop_acay[i]
    data <- dplyr::mutate(data, 
                       !!as.symbol(revenue_i) := !!as.symbol(acay_i)*!!as.symbol(price_i))
  }
  rm(price_i, revenue_i, acay_i)                   

  
  ## Standardize the revenue variable(s) needed here -> stream dj = Rj/sumRi 
  
  # To understand this line, see https://dplyr.tidyverse.org/articles/rowwise.html#row-wise-summary-functions
  data <- dplyr::mutate(data, Ri_sum = rowSums(across(.cols = starts_with("R_", ignore.case = FALSE)), na.rm = TRUE))
  
  if(qj == "max"){
  # Construct the max of all revenues
    data <- rowwise(data, all_of(c("grid_id", "year"))) %>% 
              dplyr::mutate(Ri_max = max(c_across(cols = starts_with("R_", ignore.case = FALSE)), na.rm = TRUE)) %>% 
              as.data.frame()
  }
  
  # Construct qj for all crops needed 
  for(crop_j in crops){
    # Explicit the name of the variable rj, the standardized potential revenue of j. 
    revenue_j <- paste0("R_", crop_j)
    revenue_j_std <- paste0(revenue_j, "_std")
    
    if(qj !="max"){
      data <- dplyr::mutate(data, !!as.symbol(revenue_j_std) := !!as.symbol(revenue_j)/Ri_sum)
    }
    if(qj == "max"){
      data <- dplyr::mutate(data, !!as.symbol(revenue_j_std) := if_else(condition = (!!as.symbol(revenue_j)==Ri_max), true = 1, false = 0))
    }
  }
  
  return(data)
}

# Infrastructure to store results
defo_table <- matrix(nrow = 8, ncol = 6, data = "")
colnames(defo_table) <- c("Crop", "GFR estimates", "AESI estimates driven loss", "AEAY estimates driven loss",
                          "AESI estimates first loss", "AEAY estimates first loss")
row.names(defo_table) <- c("Pasture Brazil", "Pasture global", "Soybean South America", "Soybean global", "Oil palm", "Cocoa", "Coffee", "Rubber")

# fill in the column of estimates from GFR (in same order as row.names !)
defo_table[,"GFR estimates"] <- c(21.9, 45.1, 8.25, 8.2, 10.5, 2.3, 1.9, 2.1)

### For soy and pasture, we use different data (prepared in south_america_track.R) and methods. 
# This is necessary because we want to compare to a measure of deforestation to soy in Goldman et al. 2020 that aggregates the whole continent.

if(!file.exists(here("temp_data", "glass_sgbp_southam")) | 
   !file.exists(here("temp_data", "glass_sgbp_brazil")) | 
   !file.exists(here("temp_data", "phtfloss_sgbp_indonesia"))){
  # First we need to prepare polygons of Brazil and South America 
  countries <- st_read(here("input_data", "Global_LSIB_Polygons_Detailed"))
  
  # Data need to be projected 
  crs_southam <- 31970
  crs_indonesia <- 23845
  # We simplify the continent but not Brazil, as the border effect on deforestation is significant (Burgess) 
  brazil <- countries[countries$COUNTRY_NA=="Brazil",] %>% st_geometry() %>% 
    st_transform(crs_southam) 
  
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
    st_transform(crs_southam) %>% 
    st_simplify(dTolerance = 10000) # simplifies by 10km
  
  indonesia <- countries[countries$COUNTRY_NA=="Indonesia",] %>% st_geometry() %>% 
    st_transform(crs_indonesia) %>% st_simplify(dTolerance = 1000)
  
  # Make filtering variable to keep only cells in Brazil or in South America. 
  # We do this for first loss (glass) aesi data but the sgbp variables are valid for the aeay dataset too. 
  # (Even if there are missing in rj that were not in sj, there is no filtering based on this.)
  d <- readRDS(here("temp_data", "merged_datasets", "southam_aoi", "glass_aesi_long_final.Rdata"))
  d_cs <- d[!duplicated(d$grid_id), c("grid_id", "lon", "lat")]  
  d_cs <- st_as_sf(d_cs, coords=c("lon", "lat"), crs = 4326)
  d_cs <- st_transform(d_cs, crs_southam) 
  
  sgbp_southam <- st_within(d_cs, southam) # rather fast thanks to simplifying
  sgbp_brazil <- st_within(d_cs, brazil) # takes ~1:30 hour
  
  saveRDS(sgbp_southam, here("temp_data", "glass_sgbp_southam"))
  saveRDS(sgbp_brazil, here("temp_data", "glass_sgbp_brazil"))
  rm(d, d_cs)
  
  # For Indonesia, we want primary forest data 
  phtfl <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aesi_long_final.Rdata"))
  phtfl_cs <- phtfl[!duplicated(phtfl$grid_id), c("grid_id", "lon", "lat")]  
  phtfl_cs <- st_as_sf(phtfl_cs, coords=c("lon", "lat"), crs = 4326)
  phtfl_cs <- st_transform(phtfl_cs, crs_indonesia) 
  
  sgbp_indonesia <- st_within(phtfl_cs, indonesia)
  saveRDS(sgbp_indonesia, here("temp_data", "phtfloss_sgbp_indonesia"))
  rm(phtfl_cs, phtfl)
  
} else {
  sgbp_southam <-readRDS(here("temp_data", "glass_sgbp_southam"))
  sgbp_brazil <-readRDS(here("temp_data", "glass_sgbp_brazil"))
  sgbp_indonesia <- readRDS(here("temp_data", "phtfloss_sgbp_indonesia"))
}

### GLASS AESI PREDICTIONS ####

### For soy and pasture, we use data prepared for the whole South American continent, not only its tropical region. 
d <- readRDS(here("temp_data", "merged_datasets", "southam_aoi", "glass_aesi_long_final.Rdata"))

# Build Y_j
d <- dplyr::select(d, grid_id, year, lon, lat, first_loss, nd_first_loss,
                   Soybean_std, Fodder_std)
d <- dplyr::mutate(d, Y_soybean = nd_first_loss*Soybean_std, 
                      Y_pasture = nd_first_loss*Fodder_std)

# add up annual deforestation over the period 2001-2015
d_0115 <- dplyr::filter(d, year >= 2001 & year <= 2015)
accu_0115 <- ddply(d_0115, "grid_id", summarise, 
                   accu_soy_defo = sum(Y_soybean, na.rm = TRUE),
                   accu_pasture_defo = sum(Y_pasture, na.rm = TRUE))

## SOY ##
# Restrict to cells precisely in South America 
accu_0115_southam <- accu_0115[lengths(sgbp_southam)>0,]
# Store total
defo_table["Soybean South America", "AESI estimates first loss"] <- round(sum(accu_0115_southam$accu_soy_defo)/1e6,digits=2) 

# add up annual soy deforestation over the period 2006-2015
d_0615 <- dplyr::filter(d, year >= 2006 & year <= 2015)
accu_soy_0615 <- ddply(d_0615, "grid_id", summarise, accu_soy_defo = sum(Y_soybean, na.rm = TRUE))
# Restrict to cells in Brazil 
accu_soy_0615_brazil <- accu_soy_0615[lengths(sgbp_brazil)>0,]
round(sum(accu_soy_0615_brazil$accu_soy_defo)/1e6,digits=2)


## PASTURE ##
# Restrict to cells in Brazil
accu_0115_brazil <- accu_0115[lengths(sgbp_brazil)>0,]
# Store total
defo_table["Pasture Brazil", "AESI estimates first loss"] <- round(sum(accu_0115_brazil$accu_pasture_defo)/1e6,digits=2) 



### For palm oil, cocoa and sugar, we aggregate over the whole tropics (i.e. on tropical_aoi dataset, and without restricting to specific countries)
d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aesi_long_final.Rdata"))

# further_lu_evidence <- "sbqt_direct_lu"
# ## restrict the outcome variable to evidence of further LU
# d <- dplyr::mutate(d, grassland = (!!as.symbol(further_lu_evidence) == 30)) # either direct or mode subsequent lu are grassland
# 
# d <- dplyr::mutate(d, cropland = (!!as.symbol(further_lu_evidence) == 10)) # either direct or mode subsequent lu are cropland
# 
# d <- dplyr::mutate(d, plantation = (sbqt_direct_lu != 20 & sbqt_mode_lu == 20)) # is not forest in the year directly after, but is forest again in the mode subsequent lu.


# Build Y_j
d <- dplyr::select(d, grid_id, year, lon, lat, first_loss, nd_first_loss, # plantation,
                   Fodder_std, Soybean_std, Oilpalm_std, Cocoa_std, Coffee_std, Rubber_std)
d <- dplyr::mutate(d, 
                   Y_pasture = nd_first_loss*Fodder_std,
                   Y_soybean = nd_first_loss*Soybean_std,
                   Y_oilpalm = nd_first_loss*Oilpalm_std,
                   Y_cocoa = nd_first_loss*Cocoa_std,
                   Y_coffee = nd_first_loss*Coffee_std, 
                   Y_rubber = nd_first_loss*Rubber_std)

# add up annual deforestation over the period 2001-2015
d_0115 <- dplyr::filter(d, year >= 2001 & year <= 2015)
# accu_0115 <- ddply(d_0115, "grid_id", summarise, 
#                    accu_oilpalm_defo = sum(Y_oilpalm, na.rm = TRUE), 
#                    accu_cocoa_defo = sum(Y_cocoa, na.rm = TRUE), 
#                    accu_coffee_defo = sum(Y_coffee, na.rm = TRUE), 
#                    accu_rubber_defo = sum(Y_rubber, na.rm = TRUE))

## PASTURE ## 
defo_table["Pasture global", "AESI estimates first loss"] <- round(sum(d_0115$Y_pasture, na.rm = TRUE)/1e6,digits=2) 

## OIL PALM ## 
defo_table["Soybean global", "AESI estimates first loss"] <- round(sum(d_0115$Y_soybean, na.rm = TRUE)/1e6,digits=2) 

## OIL PALM ## 
defo_table["Oil palm", "AESI estimates first loss"] <- round(sum(d_0115$Y_oilpalm, na.rm = TRUE)/1e6,digits=2) 

## COCOA ## 
defo_table["Cocoa", "AESI estimates first loss"] <- round(sum(d_0115$Y_cocoa, na.rm = TRUE)/1e6,digits=2) 

## COFFEE ## 
defo_table["Coffee", "AESI estimates first loss"] <- round(sum(d_0115$Y_coffee, na.rm = TRUE)/1e6,digits=2) 

## RUBBER ## 
defo_table["Rubber", "AESI estimates first loss"] <- round(sum(d_0115$Y_rubber, na.rm = TRUE)/1e6,digits=2) 

# And over Indonesia, for phtf loss data
phtfl <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aesi_long_final.Rdata"))
phtfl <- dplyr::select(phtfl, grid_id, year, lon, lat, country_name, phtf_loss, Oilpalm_std)
phtfl_idn <- dplyr::filter(phtfl, country_name == "Indonesia")
phtfl_idn <- dplyr::mutate(phtfl_idn, Y_oilpalm = phtf_loss*Oilpalm_std)
phtfl_idn_0116 <- dplyr::filter(phtfl_idn, year >= 2001 & year <= 2016)
phtfl_idn_accu_0116 <- ddply(phtfl_idn_0116, "grid_id", summarise, accu_oilpalm_defo = sum(Y_oilpalm, na.rm = TRUE))
# plot(st_geometry(phtfl_idn_accu_0116))
round(sum(phtfl_idn_accu_0116$accu_oilpalm_defo)/1e6,digits=2)


# And compute total phtf loss in 2008, to compare with Margono et al. and with Turubanova et al. (they both find 5000 km2)
phtfl_08 <- dplyr::filter(phtfl, year == 2008)
ddply(phtfl_08, "grid_id", summarise, accu_phtfl = sum(phtf_loss, na.rm = TRUE))
rm(phtfl, phtfl_08)

rm(d, d_0115, d_0615, accu_0115, accu_0115_brazil, accu_0115_southam, accu_soy_0615, accu_soy_0615_brazil)


### GLASS AEAY PREDICTIONS ####

### For soy and pasture, we use data prepared for the whole South American continent (i.e. southam_aoi dataset), not only its tropical region. 
d <- readRDS(here("temp_data", "merged_datasets", "southam_aoi", "glass_aeay_long_final.Rdata"))

# here we need to construct rj with function make_rj
d <- make_rj(data = d, 
             crops = c("Soybean", "Fodder"), 
             qj = "continuous")

# Build Y_j
d <- dplyr::select(d, grid_id, year, lon, lat, first_loss, nd_first_loss,
                   R_Soybean_std, R_Fodder_std)
d <- dplyr::mutate(d, Y_soybean = nd_first_loss*R_Soybean_std, 
                   Y_pasture = nd_first_loss*R_Fodder_std)

# add up annual deforestation over the period 2001-2015
d_0115 <- dplyr::filter(d, year >= 2001 & year <= 2015)
accu_0115 <- ddply(d_0115, "grid_id", summarise, 
                   accu_soy_defo = sum(Y_soybean, na.rm = TRUE),
                   accu_pasture_defo = sum(Y_pasture, na.rm = TRUE))


# add up annual soy deforestation over the period 2006-2015
d_0615 <- dplyr::filter(d, year >= 2006 & year <= 2015)
accu_soy_0615 <- ddply(d_0615, "grid_id", summarise, accu_soy_defo = sum(Y_soybean, na.rm = TRUE))

## SOY ##
# Restrict to cells in South America 
accu_0115_southam <- accu_0115[lengths(sgbp_southam)>0,]
# Store total
defo_table["Soybean South America", "AEAY estimates first loss"] <- round(sum(accu_0115_southam$accu_soy_defo)/1e6,digits=2) 
# Restrict to cells in Brazil 
accu_soy_0615_brazil <- accu_soy_0615[lengths(sgbp_brazil)>0,]
round(sum(accu_soy_0615_brazil$accu_soy_defo)/1e6,digits=2)


## PASTURE ##
# Restrict to cells in Brazil
accu_0115_brazil <- accu_0115[lengths(sgbp_brazil)>0,]
# Store total
defo_table["Pasture Brazil", "AEAY estimates first loss"] <- round(sum(accu_0115_brazil$accu_pasture_defo)/1e6,digits=2) 



### For the following commodities, we aggregate over the whole tropics (i.e. on tropical_aoi dataset, and without restricting to specific countries)
d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aeay_long_final.Rdata"))

# d <- dplyr::mutate(d, plantation = (sbqt_direct_lu != 20 & sbqt_mode_lu == 20)) # is not forest in the year directly after, but is forest again in the mode subsequent lu.

# here we need to construct rj with function make_rj
d <- make_rj(data = d, 
             crops = c("Fodder", "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber"), 
             qj = "continuous")

names(d)

# Build Y_j
d <- dplyr::select(d, grid_id, year, lon, lat, first_loss, nd_first_loss, # plantation,
                   R_Fodder_std, R_Soybean_std, R_Oilpalm_std, R_Cocoa_std, R_Coffee_std, R_Rubber_std)
d <- dplyr::mutate(d,
                   Y_pasture = nd_first_loss*R_Fodder_std,
                   Y_soybean = nd_first_loss*R_Soybean_std,
                   Y_oilpalm = nd_first_loss*R_Oilpalm_std,
                   Y_cocoa = nd_first_loss*R_Cocoa_std,
                   Y_coffee = nd_first_loss*R_Coffee_std,
                   Y_rubber = nd_first_loss*R_Rubber_std)

# add up annual deforestation over the period 2001-2015
d_0115 <- dplyr::filter(d, year >= 2001 & year <= 2015)
# accu_0115 <- ddply(d_0115, "grid_id", summarise, 
#                    accu_oilpalm_defo = sum(Y_oilpalm, na.rm = TRUE), 
#                    accu_cocoa_defo = sum(Y_cocoa, na.rm = TRUE), 
#                    accu_coffee_defo = sum(Y_coffee, na.rm = TRUE),
#                    accu_rubber_defo = sum(Y_rubber, na.rm = TRUE))

## PASTURE ## 
defo_table["Pasture global", "AEAY estimates first loss"] <- round(sum(d_0115$Y_pasture, na.rm = TRUE)/1e6,digits=2) 

## SOYBEAN ## 
defo_table["Soybean global", "AEAY estimates first loss"] <- round(sum(d_0115$Y_soybean, na.rm = TRUE)/1e6,digits=2) 

## OIL PALM ## 
defo_table["Oil palm", "AEAY estimates first loss"] <- round(sum(d_0115$Y_oilpalm, na.rm = TRUE)/1e6,digits=2) 

## COCOA ## 
defo_table["Cocoa", "AEAY estimates first loss"] <- round(sum(d_0115$Y_cocoa, na.rm = TRUE)/1e6,digits=2) 

## COFFEE ## 
defo_table["Coffee", "AEAY estimates first loss"] <- round(sum(d_0115$Y_coffee, na.rm = TRUE)/1e6,digits=2) 

## RUBBER ## 
defo_table["Rubber", "AEAY estimates first loss"] <- round(sum(d_0115$Y_rubber, na.rm = TRUE)/1e6,digits=2) 


# And over Indonesia, for phtf loss data
phtfl <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aeay_long_final.Rdata"))
phtfl_idn <- dplyr::filter(phtfl, country_name == "Indonesia")
# here we need to construct rj with function make_rj
phtfl_idn <- make_rj(data = phtfl_idn, 
                     crops = c("Oilpalm"), 
                     qj = "continuous")
phtfl_idn <- dplyr::select(phtfl_idn, grid_id, year, lon, lat, phtf_loss, R_Oilpalm_std)
phtfl_idn <- dplyr::mutate(phtfl_idn, Y_oilpalm = phtf_loss*R_Oilpalm_std)
phtfl_idn_0115 <- dplyr::filter(phtfl_idn, year >= 2001 & year <= 2016)
phtfl_idn_accu_0115 <- ddply(phtfl_idn_0115, "grid_id", summarise, accu_oilpalm_defo = sum(Y_oilpalm, na.rm = TRUE))
round(sum(phtfl_idn_accu_0115$accu_oilpalm_defo)/1e6,digits=2)


# ### AEAY MAX ###
# ### For soy and pasture, we use data prepared for the whole South American continent (i.e. southam_aoi dataset), not only its tropical region. 
# d <- readRDS(here("temp_data", "merged_datasets", "southam_aoi", "glass_aeay_long_final.Rdata"))
# 
# # here we need to construct rj with function make_rj
# d <- make_rj(data = d, 
#              crops = c("Soybean", "Fodder"), 
#              qj = "max") # this is the important line
# 
# # Build Y_j
# d <- dplyr::select(d, grid_id, year, lon, lat, first_loss, nd_first_loss,
#                    R_Soybean_std, R_Fodder_std)
# d <- dplyr::mutate(d, Y_soybean = nd_first_loss*R_Soybean_std, 
#                    Y_pasture = nd_first_loss*R_Fodder_std)
# 
# # add up annual deforestation over the period 2001-2015
# d_0115 <- dplyr::filter(d, year >= 2001 & year <= 2015)
# accu_0115 <- ddply(d_0115, "grid_id", summarise, 
#                    accu_soy_defo = sum(Y_soybean, na.rm = TRUE),
#                    accu_pasture_defo = sum(Y_pasture, na.rm = TRUE))
# 
# 
# # add up annual soy deforestation over the period 2006-2015
# d_0615 <- dplyr::filter(d, year >= 2006 & year <= 2015)
# accu_soy_0615 <- ddply(d_0615, "grid_id", summarise, accu_soy_defo = sum(Y_soybean, na.rm = TRUE))
# 
# 
# ## SOY ##
# # Restrict to cells in South America 
# accu_0115_southam <- accu_0115[lengths(sgbp_southam)>0,]
# # Store total
# round(sum(accu_0115_southam$accu_soy_defo)/1e6,digits=2) 
# # Restrict to cells in Brazil in 2006-2015
# accu_soy_0615_brazil <- accu_soy_0615[lengths(sgbp_brazil)>0,]
# round(sum(accu_soy_0615_brazil$accu_soy_defo)/1e6,digits=2)
# 
# 
# ## PASTURE ##
# # Restrict to cells in Brazil
# accu_0115_brazil <- accu_0115[lengths(sgbp_brazil)>0,]
# # Store total
# round(sum(accu_0115_brazil$accu_pasture_defo)/1e6,digits=2) 
# 
# 
# 
# ### For the following commodities, we aggregate over the whole tropics (i.e. on tropical_aoi dataset, and without restricting to specific countries)
# d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aeay_long_final.Rdata"))
# # here we need to construct rj with function make_rj
# d <- make_rj(data = d, 
#              crops = c("Oilpalm", "Cocoa", "Coffee"), 
#              qj = "max") # this is the important line
# # Build Y_j
# d <- dplyr::select(d, grid_id, year, lon, lat, first_loss, nd_first_loss,
#                    R_Oilpalm_std, R_Cocoa_std, R_Coffee_std)
# d <- dplyr::mutate(d, Y_oilpalm = nd_first_loss*R_Oilpalm_std,
#                    Y_cocoa = nd_first_loss*R_Cocoa_std,
#                    Y_coffee = nd_first_loss*R_Coffee_std)
# 
# # add up annual deforestation over the period 2001-2015
# d_0115 <- dplyr::filter(d, year >= 2001 & year <= 2015)
# accu_0115 <- ddply(d_0115, "grid_id", summarise, 
#                    accu_oilpalm_defo = sum(Y_oilpalm, na.rm = TRUE), 
#                    accu_cocoa_defo = sum(Y_cocoa, na.rm = TRUE), 
#                    accu_coffee_defo = sum(Y_coffee, na.rm = TRUE))
# 
# ## OIL PALM ## 
# round(sum(accu_0115$accu_oilpalm_defo)/1e6,digits=2) 
# 
# # And over Indonesia, for phtf loss data
# phtfl <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aeay_long_final.Rdata"))
# phtfl_idn <- dplyr::filter(phtfl, country_name == "Indonesia")
# # here we need to construct rj with function make_rj
# phtfl_idn <- make_rj(data = phtfl_idn, 
#                      crops = c("Oilpalm"), 
#                      qj = "max")
# phtfl_idn <- dplyr::select(phtfl_idn, grid_id, year, lon, lat, phtf_loss, R_Oilpalm_std)
# phtfl_idn <- dplyr::mutate(phtfl_idn, Y_oilpalm = phtf_loss*R_Oilpalm_std)
# phtfl_idn_0115 <- dplyr::filter(phtfl_idn, year >= 2001 & year <= 2016)
# phtfl_idn_accu_0115 <- ddply(phtfl_idn_0115, "grid_id", summarise, accu_oilpalm_defo = sum(Y_oilpalm, na.rm = TRUE))
# round(sum(phtfl_idn_accu_0115$accu_oilpalm_defo)/1e6,digits=2)
# 
# 
# ## COCOA ## 
# round(sum(accu_0115$accu_cocoa_defo)/1e6,digits=2) 
# 
# 
# ## COFFEE ## 
# round(sum(accu_0115$accu_coffee_defo)/1e6,digits=2) 









### DRIVEN LOSS AESI PREDICTIONS #### 
# For soy and pasture, we use data prepared for the whole South American continent, not only its tropical region. 
d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "driverloss_aesi_long_final.Rdata"))

# Build Y_j
d <- dplyr::select(d, grid_id, year, lon, lat, driven_loss, country_name,
                   Soybean_std, Fodder_std, Oilpalm_std, Cocoa_std, Coffee_std, Rubber_std)
d <- dplyr::mutate(d, Y_soybean = driven_loss*Soybean_std, 
                   Y_pasture = driven_loss*Fodder_std,
                   Y_oilpalm = driven_loss*Oilpalm_std,
                   Y_cocoa = driven_loss*Cocoa_std,
                   Y_coffee = driven_loss*Coffee_std, 
                   Y_rubber = driven_loss*Rubber_std)

# add up annual deforestation over the period 2001-2015
d_0115 <- dplyr::filter(d, year >= 2001 & year <= 2015)

## SOY ##
defo_table["Soybean global", "AESI estimates driven loss"] <- round(sum(d_0115$Y_soybean, na.rm = TRUE)/1e6,digits=2) 

# Restrict to cells precisely in South America 
southam_countries <- c("Argentina", 
                       "Bolivia", 
                       "Brazil",  
                       "Isla Brasilera (disp)",  
                       "Chile",    
                       "Colombia",  
                       "Ecuador",  
                       "French Guiana (Fr)", 
                       "Guyana", 
                       "Panama", 
                       "Paraguay", 
                       "Peru", 
                       "Suriname",  
                       "Trinidad & Tobago",  
                       "Uruguay",  
                       "Venezuela")
d_0115_southam <- dplyr::filter(d_0115, country_name %in% southam_countries)

# Store total
defo_table["Soybean South America", "AESI estimates driven loss"] <- round(sum(d_0115_southam$Y_soybean, na.rm = TRUE)/1e6,digits=2) 

# add up annual soy deforestation over the period 2006-2015
d_0615 <- dplyr::filter(d, year >= 2006 & year <= 2015)
# Restrict to cells in Brazil 
d_0615_brazil <- d_0615[d_0615$country_name == "Brazil" | d_0615$country_name == "Isla Brasilera (disp)",]
round(sum(d_0615_brazil$Y_soybean, na.rm = TRUE)/1e6,digits=2)


## PASTURE ##
defo_table["Pasture global", "AESI estimates driven loss"] <- round(sum(d_0115$Y_pasture, na.rm = TRUE)/1e6,digits=2)# Store total

# Restrict to cells in Brazil
d_0115_brazil <- d_0115[d_0115$country_name == "Brazil" | d_0115$country_name == "Isla Brasilera (disp)",]

defo_table["Pasture Brazil", "AESI estimates driven loss"] <- round(sum(d_0115_brazil$Y_pasture, na.rm = TRUE)/1e6,digits=2)# Store total


## OIL PALM ## 
defo_table["Oil palm", "AESI estimates driven loss"] <- round(sum(d_0115_brazil$Y_oilpalm, na.rm = TRUE)/1e6,digits=2) 

## COCOA ## 
defo_table["Cocoa", "AESI estimates driven loss"] <- round(sum(d_0115_brazil$Y_cocoa, na.rm = TRUE)/1e6,digits=2)          

## COFFEE ## 
defo_table["Coffee", "AESI estimates driven loss"] <- round(sum(d_0115_brazil$Y_coffee, na.rm = TRUE)/1e6,digits=2) 

## RUBBER ## 
defo_table["Rubber", "AESI estimates driven loss"] <- round(sum(d_0115_brazil$Y_rubber, na.rm = TRUE)/1e6,digits=2)



### DRIVEN LOSS AEAY PREDICTIONS #### 
d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "driverloss_aeay_long_final.Rdata"))

# here we need to construct rj with function make_rj
d <- make_rj(data = d, 
             crops = c("Fodder", "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber"), 
             qj = "continuous")

s# Build Y_j
d <- dplyr::select(d, grid_id, year, lon, lat, driven_loss, country_name,
                   R_Soybean_std, R_Fodder_std, R_Oilpalm_std, R_Cocoa_std, R_Coffee_std, R_Rubber_std)

d <- dplyr::mutate(d, Y_soybean = driven_loss*R_Soybean_std, 
                   Y_pasture = driven_loss*R_Fodder_std,
                   Y_oilpalm = driven_loss*R_Oilpalm_std,
                   Y_cocoa = driven_loss*R_Cocoa_std,
                   Y_coffee = driven_loss*R_Coffee_std, 
                   Y_rubber = driven_loss*R_Rubber_std)

# add up annual deforestation over the period 2001-2015
d_0115 <- dplyr::filter(d, year >= 2001 & year <= 2015)

## SOY ##
defo_table["Soybean global", "AEAY estimates driven loss"] <- round(sum(d_0115$Y_soybean, na.rm = TRUE)/1e6,digits=2) 

# Restrict to cells precisely in South America 
southam_countries <- c("Argentina", 
                       "Bolivia", 
                       "Brazil",  
                       "Isla Brasilera (disp)",  
                       "Chile",    
                       "Colombia",  
                       "Ecuador",  
                       "French Guiana (Fr)", 
                       "Guyana", 
                       "Panama", 
                       "Paraguay", 
                       "Peru", 
                       "Suriname",  
                       "Trinidad & Tobago",  
                       "Uruguay",  
                       "Venezuela")
d_0115_southam <- dplyr::filter(d_0115, country_name %in% southam_countries)

# Store total
defo_table["Soybean South America", "AEAY estimates driven loss"] <- round(sum(d_0115_southam$Y_soybean, na.rm = TRUE)/1e6,digits=2) 

## Count soy deforestation over the period 2006-2015 in Brazil
d_0615 <- dplyr::filter(d, year >= 2006 & year <= 2015)
d_0615_brazil <- d_0615[d_0615$country_name == "Brazil" | d_0615$country_name == "Isla Brasilera (disp)",]
round(sum(d_0615_brazil$Y_soybean, na.rm = TRUE)/1e6,digits=2)


## PASTURE ##
defo_table["Pasture global", "AEAY estimates driven loss"] <- round(sum(d_0115$Y_pasture, na.rm = TRUE)/1e6,digits=2)# Store total

# Restrict to cells in Brazil
d_0115_brazil <- d_0115[d_0115$country_name == "Brazil" | d_0115$country_name == "Isla Brasilera (disp)",]

defo_table["Pasture Brazil", "AEAY estimates driven loss"] <- round(sum(d_0115_brazil$Y_pasture, na.rm = TRUE)/1e6,digits=2)# Store total


## OIL PALM ## 
defo_table["Oil palm", "AEAY estimates driven loss"] <- round(sum(d_0115_brazil$Y_oilpalm, na.rm = TRUE)/1e6,digits=2) 

## COCOA ## 
defo_table["Cocoa", "AEAY estimates driven loss"] <- round(sum(d_0115_brazil$Y_cocoa, na.rm = TRUE)/1e6,digits=2)          

## COFFEE ## 
defo_table["Coffee", "AEAY estimates driven loss"] <- round(sum(d_0115_brazil$Y_coffee, na.rm = TRUE)/1e6,digits=2) 

## RUBBER ## 
defo_table["Rubber", "AEAY estimates driven loss"] <- round(sum(d_0115_brazil$Y_rubber, na.rm = TRUE)/1e6,digits=2)

#### EDIT LATEX TABLE #### 
colnames(defo_table) <- NULL

options(knitr.table.format = "latex")
kable(defo_table, booktabs = T, align = "c",
      caption = "Predicted deforestation for agricultural land use, 2001-2015 (Mha)") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                    "GFR estimates" = 1,
                    "AESI" = 1,
                    "AEAY" = 1,
                    "AESI" = 1,
                    "AEAY" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 2,
                     "Forest lost to agriculture" = 2,
                     "First forest loss" = 2),
                   bold = T,
                   align = "c") %>%
  column_spec(column = 1,
              width = "8em",
              latex_valign = "b") %>%
  column_spec(column = 2:ncol(defo_table),
              width = "4em",
              latex_valign = "b") 



#### DESCRIPTIVE MAP AESI ####
d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aesi_long_final.Rdata"))

# Build Y_j
d <- dplyr::mutate(d, 
                   Y_pasture = first_loss*Fodder_std,
                   Y_soybean = first_loss*Soybean_std,
                   Y_oilpalm = first_loss*Oilpalm_std,
                   Y_cocoa = first_loss*Cocoa_std,
                   Y_coffee = first_loss*Coffee_std, 
                   Y_rubber = first_loss*Rubber_std)

# add up annual deforestation over the period 2001-2015
d_0115 <- dplyr::filter(d, year >= 2001 & year <= 2015)
accu_0115 <- ddply(d_0115, "grid_id", summarise, 
                   accu_pasture_defo = sum(Y_pasture, na.rm = TRUE), 
                   accu_soybean_defo = sum(Y_soybean, na.rm = TRUE), 
                   accu_oilpalm_defo = sum(Y_oilpalm, na.rm = TRUE), 
                   accu_cocoa_defo = sum(Y_cocoa, na.rm = TRUE), 
                   accu_coffee_defo = sum(Y_coffee, na.rm = TRUE), 
                   accu_rubber_defo = sum(Y_rubber, na.rm = TRUE))


# identify the main driver and its imputed deforestation
accu_0115 <- accu_0115 %>% rowwise() %>% 
                        mutate(accu_main_driver = max(c(accu_pasture_defo,
                                                        accu_soybean_defo, 
                                                        accu_oilpalm_defo,
                                                        accu_cocoa_defo,
                                                        accu_coffee_defo,
                                                        accu_rubber_defo))) %>% 
                        as.data.frame()

accu_vars <- c("accu_pasture_defo",
               "accu_soybean_defo", 
               "accu_oilpalm_defo",
               "accu_cocoa_defo",
               "accu_coffee_defo",
               "accu_rubber_defo")
accu_0115$main_driver <- ""
for(i in 1:nrow(accu_0115)){
  accu_vec <- accu_0115[i, accu_vars]
  if(length(colnames(accu_vec)[accu_vec == accu_0115[i,"accu_main_driver"]])==1){
    accu_0115[i,"main_driver"] <- colnames(accu_vec)[accu_vec == accu_0115[i,"accu_main_driver"]]
  }else(accu_0115[i,"main_driver"] <- "none")
}
# accu_0115 <- accu_0115[!duplicated(accu_0115$grid_id),]
# accu_0115 <- accu_0115[,-c("lon", "lat")]  
d_cs <- d[!duplicated(d$grid_id),]
accu_0115 <- left_join(accu_0115, d_cs[,c("grid_id", "lon", "lat")], by = "grid_id")

accu_0115 <- st_as_sf(accu_0115, coords = c("lon", "lat"), crs = 4326) 
mercator_world_crs <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "
accu_0115 <- st_transform(accu_0115, crs = mercator_world_crs)

accu_0115 <- st_buffer(accu_0115, dist = 4000)
st_geometry(accu_0115) <- sapply(st_geometry(accu_0115), FUN = function(x){st_as_sfc(st_bbox(x))}) %>% st_sfc(crs = mercator_world_crs)


tm_shape(accu_0115[accu_0115$main_driver!="none",]) +
  tm_borders(alpha = 0) + 
  tm_fill(col = "main_driver", palette = rainbow(n = 6))

sum(accu_0115[accu_0115$main_driver == "accu_pasture_defo",]$accu_main_driver)
sum(accu_0115[accu_0115$main_driver == "accu_soybean_defo",]$accu_main_driver)
sum(accu_0115[accu_0115$main_driver == "accu_oilpalm_defo",]$accu_main_driver)
sum(accu_0115[accu_0115$main_driver == "accu_cocoa_defo",]$accu_main_driver)
sum(accu_0115[accu_0115$main_driver == "accu_coffee_defo",]$accu_main_driver)
sum(accu_0115[accu_0115$main_driver == "accu_rubber_defo",]$accu_main_driver)

nrow(accu_0115[accu_0115$main_driver=="none",])

sum(accu_0115$accu_main_driver)
sum(d$first_loss)


summary(df_cs$Fodder_std)
summary(df_cs$Soybean_std)
summary(df_cs$Oilpalm_std)
summary(df_cs$Cocoa_std)
summary(df_cs$Coffee_std)
summary(df_cs$Rubber_std)

plot(d_cs[,"Fodder_std"])

plot(d_cs[,"Coffee_std"])

summary(d$first_loss) # first loss expressed in ha 
summary(d$Soybean_std) # _std vars expressed in ratio (0:1)

d <- dplyr::mutate(d, 
                   Y_soybean = first_loss*Soybean_std, 
                   Y_pasture = first_loss*Fodder_std, 
                   Y_oilpalm = first_loss*Oilpalm_std, 
                   Y_coffee = first_loss*Coffee_std, 
                   Y_cocoa = first_loss*Cocoa_std,
                   Y_rubber = first_loss*Rubber_std)

# over 1983-2015
summary(d$Y_pasture)
summary(d$Y_soybean)
summary(d$Y_oilpalm)
summary(d$Y_cocoa)
summary(d$Y_coffee)
summary(d$Y_rubber)

# 0ver 2001-2015
d0115 <- d[d$year>=2001 & d$year <=2015,] 
summary(d0115$Y_pasture)
summary(d0115$Y_soybean)
summary(d0115$Y_oilpalm)
summary(d0115$Y_cocoa)
summary(d0115$Y_coffee)
summary(d0115$Y_rubber)

sum(d0115$Y_pasture, na.rm = TRUE)/1e6
sum(d0115$Y_soybean, na.rm = TRUE)/1e6
sum(d0115$Y_oilpalm, na.rm = TRUE)/1e6
sum(d0115$Y_cocoa, na.rm = TRUE)/1e6
sum(d0115$Y_coffee, na.rm = TRUE)/1e6
sum(d0115$Y_rubber, na.rm = TRUE)/1e6

d <- readRDS(here("temp_data", "merged_datasets", "southam_aoi", "glass_aesi_long_final.Rdata"))
d0115 <- d[d$year>=2001 & d$year <=2015,] 
d0115[lengths(sgbp_southam)>0,]


#### DESCRIPTIVE MAP AEAY ####
d <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aeay_long_final.Rdata"))

# here we need to construct rj with function make_rj
d <- make_rj(data = d, 
             crops = c("Soybean", "Fodder", "Oilpalm", "Cocoa", "Coffee", "Rubber"), 
             qj = "continuous", 
             fcr = 6)

# Build Y_j
d <- dplyr::mutate(d, 
                   Y_pasture = first_loss*R_Fodder_std,
                   Y_soybean = first_loss*R_Soybean_std,
                   Y_oilpalm = first_loss*R_Oilpalm_std,
                   Y_cocoa = first_loss*R_Cocoa_std,
                   Y_coffee = first_loss*R_Coffee_std, 
                   Y_rubber = first_loss*R_Rubber_std)

# add up annual deforestation over the period 2001-2015
d_0115 <- dplyr::filter(d, year >= 2001 & year <= 2015)
accu_0115 <- ddply(d_0115, "grid_id", summarise, 
                   accu_pasture_defo = sum(Y_pasture, na.rm = TRUE), 
                   accu_soybean_defo = sum(Y_soybean, na.rm = TRUE), 
                   accu_oilpalm_defo = sum(Y_oilpalm, na.rm = TRUE), 
                   accu_cocoa_defo = sum(Y_cocoa, na.rm = TRUE), 
                   accu_coffee_defo = sum(Y_coffee, na.rm = TRUE), 
                   accu_rubber_defo = sum(Y_rubber, na.rm = TRUE))


# identify the main driver and its imputed deforestation
accu_0115 <- accu_0115 %>% rowwise() %>% 
  mutate(accu_main_driver = max(c(accu_pasture_defo,
                                  accu_soybean_defo, 
                                  accu_oilpalm_defo,
                                  accu_cocoa_defo,
                                  accu_coffee_defo,
                                  accu_rubber_defo))) %>% 
  as.data.frame()

accu_vars <- c("accu_pasture_defo",
               "accu_soybean_defo", 
               "accu_oilpalm_defo",
               "accu_cocoa_defo",
               "accu_coffee_defo",
               "accu_rubber_defo")
accu_0115$main_driver <- ""
for(i in 1:nrow(accu_0115)){
  accu_vec <- accu_0115[i, accu_vars]
  if(length(colnames(accu_vec)[accu_vec == accu_0115[i,"accu_main_driver"]])==1){
    accu_0115[i,"main_driver"] <- colnames(accu_vec)[accu_vec == accu_0115[i,"accu_main_driver"]]
  }else(accu_0115[i,"main_driver"] <- "none")
}
# accu_0115 <- accu_0115[!duplicated(accu_0115$grid_id),]
# accu_0115 <- accu_0115[,-c("lon", "lat")]  
d_cs <- d[!duplicated(d$grid_id),]
accu_0115 <- left_join(accu_0115, d_cs[,c("grid_id", "lon", "lat")], by = "grid_id")

accu_0115 <- st_as_sf(accu_0115, coords = c("lon", "lat"), crs = 4326) 
mercator_world_crs <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "
accu_0115 <- st_transform(accu_0115, crs = mercator_world_crs)

accu_0115 <- st_buffer(accu_0115, dist = 4000)
st_geometry(accu_0115) <- sapply(st_geometry(accu_0115), FUN = function(x){st_as_sfc(st_bbox(x))}) %>% st_sfc(crs = mercator_world_crs)


tm_shape(accu_0115[accu_0115$main_driver!="none",]) +
  tm_borders(alpha = 0) + 
  tm_fill(col = "main_driver", palette = rainbow(n = 6), 
          convert2density = TRUE, area = "accu_main_driver")  


#### CONSTRUCT SPECIFICATION COMPARATIVE TABLES ####

# We build 6 tables, one for each j crop. 
# These tables compare aeay (rj) and aesi (sj) models, with and without controlling for SkPk, with different price info  

price_infoS <- c("lag1", "4pya")

# Forest definition
SY <- 2001
EY <- 2019
# with controls or not 
control_ornot <- c(FALSE, TRUE)


### BEEF
K_beef <- c("Soybean", 
            "Maize", "Wheat", "Barley", "Oat") # in prices spelling "cereal_crops" "Barley", "Oat", "Sorghum", "Maize", "Rice"

res_list_beef <- list()
elm <- 1

for(PI in price_infoS){
  for(CTRL in control_ornot){
    res_list_beef[[elm]] <- make_reg_aeay(price_info = PI,
                                               start_year = 2001, end_year = 2019,
                                               crop_j = "Fodder",
                                               price_k = K_beef,
                                               SkPk= CTRL,
                                               extra_price_k = c())#"Sheep", "Pork", "Chicken"
    names(res_list_beef)[elm] <- paste0(PI,"_",CTRL, "_acay")
    elm <- elm + 1
    
    res_list_beef[[elm]] <- make_reg_aesi(price_info = PI,
                                               start_year = SY, end_year = EY,
                                               crop_j = "Fodder",
                                               price_k = K_beef,
                                               SkPk= CTRL,
                                               extra_price_k = c())#"Sheep", "Pork", "Chicken"
    names(res_list_beef)[elm] <- paste0(PI,"_",CTRL, "_aesi")
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
      caption = paste0("Indirect effects of global commodity markets on deforestation for cattle, ",SY,"-",EY)) %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Estimates" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1, "lag1" = 4,"4pya" = 4),
                   align = "c",
                   strikeout = F) %>%
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

for(PI in price_infoS){
  for(CTRL in control_ornot){
    res_list_oilpalm[[elm]] <- make_reg_aeay(price_info = PI,
                                               start_year = SY, end_year = EY,
                                               crop_j = "Oilpalm",
                                               price_k = K_oilpalm,
                                               SkPk= CTRL,
                                               extra_price_k = c())#
    names(res_list_oilpalm)[elm] <- paste0(PI,"_",CTRL, "_acay")
    elm <- elm + 1
    
    res_list_oilpalm[[elm]] <- make_reg_aesi(price_info = PI,
                                               start_year = SY, end_year = EY,
                                               crop_j = "Oilpalm",
                                               price_k = K_oilpalm,
                                               SkPk= CTRL,
                                               extra_price_k = c())#
    names(res_list_oilpalm)[elm] <- paste0(PI,"_",CTRL, "_aesi")
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
      caption = paste0("Indirect effects of global commodity markets on deforestation for oil palm, ",SY,"-",EY)) %>% #
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Estimates" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1, "lag1" = 4,"4pya" = 4),
                   align = "c",
                   strikeout = F) %>%
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

for(PI in price_infoS){
  for(CTRL in control_ornot){
    res_list_soy[[elm]] <- make_reg_aeay(price_info = PI,
                                             start_year = SY, end_year = EY,
                                             crop_j = "Soybean",
                                             price_k = K_soy,
                                             SkPk= CTRL,
                                             extra_price_k = c())#
    names(res_list_soy)[elm] <- paste0(PI,"_",CTRL, "_acay")
    elm <- elm + 1
    
    res_list_soy[[elm]] <- make_reg_aesi(price_info = PI,
                                             start_year = SY, end_year = EY,
                                             crop_j = "Soybean",
                                             price_k = K_soy,
                                             SkPk= CTRL,
                                             extra_price_k = c())#
    names(res_list_soy)[elm] <- paste0(PI,"_",CTRL, "_aesi")
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
      caption = paste0("Indirect effects of global commodity markets on deforestation for soy, ",SY,"-",EY)) %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Estimates" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1, "lag1" = 4,"4pya" = 4),
                   align = "c",
                   strikeout = F) %>%
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

for(PI in price_infoS){
  for(CTRL in control_ornot){
    res_list_cocoa[[elm]] <- make_reg_aeay(price_info = PI,
                                             start_year = SY, end_year = EY,
                                             crop_j = "Cocoa",
                                             price_k = K_cocoa,
                                             SkPk= CTRL,
                                             extra_price_k = c())#
    names(res_list_cocoa)[elm] <- paste0(PI,"_",CTRL, "_acay")
    elm <- elm + 1
    
    res_list_cocoa[[elm]] <- make_reg_aesi(price_info = PI,
                                             start_year = SY, end_year = EY,
                                             crop_j = "Cocoa",
                                             price_k = K_cocoa,
                                             SkPk= CTRL,
                                             extra_price_k = c())#
    names(res_list_cocoa)[elm] <- paste0(PI,"_",CTRL, "_aesi")
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
      caption = paste0("Indirect effects of global commodity markets on deforestation for cocoa, ",SY,"-",EY)) %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Estimates" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1, "lag1" = 4,"4pya" = 4),
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

for(PI in price_infoS){
  for(CTRL in control_ornot){
    res_list_coffee[[elm]] <- make_reg_aeay(price_info = PI,
                                             start_year = SY, end_year = EY,
                                             crop_j = "Coffee",
                                             price_k = K_coffee,
                                             SkPk= CTRL,
                                             extra_price_k = c())#
    names(res_list_coffee)[elm] <- paste0(PI,"_",CTRL, "_acay")
    elm <- elm + 1
    
    res_list_coffee[[elm]] <- make_reg_aesi(price_info = PI,
                                             start_year = SY, end_year = EY,
                                             crop_j = "Coffee",
                                             price_k = K_coffee,
                                             SkPk= CTRL,
                                             extra_price_k = c())#
    names(res_list_coffee)[elm] <- paste0(PI,"_",CTRL, "_aesi")
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
      caption = paste0("Indirect effects of global commodity markets on deforestation for coffee, ",SY,"-",EY)) %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Estimates" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1, "lag1" = 4,"4pya" = 4),
                   align = "c",
                   strikeout = F) %>%
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

for(PI in price_infoS){
  for(CTRL in control_ornot){
    res_list_coffee[[elm]] <- make_reg_aeay(price_info = PI,
                                            start_year = SY, end_year = EY,
                                            crop_j = "Rubber",
                                            price_k = K_rubber,
                                            SkPk= CTRL,
                                            extra_price_k = c())#
    names(res_list_rubber)[elm] <- paste0(PI,"_",CTRL, "_acay")
    elm <- elm + 1
    
    res_list_coffee[[elm]] <- make_reg_aesi(price_info = PI,
                                            start_year = SY, end_year = EY,
                                            crop_j = "Rubber",
                                            price_k = K_rubber,
                                            SkPk= CTRL,
                                            extra_price_k = c())#
    names(res_list_rubber)[elm] <- paste0(PI,"_",CTRL, "_aesi")
    elm <- elm + 1
  }
}

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
  add_header_above(c("Estimates" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1,"SI" = 1,"PR" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1, "lag1" = 4,"4pya" = 4),
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
K_beef <- c("Soybean", "Soybean_meal", "Palm_oil", 
             "cereal_crops") # in prices spelling "Barley", "Oat", "Sorghum", "Maize", "Rice"

beef_res <- make_reg_aeay(outcome_variable = OV,
                             start_year = SY, end_year = EY,
                             crop_j = "Fodder",
                             price_k = K_beef,
                          SkPk= F,
                             extra_price_k = c())#"Sheep", "Pork", "Chicken"

### PALM OIL
K_palmoil <- c("Soybean_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil",
               "Sugar", "Maize") # in prices spelling "Olive_oil", "Soybean", "Soybean_meal",

oilpalm_res <- make_reg_aeay(outcome_variable = OV,
                             start_year = SY, end_year = EY,
                             crop_j = "Oilpalm",
                             price_k = K_palmoil,
                             extra_price_k = "Crude_oil")


### SOY 
K_soy <- c("Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", 
           "Sugar", "Maize") # in prices spelling "Olive_oil",


soy_res <- make_reg_aeay(outcome_variable = OV,
                         start_year = SY, end_year = EY, 
                         crop_j = "Soybean", 
                         price_k = K_soy, 
                         extra_price_k = "Crude_oil")


### COCOA 
K_cocoa <- c("Coffee", "Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", "Sugar") # in prices spelling


cocoa_res <- make_reg_aeay(outcome_variable = OV,
                           start_year = SY, end_year = EY, 
                           crop_j = "Cocoa", 
                           price_k = K_cocoa)


### COFFEE 
K_coffee <- c("Tea", "Cocoa", "Sugar", "Tobacco") # in prices spelling


coffee_res <- make_reg_aeay(outcome_variable = OV,
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



#### RUN AND PLOT AESI ####
DS <- "glass"
SY <- 1983
EY <- 2020
if(DS == "glass"){main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aesi_long_final.Rdata"))}
if(DS == "fl8320"){main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_aesi_long_final.Rdata"))}
if(DS == "phtfl"){main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_aesi_long_final.Rdata"))}


### PALM OIL
K_palmoil <- c("Soybean_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil",
               "Sugar", "Maize") # in prices spelling "Olive_oil", "Soybean", "Soybean_meal",

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



#### First loss (1983-2015) ####
main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_aesi_long_final.Rdata"))

### PALM OIL
K_palmoil <- c("Soybean_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil",
               "Sugar", "Maize") # in prices spelling "Olive_oil", "Soybean", "Soybean_meal",

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
                         Soybean = "Soybean", 
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
               "Sugar", "Maize") # in prices spelling "Olive_oil", "Soybean", "Soybean_meal",

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
                         Soybean = "Soybean", 
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



