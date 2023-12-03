##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
neededPackages <- c("data.table", "plyr", "tidyr", "dplyr",  "Hmisc", "sjmisc", "stringr",
                    "here", "readstata13", "foreign", "readxl", "writexl",
                    "raster", "sp", "spdep", "sf","gfcanalysis",  "nngeo", "stars", # "osrm", "osrmr",
                    "lubridate","exactextractr",
                    "doParallel", "foreach", "snow", "parallel",
                    "knitr", "kableExtra",
                    "DataCombine", 
                    "fixest", 
                    "boot", "fwildclusterboot", "sandwich", "MASS",
                    "ggplot2", "leaflet", "tmap",  "dotwhisker", "viridis", "hrbrthemes")
# "pglm", "multiwayvcov", "clusterSEs", "alpaca", "clubSandwich",

# Install them in their project-specific versions
renv::restore(packages = neededPackages)

# Often useful to upgrade renv to debug: 
# renv::upgrade()

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
dir.create(here("temp_data","reg_results", "rfs", "result_lists"))

### GLOBAL CRS USED ### 
mercator_world_crs <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "


### LABEL FOR PLOTS 
# relabel terms (for y axis)
# row.names(df) <- NULL
# df <- mutate(df, term = paste0(crop,".",Dynamics))
# predictors_dict <- rep("", length(df$term)) 
# names(predictors_dict) <- df$term
# predictors_dict[grepl("Cumulative", names(predictors_dict))] <- df$crop[grepl("Cumulative", names(predictors_dict))]
# df <- relabel_predictors(df,predictors_dict) 
predictors_dict <- c(eaear_Maizegrain = "Maize",
                     eaear_Cereals = "Cereals",
                     eaear_Rice = "Rice",
                     eaear_Roots = "Roots crops",
                     eaear_Soy_compo = "Soy",
                     eaear_Fodder = "Pasture crops",
                     eaear_Cotton = "Cotton",
                     eaear_Oilfeed_crops = "Other oil crops",
                     eaear_Biomass = "Biomass crops",
                     eaear_Sugarcane = "Sugarcane",
                     eaear_Tobacco = "Tobacco",
                     eaear_Citrus = "Citrus",
                     eaear_Oilpalm = "Oil palm",
                     eaear_Coconut = "Coconut",
                     eaear_Banana = "Banana",
                     eaear_Cocoa_Coffee = "Cocoa or Coffee",
                     eaear_Tea = "Tea",
                     eaear_Rubber = "Rubber" )


### OBJECTS USED IN RFS PROCESSES ###

### RFS DATA #### 
prices <- readRDS(here("temp_data", "prepared_international_prices.Rdata"))

# renewable fuel standards, as from https://www.epa.gov/renewable-fuel-standard-program/renewable-fuel-annual-standards
# remember: biodiesel are nested within advanced biofuels, which are themeselves nested within total
# conventional biofuels are the part of the total that is not "advanced"
rfs <- data.frame(year = 1994:2022, # taking from 1994 just so that even with 6-year lag there is no NA removing year 2000, which may be used as the first year with undernourishment outcome for the pre-treatment period. 
                  statute_total = 0, final_total = 0, 
                  statute_advanced = 0, final_advanced = 0, 
                  statute_biodiesel = 0, final_biodiesel = 0, 
                  statute_conv_earliest = 0)

rfs <- dplyr::arrange(rfs, year)
# this is if RFS2 takes over as soon as 2008 (9bgal). 
# The two first years of NAs represent the fact that these mandates are endogenous and should thus not be used 
# /!\ THIS IS DIFFERENT FROM THE ILUC PROJECT /!\
rfs[rfs$year >= 2006 & rfs$year <= 2022, c("statute_total")] <- c(NA, NA, 9, 11.1, 12.95,	13.95,	15.2,	16.55,	18.15,	20.5,	22.25,	24.0,	26.0,	28.0, 30.0, 33.0, 36.0)
rfs[rfs$year >= 2006 & rfs$year <= 2022, c("final_total")] <- c(NA, NA, 9, 11.1, 12.95,	13.95,	15.2,	16.55,	16.28,	16.93,	18.11,	19.28,	19.29,	19.92, 20.09, NA, NA)

rfs[rfs$year >= 2006 & rfs$year <= 2022, c("statute_advanced")] <- c(NA, NA, 0, 0.6, 0.95,	1.35,	2.0,	2.75,	3.75,	5.5,	7.25,	9.0,	11.0,	13.0, 15.0, 18.0, 21.0)
rfs[rfs$year >= 2006 & rfs$year <= 2022, c("final_advanced")] <- c(NA, NA, 0, 0.6, 0.95,	1.35,	2.0,	2.75,	2.67,	2.88,	3.61,	4.28,	4.29,	4.92, 5.09, NA, NA)

rfs[rfs$year >= 2009 & rfs$year <= 2022, c("statute_biodiesel")] <- c(0.5, 0.65, 0.8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
rfs[rfs$year >= 2009 & rfs$year <= 2022, c("final_biodiesel")] <- c(0.5, 1.15, 0.8, 1.0, 1.28, 1.63, 1.73, 1.9, 2.0, 2.1, 2.1, 2.43, 2.43, NA)

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
  for(lead in c(1:9)){
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
  
  for(lag in c(1:9)){
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
all_rfs_treatments <- grep(pattern = "statute_conv", names(prices), value = TRUE)



#### REGRESSION FUNCTION #### 
pre_process <- FALSE
# pre_processed_data <- pre_d_clean_agri

# THESE GOVERNS WHETHER IT'S A JOINT OR A DISJOINT ESTIMATION 
# for joint: all crops passed to exposure_rfs, and control_all_absolute_rfs = FALSE
exposure_rfs <- "eaear_Soy_compo"
control_all_absolute_rfs <- TRUE # whether to control for other crops
annual_rfs_controls <- TRUE # on annual or year averaged dynamics
all_exposures_rfs = crops_ctrl # which crops

trade_exposure = "trade_expo"
trade_exposure_period = "20012007"

# treatment dynamics
original_rfs_treatments <- c("statute_conv")
# These 4 arg. below determine the mandates that are interacted with the exposure(s) of the crop(s) of interest.
# Not the dynamics wpecified for control crops (this is done by annual_rfs_controls), 
# nor the dynamic effects that are to be aggregated eventually. 
rfs_lead = 2
rfs_lag = 2
rfs_fya = 0 
rfs_pya = 0
# These, below, determine the dynamic effects to aggregate
aggr_lead = 1
aggr_lag = 1

# those are always FALSE 
group_exposure_rfs <- FALSE
most_correlated_only = FALSE
exposure_pasture <- FALSE 
sjpos <- FALSE # should the sample be restricted to cells where sj is positive? 

# control heterogeneity
control_pasture <- FALSE
pasture_trend <- FALSE
fc_trend <- FALSE
exposure_quantiles <- 100
s_trend <- FALSE
s_trend_loga <- FALSE
fc_s_trend <- FALSE
fe = "grid_id + country_year" #   
preclean_level = "FE"
distribution <- "quasipoisson"
offset <- FALSE
invhypsin = FALSE
conley_cutoff <- 100
clustering = "oneway" # either "oneway" or "twoway". If oneway, it clusters on cluster_var1. 
cluster_var1 = "grid_id_10" 
cluster_var2 = "grid_id_5_year"
# dyn_tests = TRUE
rfs_rando <- ""
output = "est_object"
glm_iter <- 25 

# specified arguments
continent = "America" 
outcome_variable = "loss_cropland"
start_year = est_parameters[["start_year"]]
end_year = est_parameters[["end_year"]]
exposure_rfs = "eaear_Cereals"
control_all_absolute_rfs = TRUE
annual_rfs_controls = TRUE
all_exposures_rfs = crops_ctrl
trade_exposure = est_parameters[["trade_exposure"]]
trade_exposure_period = est_parameters[["trade_exposure_period"]]
trade_expo_spec = est_parameters[["trade_expo_spec"]]
access_exposure = est_parameters[["access_exposure"]]
control_access = est_parameters[["control_access"]]
exposure_quantiles = est_parameters[["exposure_quantiles"]]
rfs_lead = est_parameters[["leads"]] 
rfs_lag = est_parameters[["lags"]]
aggr_lead = est_parameters[["aggr_lead"]] 
aggr_lag = est_parameters[["aggr_lag"]]
rfs_fya = est_parameters[["fya"]] 
rfs_pya = est_parameters[["pya"]] 
s_trend = est_parameters[["s_trend"]]
s_trend_loga = est_parameters[["s_trend_loga"]]
s_trend_sq = est_parameters[["s_trend_sq"]]
output = "coef_table"  


rm(outcome_variable, start_year, end_year, continent, 
   original_rfs_treatments, rfs_lead, rfs_lag, exposure_rfs, group_exposure_rfs, control_absolute_rfs, control_all_absolute_rfs, remaining, sjpos, fe, distribution, invhypsin, conley_cutoff, se, boot_cluster, coefs_to_aggregate, 
   output, glm_iter)

make_main_reg <- function(pre_process = FALSE, 
                          pre_processed_data = NULL,
                          outcome_variable = "loss_commodity", # one of "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
                          start_year = 2011, 
                          end_year = 2019, 
                          continent = "all", # one of "Africa", "America", "Asia", or "all"
                          
                          # Joint / dis-joint estimation of crop-specific effects  
                          exposure_rfs = NULL,
                          control_all_absolute_rfs = NULL, # this is to be set to FALSE when doing joint estimation (i.e. all crops passed to exposure_rfs)
                          annual_rfs_controls = FALSE,
                          all_exposures_rfs = NULL, # necessary in disjoint estimation, to govern which other crops are controlled for. 
                          
                          # More exposures
                          trade_exposure = NULL, # either NULL, FALSE, or whatever, for no interaction with a trade exposure, or either "trade_expo", "export_expo", for a (X+I)/Y or X/Y additional crop-specific, country-level, time invariant exposure layer (can be "trade_expo_imp", "export_expo_imp" too)  
                          trade_exposure_period = "20012007", # either "20012007" or "20062007". Used only if trade_exposure is activated with the previous argument. 
                          trade_expo_spec = c(""), # character vectors with terms some or all of terms: "incl_maize", "incl_j", "incl_not_j_not_maize"
                          access_exposure = FALSE,

                          # Treatment dynamics                     
                          original_rfs_treatments = c("statute_conv"),
                          # These 4 arg. below determine the mandates that are interacted with the exposure(s) of the crop(s) of interest.
                          # Not the dynamics wpecified for control crops (this is done by annual_rfs_controls), 
                          # nor the dynamic effects that are to be aggregated eventually. 
                          rfs_lead = 0,
                          rfs_lag = 0,
                          rfs_fya = 0, 
                          rfs_pya = 0,
                          # These, below, determine the dynamic effects to aggregate
                          aggr_lead = rfs_lead,
                          aggr_lag = rfs_lag,

                          # old miscellaneous, not used anymore
                          group_exposure_rfs = FALSE, # This is deprecated, as it was relevant only when rfs exposures where standardized and could thus be added. If exposure_rfs is length 1, it does not matter whether this option is TRUE or FALSE.
                          most_correlated_only = FALSE, # but this restricts the controls to only interactions with the most correlated crop. 
                          sjpos = FALSE, # should the sample be restricted to cells where sj is positive? 
                          
                          # heterogeneity control
                          control_access = FALSE,
                          control_pasture = FALSE,
                          pasture_trend = FALSE,
                          remaining = FALSE, # should remaining forest be controlled for STOP DOING THIS BECAUSE IT INTRODUCES NICKELL BIAS
                          exposure_quantiles = 0, # if higher than 0, trends are not on exact levels of exposure, but on as many quantiles
                          s_trend = "NO", # either "vary_slope", "manual", or any other string (like "NO")
                          s_trend_loga = "NO", # either "vary_slope", "manual", or any other string (like "NO")
                          s_trend_sq = "NO", # either "vary_slope", "manual", or any other string (like "NO")
                          fc_trend = FALSE,
                          fc_s_trend = FALSE,
                          fe = "grid_id + country_year", 
                          preclean_level = "FE", 
                          
                          # estimation - distributional and inferential assumptions
                          distribution = "quasipoisson",#  "quasipoisson", 
                          invhypsin = TRUE, # if distribution is gaussian, should the dep. var. be transformed to inverse hyperbolic sine?
                          
                          clustering = "oneway", # either "oneway" or "twoway". If oneway, it clusters on cluster_var1. 
                          cluster_var1 = "grid_id_10", 
                          cluster_var2 = "grid_id_10",
                          # se = "twoway",# # passed to vcov argument. Currently, one of "cluster", "twoway", or an object of the form: 
                          # - "exposure2ways" for the two-way clustering to be customed as exposure_rfs + cluster_var2, and not along FE.   
                          # - vcov_conley(lat = "lat", lon = "lon", cutoff = 100, distance = "spherical")
                          # with cutoff the distance, in km, passed to fixest::vcov_conley, if se = "conley"  
                          # boot_cluster ="grid_id",
                          # old argument: cluster ="grid_id", # the cluster level if se = "cluster" (i.e. one way)
                          # coefstat = "confint", # one of "se", "tstat", "confint"
                          rfs_rando = "", # either "between", "within", or any other string. If one of the former two, randomization inference of the type is performed
                          glm_iter = 25,
                          # dyn_tests = FALSE, # should the Fisher-type panel unit root test be returned, instead of the regressions, for the outcome_variable and the first regressor
                          output = "coef_table" # one of "data", est_object, or "coef_table" 
){
  
  # Define the outcome_variable based on the crop under study (if we are not in the placebo case)
  # if(outcome_variable == "tmf_agri" & ("eaear_Oilpalm" %in% exposure_rfs | "eaear_Rubber" %in% exposure_rfs)){
  #   outcome_variable <- "tmf_plantation"
  # }
  
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
  
  # this options is to gain some time, when the data is always the same across crop exposure/regressions
  if(pre_process){
    
    d <- pre_processed_data
    
  } else {
    # manipulate a different data set so that original one can be provided to all functions and not read again every time. 
    d <- main_data
    
    # this is just useful if we introduce non-linear trends
    d <- dplyr::mutate(d, year_ln = log(year - min(year) + 1))
    d <- dplyr::mutate(d, year_sq = (year - min(year))^2)
    
    if(grepl("tmf_", outcome_variable)){
      d  <- readRDS(here("temp_data", "merged_datasets", "tmf_aoi", "tmf_aeay_pantrop_long_final_1990_2020.Rdata"))
      # release some memory upfront
      d <- dplyr::filter(d, year >= 2008, year <= 2019)
      
      d <- dplyr::mutate(d, tmf_deforestation = tmf_agri + tmf_plantation)
    }
    
    # ACCESS EXPOSURE
    if(access_exposure){
      d <- dplyr::mutate(d, hours_50kcity_inv = 1/hours_50kcity) # creates Inf values that will be removed below, it's former NAs in the sea
            
      # weight all eaear exposures by this 
      d <- d %>% 
        mutate(across(.cols = starts_with("eaear_"), 
                      .fns = ~.*hours_50kcity_inv))
    }
    
    # TRADE INTERACTION
    all_eaear_trade_exposure_rfs <- c()
    if(length(trade_exposure)>0){
      trade_expo_dat <- readRDS(here("temp_data", "processed_trade_exposures", paste0("trade_exposures_", trade_exposure_period,".Rdata")))
      
      maize_expo_dat <- trade_expo_dat[, (grepl(x = names(trade_expo_dat), pattern = "maizegrain") | names(trade_expo_dat)=="country_name") ]
      
      # if trade_exposure is not _imp, the following line will keep both expo and expo_imp variables, but it already divides by 2 the number of cols that are joined to the bigger d dataframe
      trade_expo_dat <- trade_expo_dat[, ((grepl(x = names(trade_expo_dat), pattern = trade_exposure) | 
                                          names(trade_expo_dat)=="country_name")) &
                                         !names(trade_expo_dat) %in% names(maize_expo_dat)[names(maize_expo_dat)!="country_name"]]
      
      
      # make the average of trade exposures across all crops but the focal one - ALL THIS IS A FUNCTION OF exposure_rfs
      nm_dat <- names(trade_expo_dat)
      all_but_j_and_maize <- nm_dat[!grepl(pattern = "_imp", nm_dat) & 
                                    !grepl("country_name", nm_dat) & 
                                    !grepl(str_to_lower(gsub("eaear_", "", exposure_rfs)), nm_dat) & 
                                    !grepl("maize", nm_dat)]
      
      # all_but_maize <- nm_dat[!grepl(pattern = "_imp", nm_dat) & 
      #                                 !grepl("country_name", nm_dat) & 
      #                                 !grepl("maize", nm_dat)]
      
      # name of that variable 
      trade_expo_not_j_and_maize <- paste0(trade_exposure,"_avg_all_but_j_and_maize")
      # average 
      trade_expo_dat <- 
        trade_expo_dat %>% 
        rowwise() %>% 
        mutate(
          !!as.symbol(trade_expo_not_j_and_maize) := mean(c_across(all_of(all_but_j_and_maize)), na.rm = TRUE)) %>%
        ungroup()
      
      d <- left_join(d, trade_expo_dat, by = "country_name")
      d <- left_join(d, maize_expo_dat, by = "country_name")
      
      if(any(grepl("[.]",names(d)))){stop("merger with export data not as expected")}
      
      # make ALL the eaear-trade exposures
      names(d)
      for(exp_rfs in unique(c(exposure_rfs, all_exposures_rfs))){
        # identify the corresponding trade_expo 
        trade_expo_j <- names(d)[names(d) == paste0(trade_exposure,"_",str_to_lower(gsub(pattern = "eaear_", replacement = "", x = exp_rfs) ) ) ]
        eaear_trade_expo_j <- paste0("eaear_",trade_expo_j)
        d <- mutate(d, !!as.symbol(eaear_trade_expo_j) := !!as.symbol(exp_rfs) * !!as.symbol(trade_expo_j))
        all_eaear_trade_exposure_rfs <- c(all_eaear_trade_exposure_rfs, eaear_trade_expo_j)
      }
      # produce exposure to trade of maize 
      trade_expo_maize <- names(d)[names(d) == paste0(trade_exposure,"_maizegrain")]
      eaear_trade_expo_maize <- paste0("eaear_",trade_expo_maize)
      d <- mutate(d, !!as.symbol(eaear_trade_expo_maize) := eaear_Maizegrain * !!as.symbol(trade_expo_maize))
      all_eaear_trade_exposure_rfs <- unique(c(all_eaear_trade_exposure_rfs, eaear_trade_expo_maize))
      
      # exposure to trade of other crops than j and maize
      # build the variable only for j (exposure_rfs), interaction with all is not a thing anymore
      eaear_trade_expo_not_j_and_maize <- paste0("eaear_",trade_expo_not_j_and_maize)
      d <- mutate(d, !!as.symbol(eaear_trade_expo_not_j_and_maize) := !!as.symbol(exposure_rfs) * !!as.symbol(trade_expo_not_j_and_maize))
      all_eaear_trade_exposure_rfs <- c(all_eaear_trade_exposure_rfs, eaear_trade_expo_not_j_and_maize)
      
      # exposure to import and export of maize (not done above in oil palm process, and not done for both export and import separately)
      d <- mutate(d, eaear_export_expo_maizegrain = eaear_Maizegrain * export_expo_maizegrain)
      d <- mutate(d, eaear_import_expo_maizegrain = eaear_Maizegrain * import_expo_maizegrain)
      
      
      all_eaear_trade_exposure_rfs <- unique(c(all_eaear_trade_exposure_rfs,
                                              "eaear_export_expo_maizegrain", 
                                              "eaear_import_expo_maizegrain"))
      
      rm(trade_expo_dat)
    }
    
    ## ACCESSIBILITY INTERACTION
    # if(access_interaction){
    #   d <- dplyr::mutate(d, hours_50kcity_inv = 1/hours_50kcity) # creates Inf values that will be removed below, it's former NAs in the sea
    #   all_eaear_trade_exposure_rfs <- c(all_eaear_trade_exposure_rfs, "hours_50kcity_inv")
    #   for(exp_rfs in unique(c(exposure_rfs, all_exposures_rfs))){
    #     eaear_acc_expo_j <- paste0("eaear_acc_",exp_rfs)
    #     d <- mutate(d, !!as.symbol(eaear_acc_expo_j) := !!as.symbol(exp_rfs) * hours_50kcity_inv) 
    #     all_eaear_trade_exposure_rfs <- c(all_eaear_trade_exposure_rfs, eaear_acc_expo_j)
    #   }
    # }
    
    
    # code below should work whether d is from tmf or losscommo 
    # Keep only in data the useful variables 
    d <- dplyr::select(d, all_of(unique(c("grid_id", "year", "year_ln", "year_sq", "lat", "lon", "continent_name", "country_name", "country_year",
                                          # "fc_2000", "fc_2009", "remaining_fc", # "accu_defo_since2k",
                                          "grid_id_5", "grid_id_10", "grid_id_20", "grid_id_5_year", "grid_id_10_year", "grid_id_20_year",
                                          outcome_variable,# "tmf_agri", "tmf_flood", "tmf_plantation",
                                          "pasture_share_2000",
                                          "hours_50kcity",
                                          exposure_rfs, all_eaear_trade_exposure_rfs, all_exposures_rfs )))) #sj, 
    
    # Merge only the "prices" needed, not the whole dataframe
    d <- left_join(d, prices[,c("year", unique(c(rfs_treatments, all_rfs_treatments)))], by = c("year"))#, all_treatments
  } 
  
  if((distribution == "gaussian") & invhypsin){
    # transform dependent variable, if gaussian GLM 
    d <- dplyr::mutate(d, !!as.symbol(outcome_variable) := asinh(!!as.symbol(outcome_variable)))
  }
  
  
  ### SPECIFICATION #### 
  
  ### Regressors 
  # if(any(grepl("Fodder", exposure_rfs)) & exposure_pasture){
  #   exposure_rfs[grepl("Fodder", exposure_rfs)] <- "pasture_share_2000"
  # }
  
  # Define which are the eaear_trade exposures in the present regression
  eaear_trade_exposure_rfs <- c()
  if(length(trade_exposure)>0){
    # the one of the same crop
    for(exp_rfs in exposure_rfs){
      # identify the corresponding trade_expo 
      trade_expo_j <- paste0("eaear_",trade_exposure,"_",str_to_lower(gsub(pattern = "eaear_", replacement = "", x = exp_rfs) ) ) 
      if("incl_j" %in% trade_expo_spec){
        eaear_trade_exposure_rfs <- c(eaear_trade_exposure_rfs, trade_expo_j)
      }
    }
    # the one that is averaged over all but the same crop and maize
    if("incl_not_j_and_maize" %in% trade_expo_spec){
      eaear_trade_exposure_rfs <- c(eaear_trade_exposure_rfs, eaear_trade_expo_not_j_and_maize)
    }
    if("incl_maize" %in% trade_expo_spec){
      eaear_trade_exposure_rfs <- c(eaear_trade_exposure_rfs, "eaear_export_expo_maizegrain", "eaear_import_expo_maizegrain")
      
    }
    # # maize, if focal crop isn't already maize
    # if(exposure_rfs!="eaear_Maizegrain"){
    #   eaear_trade_exposure_rfs <- c(eaear_trade_exposure_rfs, paste0("eaear_", trade_exposure, "_maizegrain"))
    # }
  }
  
  eaear_acc_exposure_rfs <- c()
  acc_exposure_rfs <- c() # have this one alone, to interact with mandates but not with AEAY
  # if(access_exposure){
  #   eaear_acc_exposure_rfs <- c(eaear_acc_exposure_rfs, paste0("eaear_acc_",exposure_rfs))
  #   acc_exposure_rfs <- "hours_50kcity_inv"
  # }
  
  
  regressors <- c()
  for(rfs_var in rfs_treatments){
    if(group_exposure_rfs){
      # group (sum) exposures
      d <- dplyr::mutate(d, grp_exp = rowSums(across(.cols = (any_of(exposure_rfs)))))
      # make the regressor of interest
      varname <- paste0("grp_exp_X_", rfs_var)
      regressors <- c(regressors, varname)
      d <- mutate(d, !!as.symbol(varname) := grp_exp * !!as.symbol(rfs_var))
    }else{
      for(exp_rfs in unique(c(exposure_rfs, eaear_trade_exposure_rfs, eaear_acc_exposure_rfs, acc_exposure_rfs))){
        # make regressors of interest
        varname <- paste0(exp_rfs, "_X_", rfs_var)
        regressors <- c(regressors, varname)
        d <- mutate(d, !!as.symbol(varname) := !!as.symbol(exp_rfs) * !!as.symbol(rfs_var))
      }
    }
  }
  
  ### Controls
  # it's important that this is not conditioned on anything so these objects exist
  controls <- c()
  
  # IT IS NORMAL THAT ABS. EXPOSURE CONTROLS ARE INTERACTED WITH 2fya AND 1pya: BOTH ARE AVERAGES OF TWO YEARS. 
  # 1pya is the avg of current and lag1 values, which are grouped together because they reflect the same kind of mechanism. 
  if(control_all_absolute_rfs){
    all_abs <- all_exposures_rfs[!(all_exposures_rfs %in% exposure_rfs)]
    
    # # if we are regressing tmf_plantation, then there is no need to control for all the crops that are not potential drivers (by construction of the product) 
    # if(outcome_variable == "tmf_plantation"){
    #   all_abs <- all_abs[grepl(pattern = "Rubber", all_abs) | grepl(pattern = "Oilpalm", all_abs)]
    #   
    #   # in this case, we control for annual effects on the other one. 
    #   annual_rfs_controls <- TRUE
    # }
    
    # add the controls for the 
    # all_abs <- c(all_abs, eaear_trade_exposure_rfs)
    
    if(most_correlated_only){
      all_abs <- all_abs[all_abs %in% corr_mapmat[corr_mapmat[,"Crops"]==exposure_rfs,"fst_corr"]]
    }
    
    # add the share of pasture in grid cell in 2000 as an exposure to pastures (if the exposure of interest is not pasture)
    if(control_pasture & !grepl("Fodder", exposure_rfs)){
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
  
  if(control_access){
    for(rfs_var in rfs_treatments){
      varname <- paste0("hours_50kcity_X_", rfs_var)
      controls <- c(controls, varname)
      d <- mutate(d, !!as.symbol(varname) := hours_50kcity * !!as.symbol(rfs_var))
    }
  }
  
  # from here, if we are in rfs process, sj may refer to the group of exposures
  if(group_exposure_rfs){
    sj <- "grp_exp"
  }
  

  ### TRENDS ###
  # VARYING SLOPES makes computations much longer (and may not converge)
  # We specify heterogeneous trends using fixest. 
  # as inspired by the answer from Laurent Bergé here https://stackoverflow.com/questions/34232834/fixed-effects-regression-with-state-specific-trends
  # and by the vignettes https://cran.r-project.org/web/packages/fixest/vignettes/fixest_walkthrough.html#41_Interactions_involving_fixed-effects
  # and documentation of fixest https://lrberge.github.io/fixest/reference/feglm.html#varying-slopes-1
  # the point is: we want to allow the year variable to have a different slope (coefficient) for every level of exposure variable
  
  categorized_exposures <- c()
  # eaear_exp_rfs <- c(exposure_rfs, eaear_trade_exposure_rfs)[1]
  for(eaear_exp_rfs in unique(c(exposure_rfs, eaear_trade_exposure_rfs, eaear_acc_exposure_rfs, acc_exposure_rfs))){
    if(exposure_quantiles>0){
      d <- mutate(d, !!as.symbol(paste0(eaear_exp_rfs,"_",exposure_quantiles,"tiles")) := cut(!!as.symbol(eaear_exp_rfs), breaks = exposure_quantiles, labels = FALSE)) # paste0(eaear_exp_rfs,"_Q", 1:exposure_quantiles)
      eaear_exp_rfs <- paste0(eaear_exp_rfs,"_",exposure_quantiles,"tiles")
      categorized_exposures <- c(categorized_exposures, eaear_exp_rfs)
    }
    # d[,c("grid_id", "year", exposure_rfs, eaear_exp_rfs)]
    if(s_trend == "vary_slope"){
      fe <- paste0(fe, paste0(" + ",eaear_exp_rfs,"[[year]]"))
    }  
    if(s_trend_loga == "vary_slope"){
      fe <- paste0(fe, paste0(" + ",eaear_exp_rfs,"[[year_ln]]"))
    }
    if(s_trend_sq == "vary_slope"){
      fe <- paste0(fe, paste0(" + ",eaear_exp_rfs,"[[year_sq]]"))
    }
    
    if(s_trend=="manual"){
      varname <- paste0(eaear_exp_rfs, "_trend")
      controls <- c(controls, varname)
      d <- mutate(d, !!as.symbol(varname) := !!as.symbol(eaear_exp_rfs) * year) # for the linear trend it does not matter
    }
    if(s_trend_loga=="manual"){
      varname <- paste0(eaear_exp_rfs, "_trend_loga")
      controls <- c(controls, varname)
      #make the logarithmic trend start when RFS actually starts, such that it fits best
      # so it should be log(1) in 2005 (although never occurring in the data used for analysis but the log trend starts from there)
      d <- mutate(d, !!as.symbol(varname) := !!as.symbol(eaear_exp_rfs) * year_ln)
    }
    if(s_trend_sq=="manual"){
      varname <- paste0(eaear_exp_rfs, "_trend_sq")
      controls <- c(controls, varname)
      d <- mutate(d, !!as.symbol(varname) := !!as.symbol(eaear_exp_rfs) * year_sq)
    }
  }

  
  # other miscellaneous, deprecated trends 
  # add pasture share trend
  if(pasture_trend){
    d <- mutate(d, pasture_share_trend = pasture_share_2000 * year)
    controls <- c(controls, "pasture_share_trend")
  }
  
  if(fc_trend){
    # # the conditions are just for the sake of computation efficiency
    # if(start_year != 2001 & start_year != 2010){
    #   fc_year <- d[d$year == start_year - 1, c("grid_id", "remaining_fc")]
    #   names(fc_year) <- c("grid_id", paste0("fc_",start_year-1))
    #   d <- left_join(d, fc_year, by = "grid_id")
    #   d <- mutate(d, forest_cover_trend := !!as.symbol(paste0("fc_",start_year-1))  * year)
    # }
    # if(start_year == 2001){
    #   d <- mutate(d, forest_cover_trend = fc_2000 * year)
    # }
    # if(start_year == 2010){
    #   d <- mutate(d, forest_cover_trend = fc_2009 * year)
    # }
    
    fc_year <- d[d$year == start_year - 1, c("grid_id", "tmf_ext")]
    names(fc_year) <- c("grid_id", paste0("fc_",start_year-1))
    d <- left_join(d, fc_year, by = "grid_id")
    d <- mutate(d, forest_cover_trend := !!as.symbol(paste0("fc_",start_year-1))  * year)
    
    controls <- c(controls, "forest_cover_trend")
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
  
  
  #### FORMULA ####  
  if(length(controls) > 0){ 
    alpha_model <- as.formula(paste0(outcome_variable,
                                     " ~ ",
                                     paste0(regressors, collapse = "+"),
                                     " + ",
                                     paste0(controls, collapse = "+"),
                                     " | ",
                                     fe))
  }else{
    alpha_model <- as.formula(paste0(outcome_variable,
                                     " ~ ",
                                     paste0(regressors, collapse = "+"),
                                     " | ",
                                     fe))
  }
  
  ### KEEP OBSERVATIONS THAT: ####
  
  # if pre-processed data is supplied, this does not need to be done
  if(pre_process){
    d_clean <- d
    rm(d)
  } else {
    # - are suitable to crop j 
    if(sjpos){
      d <- dplyr::filter(d, !!as.symbol(exposure_rfs) > 0)  
    }
    
    # - are in study period
    # if(start_year != 2011 | end_year != 2019){
    d <- dplyr::filter(d, year >= start_year)
    d <- dplyr::filter(d, year <= end_year)
    # }
    
    # - are in study area
    if(continent != "all"){
      d <- dplyr::filter(d, continent_name == continent)
    }
    
    # have remaining forest
    # d <- dplyr::filter(d, remaining_fc > 0)
    
    # remove units with no country name (in the case of all_drivers data set currently, because nearest_feature function has not been used in this case, see add_variables.R)
    # d <- dplyr::filter(d, !is.na(country_name))
    
    used_vars <- unique(c("grid_id", "year", "year_ln", "year_sq", "lat", "lon","continent_name",  "country_year",  #"country_name",  "remaining_fc", "accu_defo_since2k", # "sj_year",
                          "grid_id_5", "grid_id_10", "grid_id_20", "grid_id_5_year", "grid_id_10_year", "grid_id_20_year",
                          outcome_variable,# "tmf_agri", "tmf_flood", "tmf_plantation",
                          regressors, controls, 
                          categorized_exposures,
                          exposure_rfs, original_rfs_treatments, rfs_treatments))    # this is necessary to reconstruct variables in randomization inference processes
    
    # DO NOT INCLUDE all_eaear_trade_exposure_rfs IN used_vars BECAUSE IT MAY HAVE MANY MISSINGS, 
    # and we won't actually need all its elements
    
    
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
    rm(temp_est)
    # if we don't remove obs that have no remaining forest, and if we don't remove always 0 outcome obs. in previous step, d IS A BALANCED PANEL !!
    # this is a matter only if we do beta process, with cluster bootstrap
    # if(!nrow(d) == length(unique(d$grid_id))*length(unique(d$year))){
    #   warning("data is not balanced")
    # }
    rm(d)
  }
  
  
  ### Fisher-type panel unit root test  
  # We need to convert the data to a pdata.frame format, to use this interface in purtest 
  # pd_clean <- pdata.frame(d_clean[,c("grid_id", "year", outcome_variable, regressors[1])], index = c("grid_id", "year"))
  # wide_d_clean <- reshape(pd_clean, direction = "wide", v.names = c(outcome_variable, regressors[1]), 
  #                         timevar = "grid_id", idvar = "year", times = "grid_id")
  # purt <- list()
  # for(correction in c("none", "intercept", "trend")){
  #   purt[[correction]] <- plm::purtest(pd_clean[,outcome_variable],
  #                                      formula = as.formula("y ~ trend"),
  #                                      pmax = length(unique(d_clean$year)) - 1, 
  #                                      exo = correction,
  #                                      test = "levinlin")
  # }
  # note that the results of the tests are independent of the crop chosen for exposure. 
  
  ### REGRESSIONS ####
  
  # Store only information necessary, in a dataframe. otherwise the output of fixest estimation is large and we can't collect too many at the same time (over loops)  
  # either there are several elemnts in regressors, and then we want to aggregate them, or there is only one. 
  # In both cases, we are interested in a one-line output
  
  # handle SE computation flexibly within feglm now, through argument vcov
  
  # in case of unit fe, adjust also the vcov argument, because two way wont work, but we might still want to cluster on country year
  if(clustering == "exposure2ways"){
    se <- as.formula(paste0("~ ", exposure_rfs, " + ", cluster_var2))
  }
  
  if(is.numeric(cluster_var1)){
    d_clean <- mutate(d_clean, !!as.symbol(paste0(exposure_rfs[1],"_",cluster_var1,"tiles")) := cut(!!as.symbol(exposure_rfs[1]), breaks = cluster_var1, labels = paste0(exposure_rfs[1],"_Q", 1:cluster_var1)))
    cluster_var1 <- paste0(exposure_rfs[1],"_",cluster_var1,"tiles")
  }
  
  if(clustering =="twoway"){
    se <- as.formula(paste0("~ ", paste0(c(cluster_var1, cluster_var2), collapse = "+")))
  }
  if(clustering =="oneway"){
    se <- as.formula(paste0("~ ", cluster_var1))
  }
  
  
  reg_res <- fixest::feglm(alpha_model,
                           data = d_clean, 
                           family = distribution,# "gaussian",#  # "poisson" ,
                           vcov = se,
                           # this is just to get the same p value by recomputing by hand below. 
                           # see https://cran.r-project.org/web/packages/fixest/vignettes/standard_errors.html
                           # ssc = ssc(cluster.df = "conventional", t.df = "conventional"),
                           # glm.iter = 25,
                           # fixef.iter = 1000,
                           nthreads = 3,
                           fixef.rm = "perfect",
                           glm.iter = glm_iter,
                           notes = TRUE) #, verbose = 3
  
  df_res <- summary(reg_res)$coeftable#[paste0(original_sj, "_X_", original_Pk), ]
  
  fixest_df <- degrees_freedom(reg_res, type = "t")
  
  ## MAKE AGGREGATE RESULTS 
  
  ## Annual mandates
  # joint OR DISJOINT estimation
  if((rfs_lead > 0 | rfs_lag > 0) & rfs_fya == 0 & rfs_pya == 0 ){
    for(EOI in unique(c(exposure_rfs, eaear_trade_exposure_rfs, eaear_acc_exposure_rfs, acc_exposure_rfs))){
      # In this case, we are interested in LEAD AND LAG effects, aggregated separately, and all together.
      df_res <- rbind(rep(NA, ncol(df_res)), rep(NA, ncol(df_res)), df_res)
      
      # ORDER MATTERS
      aggr_names <- paste0(EOI, c("_X_aggrleads", "_X_aggrlags"))
      row.names(df_res)[1:2] <- aggr_names
      # regressors of interest
      base_reg_name <- paste0(EOI,"_X_",original_rfs_treatments)
      # Contemporaneous value is NOT in leads
      lead_roi <- c(grep(pattern = paste0(base_reg_name,"_lead"), 
                         regressors, value = TRUE))
      # Contemporaneous value is in lags
      lag_roi <- c(base_reg_name, 
                   grep(pattern = paste0(base_reg_name,"_lag"), 
                        regressors, value = TRUE))
      
      df_res[aggr_names[1],"Estimate"] <- reg_res$coefficients[lead_roi] %>% sum()
      df_res[aggr_names[2],"Estimate"] <- reg_res$coefficients[lag_roi] %>% sum()
      
      # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
      # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
      df_res[aggr_names[1],"Std. Error"] <- reg_res$cov.scaled[lead_roi, lead_roi] %>% as.matrix() %>% sum() %>% sqrt()
      df_res[aggr_names[2],"Std. Error"] <- reg_res$cov.scaled[lag_roi, lag_roi] %>% as.matrix() %>% sum() %>% sqrt()
      
      df_res[aggr_names[1],"t value"]  <- (df_res[aggr_names[1],"Estimate"] - 0)/(df_res[aggr_names[1],"Std. Error"])
      df_res[aggr_names[2],"t value"]  <- (df_res[aggr_names[2],"Estimate"] - 0)/(df_res[aggr_names[2],"Std. Error"])
      
      # use t distribution with degrees of freedom equal to that used by fixest, i.e. after two way cluster adjustment.
      # does not make a significant difference given sample size
      df_res[aggr_names[1],"Pr(>|t|)"]  <- (2*pt(abs(df_res[aggr_names[1],"t value"]), 
                                                 lower.tail = FALSE, 
                                                 df = fixest_df)) 
      df_res[aggr_names[2],"Pr(>|t|)"]  <- (2*pt(abs(df_res[aggr_names[2],"t value"]), 
                                                 lower.tail = FALSE, 
                                                 df = fixest_df)) 
      
      
      ## Then (on top of dataframe), aggregate all leads and lags together
      # it's important that overall aggregate comes after, for at least two reasons in current code:
      # 1. because selection of estimate to plot is based on order: it takes the first row of df_res 
      # 2. so that aggr_names corresponds to *_X_aggrall in randomization inference below
      
      df_res <- rbind(rep(NA, ncol(df_res)), df_res)
      
      # ORDER MATTERS
      aggr_names <- paste0(EOI, c("_X_aggrall"))
      row.names(df_res)[1] <- aggr_names
      # regressors of interest
      base_reg_name <- paste0(EOI,"_X_",original_rfs_treatments)
      
      # Contemporaneous, lead, and lag values, up to those we want to account for 
      all_roi <- c(base_reg_name, 
                   paste0(base_reg_name,"_lead",0:aggr_lead),
                   paste0(base_reg_name,"_lag",0:aggr_lag))
      
      # Check that they were in the estimation 
      all_roi <- all_roi[sapply(all_roi, function(roi){(roi %in% regressors)})]
      
      # Keep only 
      df_res[aggr_names[1],"Estimate"] <- reg_res$coefficients[all_roi] %>% sum()
      
      # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
      # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
      df_res[aggr_names[1],"Std. Error"] <- reg_res$cov.scaled[all_roi, all_roi] %>% as.matrix() %>% sum() %>% sqrt()
      
      df_res[aggr_names[1],"t value"]  <- (df_res[aggr_names[1],"Estimate"] - 0)/(df_res[aggr_names[1],"Std. Error"])
      
      # use t distribution with degrees of freedom equal to that used by fixest, i.e. after two way cluster adjustment.
      # does not make a significant difference given sample size
      df_res[aggr_names[1],"Pr(>|t|)"]  <- (2*pt(abs(df_res[aggr_names[1],"t value"]), 
                                                 lower.tail = FALSE, 
                                                 df = fixest_df)) 
    }
  }
  
  
  ## Averaged mandates
  if(control_all_absolute_rfs & rfs_lead == 0 & rfs_lag == 0 & rfs_fya > 0 & rfs_pya > 0){
    # In this case, we are interested in AVERAGE LEAD AND LAG effects, separately and aggregated
    df_res <- rbind(rep(NA, ncol(df_res)), df_res)
    
    # ORDER MATTERS
    aggr_names <- paste0(exposure_rfs, c("_X_aggrall"))
    row.names(df_res)[1] <- aggr_names
    # regressors of interest
    base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
    
    # Contemporaneous, lead, and lag values
    all_roi <- c(grep(pattern = paste0(base_reg_name,"_",rfs_fya,"fya"), 
                      regressors, value = TRUE),
                 grep(pattern = paste0(base_reg_name,"_",rfs_pya,"pya"), 
                      regressors, value = TRUE))
    
    df_res[aggr_names[1],"Estimate"] <- reg_res$coefficients[all_roi] %>% sum()
    
    # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
    # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
    df_res[aggr_names[1],"Std. Error"] <- reg_res$cov.scaled[all_roi, all_roi] %>% as.matrix() %>% sum() %>% sqrt()
    
    df_res[aggr_names[1],"t value"]  <- (df_res[aggr_names[1],"Estimate"] - 0)/(df_res[aggr_names[1],"Std. Error"])
    
    # use t distribution with degrees of freedom equal to that used by fixest, i.e. after two way cluster adjustment.
    # does not make a significant difference given sample size
    df_res[aggr_names[1],"Pr(>|t|)"]  <- (2*pt(abs(df_res[aggr_names[1],"t value"]), 
                                               lower.tail = FALSE, 
                                               df = fixest_df)) 
  }
  if(control_all_absolute_rfs & rfs_lead > 0 & rfs_lag == 0 & rfs_fya == 0 & rfs_pya > 0){
    # In this case, we are interested in aggregated LEAD effects, and all together.
    
    # first aggregate leads
    df_res <- rbind(rep(NA, ncol(df_res)), df_res)
    
    # ORDER MATTERS
    aggr_names <- paste0(exposure_rfs, c("_X_aggrleads"))
    row.names(df_res)[1] <- aggr_names
    # regressors of interest
    base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
    # Contemporaneous value is NOT in leads
    lead_roi <- c(grep(pattern = paste0(base_reg_name,"_lead"), 
                       regressors, value = TRUE))
    
    df_res[aggr_names[1],"Estimate"] <- reg_res$coefficients[lead_roi] %>% sum()
    
    # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
    # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
    df_res[aggr_names[1],"Std. Error"] <- reg_res$cov.scaled[lead_roi, lead_roi] %>% as.matrix() %>% sum() %>% sqrt()
    
    df_res[aggr_names[1],"t value"]  <- (df_res[aggr_names[1],"Estimate"] - 0)/(df_res[aggr_names[1],"Std. Error"])
    
    # use t distribution with degrees of freedom equal to that used by fixest, i.e. after two way cluster adjustment.
    # does not make a significant difference given sample size
    df_res[aggr_names[1],"Pr(>|t|)"]  <- (2*pt(abs(df_res[aggr_names[1],"t value"]), 
                                               lower.tail = FALSE, 
                                               df = fixest_df)) 
    
    ## Then (on top of dataframe), aggregate all leads and lags together
    # it's important that overall aggregate comes after, for at least two reasons in current code:
    # 1. because selection of estimate to plot is based on order: it takes the first row of df_res 
    # 2. so that aggr_names corresponds to *_X_aggrall in randomization inference below
    
    df_res <- rbind(rep(NA, ncol(df_res)), df_res)
    
    # ORDER MATTERS
    aggr_names <- paste0(exposure_rfs, c("_X_aggrall"))
    row.names(df_res)[1] <- aggr_names
    # regressors of interest
    base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
    
    # Contemporaneous, lead, and lag values
    all_roi <- c(grep(pattern = paste0(base_reg_name,"_lead"), 
                      regressors, value = TRUE),
                 grep(pattern = paste0(base_reg_name,"_",rfs_pya,"pya"), 
                      regressors, value = TRUE))
    
    df_res[aggr_names[1],"Estimate"] <- reg_res$coefficients[all_roi] %>% sum()
    
    # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
    # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
    df_res[aggr_names[1],"Std. Error"] <- reg_res$cov.scaled[all_roi, all_roi] %>% as.matrix() %>% sum() %>% sqrt()
    
    df_res[aggr_names[1],"t value"]  <- (df_res[aggr_names[1],"Estimate"] - 0)/(df_res[aggr_names[1],"Std. Error"])
    
    # use t distribution with degrees of freedom equal to that used by fixest, i.e. after two way cluster adjustment.
    # does not make a significant difference given sample size
    df_res[aggr_names[1],"Pr(>|t|)"]  <- (2*pt(abs(df_res[aggr_names[1],"t value"]), 
                                               lower.tail = FALSE, 
                                               df = fixest_df)) 
  }
  
  # eyeball results
  # df_res[sort(row.names(df_res)), ]
  # df_res[grepl("aggrall", row.names(df_res)),]
  
  # take data set as exactly used in estimation - NOT NECESSARY anymore, given the precleaning
  # if(length(reg_res$obs_selection) > 0){
  #   d_clean <- d_clean[reg_res$obs_selection[[1]], ]
  # }

  
  # store info on SE method
  if(clustering == "conley"){
    se_info <- paste0("Conley (",conley_cutoff,"km)")
  }
  if(clustering == "twoway"){
    # " If the two variables were used as fixed-effects in the estimation, 
    # you can leave it blank with vcov = "twoway""
    se_info <- "Two-way clustered (grid cell - year)"
  }
  if(clustering == "oneway"){
    se_info <- paste0("Clustered (",cluster_var1,")")
  }
  
  #df_res$inference <- se_info  
  
  
  # output wanted
  # We can't estimate VCOV if we reduce memory size of fixest object with lean
  # reg_res <- summary(reg_res, lean = TRUE)
  
  if(output == "everything"){
    toreturn <- list(reg_res, d_clean, df_res) 
  }
  if(output == "est_object"){ 
    toreturn <- reg_res
  }
  if(output == "data"){
    toreturn <- list(df_res, d_clean)
  }
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
  
  
  rm(d_clean, df_res, reg_res)
  return(toreturn)
  rm(toreturn)
}

# POST REG FUNCTION ----------------------
# Post estimation function to aggregate annual effects and compute t stats 
# this is data less for now: does not recompute SE 

# For debugging
est_obj <- all_tests_est[[CNT]][["eaear_Maizegrain"]]
d_clean <- d_clean_list[[CNT]][["loss_cropland"]]
base_exposure = "eaear_Maizegrain"
trade_exposure = est_parameters[["trade_exposure"]]
# access_exposure = est_parameters[["access_exposure"]] 
aggr_lag = est_parameters[["aggr_lag"]]
aggr_lead = est_parameters[["aggr_lead"]]
clustering = est_parameters[["clustering"]]
cluster_var1 = est_parameters[["cluster_var1"]]
cluster_var2 = est_parameters[["cluster_var2"]]
output = "tstat"

rm(est_obj, d_clean, base_exposure, aggr_lag, aggr_lead, clustering, cluster_var1, cluster_var2, output)

post_est_fnc <- function(est_obj, # a fixest estmation object
                         # d_clean,
                         base_exposure = CROP,
                         trade_exposure = NULL, # default
                         # access_exposure = FALSE, # default
                         aggr_lag = est_parameters[["aggr_lag"]], 
                         aggr_lead = est_parameters[["aggr_lead"]], 
                         clustering = est_parameters[["clustering"]],
                         cluster_var1 = est_parameters[["cluster_var1"]],
                         cluster_var2 = est_parameters[["cluster_var2"]],
                         output = "tstat"){
  
  df_res <- summary(est_obj)$coeftable#[paste0(original_sj, "_X_", original_Pk), ]
  
  regressors <- names(coef(est_obj))
  
  fixest_df <- degrees_freedom(est_obj, type = "t")
  
  df_res <- df_res[grepl(base_exposure, row.names(df_res)), ]
  
  exposures_of_interest <- base_exposure
  if(length(trade_exposure)>0){
    exposures_of_interest <- c(exposures_of_interest, 
                               paste0("eaear_",trade_exposure,"_",str_to_lower(gsub(pattern = "eaear_", replacement = "", x = base_exposure) ) ) )
  } 
  # if(access_exposure){
  #   exposures_of_interest <- c(exposures_of_interest, 
  #                              paste0("eaear_acc_",base_exposure))
  # }
  
  # RE-COMPUTE THE VCOV, ACCORDING TO THE SPECIFIED CLUSTERING 
  # cluster var on percentile is constructed outside this function 
  if(is.numeric(cluster_var1)){
    # d_clean <- mutate(d_clean, !!as.symbol(paste0(CROP[1],"_",est_parameters[["cluster_var1"]],"tiles")) := cut(!!as.symbol(CROP[1]), breaks = est_parameters[["cluster_var1"]], labels = paste0(CROP[1],"_Q", 1:est_parameters[["cluster_var1"]])))
    cluster_var1 <- paste0(base_exposure[1],"_",cluster_var1,"tiles")
  }
  if(clustering =="twoway"){
    se <- as.formula(paste0("~ ", paste0(c(cluster_var1, cluster_var2), collapse = "+")))
  }
  if(clustering =="oneway"){
    se <- as.formula(paste0("~ ", cluster_var1))
  }
  
  est_obj <- summary(est_obj, cluster = se)
  
  for(EOI in exposures_of_interest){
    # regressors of interest
    base_reg_name <- paste0(EOI,"_X_statute_conv")
    
    ## MAKE AGGREGATE RESULTS 
    
    ## On top of dataframe, aggregate all leads and lags together
    
    df_res <- rbind(rep(NA, ncol(df_res)), df_res)
    
    aggr_names <- paste0(EOI, c("_X_aggrall"))
    row.names(df_res)[1] <- aggr_names
    
    
    # Contemporaneous, lead, and lag values, up to those we want to account for 
    all_roi <- c(base_reg_name, 
                 paste0(base_reg_name,"_lead",0:aggr_lead),
                 paste0(base_reg_name,"_lag",0:aggr_lag))
    
    # Check that they were in the estimation (i.e. remove lead0 or lag0 in particular)
    all_roi <- all_roi[sapply(all_roi, function(roi){(roi %in% regressors)})]
    
    # Keep only 
    df_res[aggr_names[1],"Estimate"] <- est_obj$coefficients[all_roi] %>% sum()
    
    # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
    # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
    df_res[aggr_names[1],"Std. Error"] <- est_obj$cov.scaled[all_roi, all_roi] %>% as.matrix() %>% sum() %>% sqrt()
    
    df_res[aggr_names[1],"t value"]  <- (df_res[aggr_names[1],"Estimate"] - 0)/(df_res[aggr_names[1],"Std. Error"])
    
    # use t distribution with degrees of freedom equal to that used by fixest, i.e. after two way cluster adjustment.
    # does not make a significant difference given sample size
    df_res[aggr_names[1],"Pr(>|t|)"]  <- (2*pt(abs(df_res[aggr_names[1],"t value"]), 
                                               lower.tail = FALSE, 
                                               df = fixest_df)) 
  }
  
  
  if(output=="tstat"){
    post_est_return <- data.frame(tstat = df_res[,"t value"]) 
  } else {
    post_est_return <- df_res
  }
  
  return(post_est_return)
}


#### Main specification and objects -----------------------------------------------------------------------
est_parameters <- list(start_year = 2008, 
                       end_year = 2016, 
                       trade_exposure = NULL, # "export_expo_imp0", # "trade_expo", #  # #  "trade_expo_imp", # NULL, # "trade_expo_imp",
                       trade_exposure_period = "20012007",
                       trade_expo_spec = c(""),
                       access_exposure = FALSE,
                       control_access = FALSE,
                       annual_rfs_controls = TRUE,
                       leads = 2,
                       lags = 2,
                       fya = 0, 
                       pya = 0,
                       aggr_lead = 1,
                       aggr_lag = 1, 
                       sjpos= FALSE,
                       exposure_quantiles = 0,
                       s_trend = "NO", # either "vary_slope", "manual" or "NO" / whatever
                       s_trend_loga = "manual",
                       s_trend_sq = "NO",
                       clustering = "oneway",
                       cluster_var1 = "grid_id_10", 
                       cluster_var2 = "grid_id_10", 
                       glm_iter = 50,
                       output = "est_object") # !!! MIND THIS ARG "est_object", "everything" "coef_table"

# This sets the crops that are controled for 
crops_groups <- list(` AMaize` = "eaear_Maizegrain",
                    # marg_land = c("eaear_Biomass"), #"eaear_Fodder",  this only serves as a control
                    `Type 1` = c("eaear_Cereals", "eaear_Rice", "eaear_Roots", "eaear_Fodder", "eaear_Soy_compo"), #
                    `Type 2` = c("eaear_Cotton",  "eaear_Oilfeed_crops"),
                    `Type 4` = c("eaear_Biomass", "eaear_Sugarcane"), #
                    `Type 5` = c("eaear_Tobacco"),
                    `Type 3` = c("eaear_Oilpalm"),
                    `Type 6` = c("eaear_Citrus", "eaear_Coconut", "eaear_Banana", "eaear_Cocoa_Coffee", "eaear_Tea", "eaear_Rubber"))

# crops_groups <- list(Maize = "eaear_Maizegrain",
#                      # marg_land = c("eaear_Biomass"), #"eaear_Fodder",  this only serves as a control 
#                      `Group 1 - 2` = c("eaear_Fodder", "eaear_Soy_compo", "eaear_Cereals", "eaear_Rice", "eaear_Cotton",  "eaear_Oilfeed_crops"), 
#                      `Group 4` = c("eaear_Biomass", "eaear_Sugarcane"), # 
#                      `Group 5` = c("eaear_Tobacco"), 
#                      `Group 3` = c("eaear_Oilpalm"), 
#                      `Group 6` = c("eaear_Citrus", "eaear_Coconut", "eaear_Banana", "eaear_Cocoa_Coffee", "eaear_Tea", "eaear_Rubber"))

all_crops <- unlist(unname(crops_groups)) 
# fodder appears in g1, but actually we decide below in the loop whether to keep it or not in cropland reg.
# biomass is to be moved from g4 to marg_land in order to estimate sugarcane conditional to it too. 


# All crops, by common outcome that can be regressed on them
gentest_crops <- c("eaear_Maizegrain", 
                "eaear_Cereals", "eaear_Rice", "eaear_Roots", "eaear_Soy_compo", "eaear_Fodder",
                "eaear_Cotton",  "eaear_Oilfeed_crops")

displtest_crops <- c("eaear_Biomass", "eaear_Sugarcane", "eaear_Tobacco")#  
if(length(est_parameters[["trade_exposure"]])>0){
  displtest_crops <- displtest_crops[displtest_crops!="eaear_Biomass"]
}
# oil palm, coconut, and citrus are particular cases

g6_crops <- c("eaear_Banana", "eaear_Cocoa_Coffee", "eaear_Tea", "eaear_Rubber")

continents <- c("all", "America", "Africa", "Asia") # "all"#  


### DATA 
main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "loss_cropcommo_opindusnotrans_residu_aeaycompo_long_final.Rdata"))
# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "loss_cropcommo_op_aeaycompo_long_final.Rdata"))
# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "loss_cropcommo_opcommo_aeaycompo_long_final.Rdata"))
# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "loss_crop_opindus_aeaycompo_long_final.Rdata"))
# this is with cropland outside commo, and without pasture
#main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "loss_aeaycompo_long_final.Rdata"))

#main_data <- dplyr::mutate(main_data, loss_commodity = loss_cropland + loss_oilpalm_both + loss_pasture)
# release some memory upfront
main_data <- dplyr::filter(main_data, year >= 2008, year <= 2019)

# this may be used to remove some controls in loss_cropland regressions
woody_perrenials <- c("eaear_Citrus", "eaear_Coconut", "eaear_Oilpalm", 
                      "eaear_Banana", "eaear_Cocoa_Coffee", "eaear_Tea", "eaear_Rubber")

# this may be used to remove some controls in loss residual regressions (called loss_pasture)
cropland_crops <- c(all_crops[!(all_crops %in% c(woody_perrenials))])

# loss_type <- "loss_cropland"
# CNT <- "America"
# CROP <- "eaear_Soy_compo"


#### LOSS COMMO CATEGORIES ------------------------------------------
all_tests_est <- list()
d_clean_list <- list()
for(CNT in continents){
  
  ### CROPLAND ### 
  loss_type <- "loss_cropland"
  
  # GROUP 1-2
  for(CROP in c(gentest_crops, displtest_crops)){
    ## Determine the crops to control for
    #crops_ctrl <- crops_groups[sapply(crops_groups, match, x = CROP, nomatch = FALSE)==0] %>% unlist() %>% unname()  # this removes from the set of controls, the crops that are in the same group as CROP
    crops_ctrl <- all_crops

    # **NOPE** (control transitory LUC)
    # further remove woody perennials, if we regress cropland anyway
    crops_ctrl <- crops_ctrl[sapply(crops_ctrl, function(c){!(c %in% c(woody_perrenials))})]# , "eaear_Fodder"

    # the order of the strata in the list are determined by how we want to plot results
    temp_list_cropland <- make_main_reg(continent = CNT,
                                         outcome_variable = loss_type,

                                         start_year = est_parameters[["start_year"]],
                                         end_year = est_parameters[["end_year"]],

                                         # regress on one crop at a time, and control for other crops
                                         exposure_rfs = CROP, # "eaear_Maizegrain",
                                         control_all_absolute_rfs = TRUE,
                                         annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                         all_exposures_rfs = crops_ctrl,

                                         trade_exposure = est_parameters[["trade_exposure"]],
                                         trade_exposure_period = est_parameters[["trade_exposure_period"]],

                                         rfs_lead = est_parameters[["leads"]],
                                         rfs_lag = est_parameters[["lags"]],
                                         aggr_lead = est_parameters[["aggr_lead"]],
                                         aggr_lag = est_parameters[["aggr_lag"]],
                                         sjpos = est_parameters[["sjpos"]],
                                         exposure_quantiles = est_parameters[["exposure_quantiles"]],
                                         s_trend = est_parameters[["s_trend"]],
                                         s_trend_loga = est_parameters[["s_trend_loga"]],
                                         s_trend_sq = est_parameters[["s_trend_sq"]],

                                         glm_iter = est_parameters[["glm_iter"]],

                                         output = "everything")
    
    all_tests_est[[CNT]][[CROP]] <- temp_list_cropland[[1]] # this is the fixest object
    d_clean_list[[CNT]][[loss_type]] <- temp_list_cropland[[2]] # this is the data (always the same at every CROP iteration)
    rm(temp_list_cropland)
    rm(CROP)
  }
  
  ### OIL PALM (GROUP 3) ###
  oilpalm_loss <- "loss_oilpalm_indus"
  # oilpalm_ctrl <- cropland_crops
  # oilpalm_ctrl <- all_crops # [all_crops != "eaear_Coconut"]
  oilpalm_ctrl <- c()
  #for(oilpalm_loss in c("loss_oilpalm_both", "loss_oilpalm_indus")){
  temp_list_oilpalm <- make_main_reg(continent = CNT,
                                      outcome_variable = oilpalm_loss,
                                      start_year = est_parameters[["start_year"]],
                                      end_year = est_parameters[["end_year"]],
  
                                      # regress on crop at a time, but control for annual interactions with other crops
                                      exposure_rfs = "eaear_Oilpalm",
                                      control_all_absolute_rfs = TRUE,
                                      annual_rfs_controls = TRUE,
                                      all_exposures_rfs = oilpalm_ctrl,
  
                                     trade_exposure = est_parameters[["trade_exposure"]],
                                     trade_exposure_period = est_parameters[["trade_exposure_period"]],
  
                                      rfs_lead = est_parameters[["leads"]],
                                      rfs_lag = est_parameters[["lags"]],
                                      aggr_lead = est_parameters[["aggr_lead"]],
                                      aggr_lag = est_parameters[["aggr_lag"]],
                                      sjpos = est_parameters[["sjpos"]],
                                     exposure_quantiles = est_parameters[["exposure_quantiles"]],
                                     s_trend = est_parameters[["s_trend"]],
                                     s_trend_loga = est_parameters[["s_trend_loga"]],
                                     s_trend_sq = est_parameters[["s_trend_sq"]],
                                      glm_iter = est_parameters[["glm_iter"]],
  
                                      output = "everything")
  
  all_tests_est[[CNT]][["eaear_Oilpalm"]] <- temp_list_oilpalm[[1]] # this is the fixest object
  d_clean_list[[CNT]][[oilpalm_loss]] <- temp_list_oilpalm[[2]] # this is the data 
  rm(temp_list_oilpalm)
  
}

# est_obj <- gentest_crops_est[[loss_type]][[CNT]][[CROP]][[1]]
# est_dat <- gentest_crops_est[[loss_type]][[CNT]][[CROP]][[2]]
# est_df <- gentest_crops_est[[loss_type]][[CNT]][[CROP]][[3]]
# unique(est_dat$year)
for(CNT in continents){
  sapply(all_tests_est[[CNT]], FUN = function(eo){eo$convStatus}) %>% print()
}


# all_tests_est[["all"]][["loss_oilpalm_both"]] %>% post_est_fnc(base_exposure = "eaear_Oilpalm")
all_tests_est[["America"]][["eaear_Biomass"]] %>% summary()
all_tests_est[["America"]][["eaear_Soy_compo"]] %>% summary()
all_tests_est[["Asia"]][["eaear_Oilpalm"]]$convStatus
all_tests_est[["Asia"]][["eaear_Oilpalm"]] %>% summary()

### PLOT COEFFICIENTS ### 
# prepare regression outputs in a tidy data frame readable by dwplot
# this works annual_controls is TRUE and we add (as controls) all crops that we are interested in. 
rm(df, d_clean)
cnt_crop_list <- list()
cnt_list <- list()
for(CNT in continents){# "all", , "Africa", "Asia"
  for(CROP in c(gentest_crops, displtest_crops, "eaear_Oilpalm")){
    if(CROP == "eaear_Oilpalm"){
      loss_type <- oilpalm_loss
    } else { 
      loss_type <- "loss_cropland"  
    }
    d_clean <- d_clean_list[[CNT]][[loss_type]]
    
    # store post-estimation outputs of all clustering levels
    postest_list <- list()
    # THREE LEVELS OF SPATIAL CLUSTERING
    # Main 
    postest_list[["grid_id_10"]] <- post_est_fnc(est_obj = all_tests_est[[CNT]][[CROP]], # loss_type
                                                base_exposure = CROP,
                                                trade_exposure = NULL,
                                                aggr_lag = est_parameters[["aggr_lag"]], 
                                                aggr_lead = est_parameters[["aggr_lead"]], 
                                                clustering = "oneway",
                                                cluster_var1 = "grid_id_10",
                                                output = "tstat")
    # 5 times larger
    postest_list[["grid_id_5"]] <- post_est_fnc(est_obj = all_tests_est[[CNT]][[CROP]], # loss_type
                                                base_exposure = CROP,
                                                trade_exposure = NULL,
                                                aggr_lag = est_parameters[["aggr_lag"]], 
                                                aggr_lead = est_parameters[["aggr_lead"]], 
                                                clustering = "oneway",
                                                cluster_var1 = "grid_id_5",
                                                output = "tstat")
    # 20 times larger
    postest_list[["grid_id_20"]] <- post_est_fnc(est_obj = all_tests_est[[CNT]][[CROP]], # loss_type
                                                base_exposure = CROP,
                                                trade_exposure = NULL,
                                                aggr_lag = est_parameters[["aggr_lag"]], 
                                                aggr_lead = est_parameters[["aggr_lead"]], 
                                                clustering = "oneway",
                                                cluster_var1 = "grid_id_20",
                                                output = "tstat")
    
    ## TWO-WAY  CLUSTERING
    # SET THE QUANTILE LEVEL HERE 
    q_expo <- 10
    # need to do this here and not in post_est_fnc for fixest to find the data to construct the new vcov 
    # make sure that the exposure var is not there already
    d_clean <- left_join(d_clean[,names(d_clean)!=CROP[1]], main_data[,c("grid_id", "year", CROP)], by = c("grid_id", "year"))
    
    d_clean <- mutate(d_clean, !!as.symbol(paste0(CROP[1],"_",q_expo,"tiles")) := cut(!!as.symbol(CROP[1]), breaks = q_expo, labels = paste0(CROP[1],"_Q", 1:q_expo)))
    
    # extract t stats of coef of interest 
    postest_list[["twoway_percentile"]] <- post_est_fnc(est_obj = all_tests_est[[CNT]][[CROP]], # loss_type
                                                        base_exposure = CROP,
                                                        trade_exposure = NULL,
                                                        aggr_lag = est_parameters[["aggr_lag"]], 
                                                        aggr_lead = est_parameters[["aggr_lead"]],
                                                        clustering = "twoway",
                                                        cluster_var1 = q_expo,
                                                        cluster_var2 = "grid_id_10",
                                                        output = "tstat")
    rm(d_clean)
    
    # save row names 
    estimate_names <- row.names(postest_list[[1]])
      
    x <- bind_cols(postest_list, .name_repair = "minimal")
    names(x) <- names(postest_list)
    x <- dplyr::mutate(x, across(.cols = everything(), 
                                 ~abs(.)>1.645, 
                                 .names = paste0("{.col}","_isprecise")))
    x <- dplyr::mutate(x, rob_preci = rowMeans(across(.cols = contains("_isprecise")) ))
    # keep only the main t-stat
    x <- dplyr::select(x, -grid_id_5, -grid_id_20, -twoway_percentile)
    names(x)[names(x)=="grid_id_10"] <- "tstat"
    
    x$model <- CNT
    x$via_trade <- grepl("_expo_", estimate_names)    
    x$via_access <- grepl("_acc_", estimate_names)    
    x$crop <- CROP
    x$Dynamics <- ""
    x$Dynamics[grepl("aggrall", estimate_names)] <- "Cumulative"
    x$Dynamics[grepl("lead1", estimate_names)] <- "t+1"
    x$Dynamics[grepl("lead2", estimate_names)] <- "t+2"
    x$Dynamics[estimate_names==paste0(CROP,"_X_statute_conv")] <- "t"
    x$Dynamics[grepl("lag1", estimate_names)] <- "t-1"
    x$Dynamics[grepl("lag2", estimate_names)] <- "t-2"
    x$Dynamics[grepl("lag3", estimate_names)] <- "t-3"
    
    # order dynamics 
    # estimate_names <- x$Dynamics  
    # x <- x[c("Cumulative", "t+1", "t", "t-1", "t-2", "t-3", "t-4"),] 
    # # remove NA rows generated if some dynamics are absent for this crop 
    # x <- x[!is.na(x$tstat), ]
    row.names(x) <- NULL # otherwise they duplicate while stacking
    cnt_crop_list[[CNT]][[CROP]] <- x
    
  }
  cnt_list[[CNT]] <- cnt_crop_list[[CNT]] %>% bind_rows()
}
df <- cnt_list %>% bind_rows() # "loss_cropland"

df <- df[df$Dynamics=="Cumulative",] ## CHANGE HERE TO "t+2" to display only anticipation effects

df$significant01 <- ""
df[abs(df[,"tstat"]) > 1.645, "significant01"] <- "p-value < .1"
# df$highlight <- ""
# df[df$significant01 == "p-value < .1" & df$rob_preci==1, "highlight"] <- "p-value < .1 & robust" 
df <- mutate(df, highlight = rob_preci)
# df <- mutate(df, highlight = as.discrete(round((rob_preci+0.5)/1.5, 2)))


# attribute crop Group 
# horrible code mais j'ai pas réussi à faire mieux ... 
df_gn <- sapply(names(crops_groups), function(gn){if_else(df$crop %in% crops_groups[[gn]], gn, "")}) 
df$Group <- ""
for(i in 1:nrow(df_gn)){
 row <- df_gn[i,]
 df$Group[i] <- row[row!=""] %>% unname()
}

# make some variable name changes necessary for dotwhisker
# df <- df[df$crop != "eaear_Maizegrain",]
df <- mutate(df, sizes = if_else(model == "all", true = 4, false = 2.5))
names(df)[names(df)=="tstat"] <- "estimate"
# names(df)[names(df)=="model"] <- "continent"
# names(df)[names(df)=="Dynamics"] <- "model"
names(df)[names(df)=="crop"] <- "term"
# df$cumulative <- factor(df$Dynamics=="Cumulative")
head(df)


dwplot(df,
       dot_args = list(aes(color = model, shape = model, size = model, alpha = highlight)) ) %>%  #size = 2.5,  shape = Dynamics, size = cumulative
  relabel_predictors(predictors_dict) + 
  # # critical values  : 1.645   1.960  2.576
  geom_vline(xintercept = 1.96, colour = "grey60", linetype = "dotdash", alpha = 1) +
  geom_vline(xintercept = -1.96, colour = "grey60", linetype = "dotdash", alpha = 1) +
  geom_vline(xintercept = 2.576, colour = "grey60", linetype = "dotted", alpha = 1) +
  geom_vline(xintercept = -2.576, colour = "grey60", linetype = "dotted", alpha = 1) +
    
  facet_wrap(facets = ~Group, scales = "free_y", nrow = 3, ncol = 2, dir = "v") + 
  
  scale_color_brewer(type = "qual",palette="Set1", #  Accent  
                     breaks=c("all", "America", "Africa", "Asia"),
                     labels = c("Pan-tropical", "America (tropical lat.)", "Africa (tropical lat.)", "Asia (tropical lat.)"),
                     name="") +
  
  scale_shape_manual(values = c(18, 17, 15, 16), # c(21, 24, 22, 25),# c(1, 2, 0, 3), # c(16, 17, 15, 18)
                     breaks=c("all", "America", "Africa", "Asia"),
                     labels = c("Pan-tropical", "America (tropical lat.)", "Africa (tropical lat.)", "Asia (tropical lat.)"),
                     name="") +
  
  scale_size_manual(values = c(3.5, 2.5, 2.5, 2.5),
                    breaks=c("all", "America", "Africa", "Asia"),
                    labels = c("Pan-tropical", "America (tropical lat.)", "Africa (tropical lat.)", "Asia (tropical lat.)"),
                    name="") +
  
  scale_alpha_continuous(# values = sort(unique(df$highlight)),# c(0.33, 0.50, 0.67, 0.83, 1.00), # c(0.2, 0.4, 0.6, 0.8, 1),
                    breaks= c(0.33, 0.50, 0.67, 0.83, 1.00),# sort(unique(df$highlight)),# 
                    labels = c("0", "1/4", "2/4", "3/4", "4/4"),
                    name="Precision robustness") +

  scale_x_continuous(breaks = c(-1.96, 1.96), 
                     labels = c("-1.96","1.96")) +
  
  theme_bw() + xlab("t-statistics") + ylab("") +  
  theme(plot.title = element_text(face="bold", size=c(10)),
        legend.position = "bottom",# c(0.8, 0.05),
        legend.justification = c(0, 0), 
        legend.background = element_rect(colour="grey80")
        ) # legend.title = element_blank()






#### LOSS COMMO CATEGORIES - TRADE -------------------------------------
# In this case, all crop estimations have different samples d_clean because they require their 
# own trade exposure to not be NA, and this varies by country. 
# So don't make post estimation, because this requires to store the data used for estimation, and would be too memory intensive in this case. 
est_parameters[["trade_exposure"]] <- "trade_expo"
est_parameters[["trade_expo_spec"]] <- c("incl_maize")
temp_list_trade <- list()
for(CNT in continents[continents!="all"]){ # 
  
  ### CROPLAND ### 
  loss_type <- "loss_cropland"
  # keep ony those crops for which there is trade data 
  commercial_crops <- c(gentest_crops, displtest_crops)
  commercial_crops <- commercial_crops[!(commercial_crops %in% c("eaear_Biomass", "eaear_Roots"))] # 
  
  # this does not work
  if(CNT=="all"){commercial_crops <- commercial_crops[commercial_crops!="eaear_Oilfeed_crops"]}
  
  # GROUP 1-2
  for(CROP in commercial_crops){# c(gentest_crops, displtest_crops)
    ## Determine the crops to control for
    #crops_ctrl <- crops_groups[sapply(crops_groups, match, x = CROP, nomatch = FALSE)==0] %>% unlist() %>% unname()  # this removes from the set of controls, the crops that are in the same group as CROP
    crops_ctrl <- commercial_crops
    
    # **NOPE** (control transitory LUC)
    # further remove woody perennials, if we regress cropland anyway
    crops_ctrl <- crops_ctrl[sapply(crops_ctrl, function(c){!(c %in% c(woody_perrenials))})]# , "eaear_Fodder"
    
    # the order of the strata in the list are determined by how we want to plot results
    temp_list_trade[[CNT]][[CROP]] <- make_main_reg(continent = CNT,
                                            outcome_variable = loss_type,
                                            
                                            start_year = est_parameters[["start_year"]],
                                            end_year = est_parameters[["end_year"]],
                                            
                                            # regress on one crop at a time, and control for other crops
                                            exposure_rfs = CROP, # "eaear_Maizegrain",
                                            control_all_absolute_rfs = TRUE,
                                            annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                            all_exposures_rfs = crops_ctrl,
                                            
                                            trade_exposure = est_parameters[["trade_exposure"]],
                                            trade_exposure_period = est_parameters[["trade_exposure_period"]],
                                            trade_expo_spec = est_parameters[["trade_expo_spec"]],
                                            
                                            rfs_lead = est_parameters[["leads"]],
                                            rfs_lag = est_parameters[["lags"]],
                                            aggr_lead = est_parameters[["aggr_lead"]],
                                            aggr_lag = est_parameters[["aggr_lag"]],
                                            sjpos = est_parameters[["sjpos"]],
                                            exposure_quantiles = est_parameters[["exposure_quantiles"]],
                                            s_trend = est_parameters[["s_trend"]],
                                            s_trend_loga = est_parameters[["s_trend_loga"]],
                                            s_trend_sq = est_parameters[["s_trend_sq"]],
                                            
                                            glm_iter = est_parameters[["glm_iter"]],
                                            
                                            output = "coef_table")# 
    
    rm(CROP)
  }
  
  ### OIL PALM (GROUP 3) ###
  oilpalm_loss <- "loss_oilpalm_indus"
  # oilpalm_ctrl <- cropland_crops
  # oilpalm_ctrl <- all_crops # [all_crops != "eaear_Coconut"]
  oilpalm_ctrl <- c()
  #for(oilpalm_loss in c("loss_oilpalm_both", "loss_oilpalm_indus")){
  temp_list_trade[[CNT]][["eaear_Oilpalm"]] <- make_main_reg(continent = CNT,
                                             outcome_variable = oilpalm_loss,
                                             start_year = est_parameters[["start_year"]],
                                             end_year = est_parameters[["end_year"]],
                                             
                                             # regress on crop at a time, but control for annual interactions with other crops
                                             exposure_rfs = "eaear_Oilpalm",
                                             control_all_absolute_rfs = TRUE,
                                             annual_rfs_controls = TRUE,
                                             all_exposures_rfs = oilpalm_ctrl,
                                             
                                             trade_exposure = est_parameters[["trade_exposure"]],
                                             trade_exposure_period = est_parameters[["trade_exposure_period"]],
                                             trade_expo_spec = est_parameters[["trade_expo_spec"]],
                                             
                                             rfs_lead = est_parameters[["leads"]],
                                             rfs_lag = est_parameters[["lags"]],
                                             aggr_lead = est_parameters[["aggr_lead"]],
                                             aggr_lag = est_parameters[["aggr_lag"]],
                                             sjpos = est_parameters[["sjpos"]],
                                             exposure_quantiles = est_parameters[["exposure_quantiles"]],
                                             s_trend = est_parameters[["s_trend"]],
                                             s_trend_loga = est_parameters[["s_trend_loga"]],
                                             s_trend_sq = est_parameters[["s_trend_sq"]],
                                             glm_iter = est_parameters[["glm_iter"]],
                                             
                                             output = "coef_table")
  
  

}

# Extract the term on TRADE exposure 
cnt_crop_list <- list()
cnt_list <- list()
for(CNT in continents[continents!="all"]){# "all", , "Africa", "Asia" [continents!="all"]
  if(CNT=="all"){
    commercial_crops <- commercial_crops[commercial_crops!="eaear_Oilfeed_crops"]
  } else{
    commercial_crops <- c(gentest_crops, displtest_crops)
    commercial_crops <- commercial_crops[!(commercial_crops %in% c("eaear_Biomass", "eaear_Roots"))] # 
  }
  
  for(CROP in c(commercial_crops, "eaear_Oilpalm")){

    x <- temp_list_trade[[CNT]][[CROP]] %>% as.data.frame()
    
    # save row names 
    estimate_names <- row.names(x)
    x$ESTIMATE_NAMES <- estimate_names
    
    x$model <- CNT
    x <- 
      x %>% 
      mutate(via_trade = case_when(
        grepl("_expo_", ESTIMATE_NAMES) & 
          !grepl("_all_but_j_and_maize", ESTIMATE_NAMES) & 
          !grepl("_expo_maize", ESTIMATE_NAMES) ~ "via_trade_j",
        grepl("_all_but_j_and_maize", ESTIMATE_NAMES) ~ "via_trade_not_j_not_maize",
        grepl("export_expo_maize", ESTIMATE_NAMES) ~ "via_export_maize",
        grepl("import_expo_maize", ESTIMATE_NAMES) ~ "via_import_maize",
        TRUE ~ "not_via_trade"
    )) 
    x$via_access <- grepl("_acc_", estimate_names)    
    x$crop <- CROP
    x$Dynamics <- ""
    x$Dynamics[grepl("aggrall", estimate_names)] <- "Cumulative"
    x$Dynamics[grepl("lead1", estimate_names)] <- "t+1"
    x$Dynamics[grepl("lead2", estimate_names)] <- "t+2"
    x$Dynamics[estimate_names==paste0(CROP,"_X_statute_conv")] <- "t"
    x$Dynamics[grepl("lag1", estimate_names)] <- "t-1"
    x$Dynamics[grepl("lag2", estimate_names)] <- "t-2"
    x$Dynamics[grepl("lag3", estimate_names)] <- "t-3"
    
    # order dynamics 
    # estimate_names <- x$Dynamics  
    # x <- x[c("Cumulative", "t+1", "t", "t-1", "t-2", "t-3", "t-4"),] 
    # # remove NA rows generated if some dynamics are absent for this crop 
    # x <- x[!is.na(x$tstat), ]
    row.names(x) <- NULL # otherwise they duplicate while stacking
    cnt_crop_list[[CNT]][[CROP]] <- x
  }
  cnt_list[[CNT]] <- cnt_crop_list[[CNT]] %>% bind_rows()
}
df <- cnt_list %>% bind_rows() %>% dplyr::select(-ESTIMATE_NAMES) # "loss_cropland"

# extract both the aggregated effect of AEAY*M and AEAY*M*TRADE
df <- df[df$Dynamics=="Cumulative",] ## CHANGE HERE TO "t+2" to display only anticipation effects

names(df)[names(df)=="t value"] <- "tstat"

df$significant01 <- FALSE
df[abs(df[,"tstat"]) > 1.645, "significant01"] <- TRUE

# attribute crop Group 
# horrible code mais j'ai pas réussi à faire mieux ... 
df_gn <- sapply(names(crops_groups), function(gn){if_else(df$crop %in% crops_groups[[gn]], true = gn, false = NULL)}) 
df$Group <- ""
for(i in 1:nrow(df_gn)){
  row <- df_gn[i,]
  df$Group[i] <- row[!is.na(row)] %>% unname()
}

# make some variable name changes necessary for dotwhisker
# df <- df[df$crop != "eaear_Maizegrain",]
df <- mutate(df, sizes = if_else(model == "all", true = 4, false = 2.5))
names(df)[names(df)=="tstat"] <- "estimate"
# names(df)[names(df)=="model"] <- "continent"
# names(df)[names(df)=="Dynamics"] <- "model"
names(df)[names(df)=="crop"] <- "term"
# df$cumulative <- factor(df$Dynamics=="Cumulative")
head(df)


dwplot(df,
       dot_args = list(aes(color = via_trade, shape = model, size = model, alpha = significant01)) ) %>%  #size = 2.5,  shape = Dynamics, size = cumulative
  relabel_predictors(predictors_dict) + 
  # # critical values  : 1.645   1.960  2.576
  geom_vline(xintercept = 1.96, colour = "grey60", linetype = "dotdash", alpha = 1) +
  geom_vline(xintercept = -1.96, colour = "grey60", linetype = "dotdash", alpha = 1) +
  geom_vline(xintercept = 2.576, colour = "grey60", linetype = "dotted", alpha = 1) +
  geom_vline(xintercept = -2.576, colour = "grey60", linetype = "dotted", alpha = 1) +
  
  facet_wrap(facets = ~Group, scales = "free_y") + 
  
  scale_color_brewer(type = "qual",palette="Dark2", #  Accent  
                     breaks=c("via_trade_j", "via_export_maize", "via_import_maize", "via_trade_not_j_not_maize", "not_via_trade"), 
                     labels=c("export of the crop", 
                              "export of maize",
                              "import of maize",
                              "export of other crops", 
                              "none"), # "Moderated by trade exposure", "Irrespective of trade exposure"
                     name="Moderated by 2001-2007 exposure to") +
  
  scale_shape_manual(values = c(17, 15, 16), #18,  c(21, 24, 22, 25),# c(1, 2, 0, 3), # c(16, 17, 15, 18)
                     breaks=c("America", "Africa", "Asia"),# "all", 
                     labels = c("America (tropical lat.)", "Africa (tropical lat.)", "Asia (tropical lat.)"),# "Pan-tropical", 
                     name="") +
  
  scale_size_manual(values = c(3.5, 2.5, 2.5, 2.5),
                    breaks=c("all", "America", "Africa", "Asia"),
                    labels = c("Pan-tropical", "America (tropical lat.)", "Africa (tropical lat.)", "Asia (tropical lat.)"),
                    name="") +
  
  scale_alpha_manual(# values = sort(unique(df$highlight)),# c(0.33, 0.50, 0.67, 0.83, 1.00), # c(0.2, 0.4, 0.6, 0.8, 1),
    values = c(0.3, 1) ,
    breaks= c(FALSE, TRUE),# sort(unique(df$highlight)),# 
     labels = c(" ", "p-value < .1"),
  name=" ") +
  
  scale_x_continuous(breaks = c(-1.96, 1.96), 
                     labels = c("-1.96","1.96")) +
  
  theme_bw() + xlab("t-statistics") + ylab("") +  
  theme(plot.title = element_text(face="bold", size=c(10)),
        legend.position = "bottom",# c(0.8, 0.05),
        legend.justification = c(0, 0), 
        legend.background = element_rect(colour="grey80")
  ) # legend.title = element_blank()


rm(cnt_list, cnt_crop_list)

#### EAEAR CORRELATION MATRIX ####

# Work on rasters directly? NOP because we want to compare with crop groups and derive tests. 
# croped at global scale (under 60°)

# Irrigated only
gaez_dir_global <- here("temp_data", "GAEZ", "v4", "AEAY_out_density", "Rain-fed-all-phases", "global_AOI")
gaez_crops <- list.files(path = here(gaez_dir_global, "High-input"), 
                         pattern = "", 
                         full.names = FALSE)
gaez_crops <- gsub(pattern = ".tif", replacement = "", x = gaez_crops)

gaez_global <- brick(here(gaez_dir_global, "high_input_all.tif"))
names(gaez_global) <- gaez_crops

gaez_crops <- names(gaez_global)

# Make a data frame
wide_df <- raster::as.data.frame(gaez_global, na.rm = TRUE, xy = TRUE, centroids = TRUE, long = FALSE) # ~700s. 

# Rename coordinate variables
names(wide_df)
head(wide_df[,c("x", "y")])
wide_df <- dplyr::rename(wide_df, lon = x, lat = y)

# Group
mapmat_data <- c(
  "Banana","Banana",
  "Barley", "Barley",
  "Crude_oil", "Biomass", # these crop categories are gonna be created in the present script
  "Orange", "Citrus", # Citrus sinensis in both GAEZ and FAO
  "Cocoa", "Cocoa",
  "Coconut_oil", "Coconut", # Coconut not excluded as we don't use prices anymore. See below, in conversion part, why we would exclude it if we needed price scaling 
  "Coffee", "Coffee",
  "Cotton", "Cotton",
  "Beef", "Fodder", # these crop categories are gonna be created in the present script 
  "Groundnuts", "Groundnut",
  "Maize", "Maizegrain",
  "Oat", "Oat",
  "Olive_oil", "Olive",  
  "Palm_oil", "Oilpalm",
  "Rapeseed_oil", "Rapeseed",
  "Rice", "Rice",
  "Roots", "Roots", # these crop categories are gonna be created in the present script 
  "Rubber", "Rubber",
  "Sorghum", "Sorghum", # this will be matched with barley and wheat, i.e. grains we have price data on.
  "Soybean", "Soybean",
  "Soybean_meal", "Soybean_meal", # these crop categories are gonna be created in the present script
  "Soybean_oil", "Soybean_oil", # these crop categories are gonna be created in the present script
  "Sugar", "Sugar", # these crop categories are gonna be created in the present script
  "Sugar", "Sugarbeet",
  "Sugar", "Sugarcane",
  "Sunflower_oil", "Sunflower",
  "Tea", "Tea",
  "Tobacco", "Tobacco", 
  "Wheat", "Wheat")

mapmat <- matrix(data = mapmat_data, 
                 nrow = length(mapmat_data)/2,
                 ncol = 2, 
                 byrow = TRUE)

colnames(mapmat) <- c("Prices", "Crops")

# crops to group based on potential REVENUE
crops2grp <- c("Barley", "Sorghum", "Wheat", "Cocoa", "Coffee", "Groundnut", "Rapeseed", "Sunflower")

# crops to standardize. There is not fodder, rubber, citrus, banana, cocoa, coffee, olive and tea
eaear2std <- paste0("eaear_", c("Cereals", "Oilfeed_crops", "Cotton", "Maizegrain", "Oat", "Oilpalm", "Rice",
                                "Soy_compo", "Sugar", "Tobacco")) 
# add cocoa, coffee and tea for std2 
eaear2std_bis <- paste0("eaear_", c("Banana", "Biomass", "Cereals", "Oilfeed_crops", "Cocoa_Coffee", "Cotton", 
                                    "Maizegrain", "Oat", "Olive", "Oilpalm", "Rice",
                                    "Soy_compo", "Sugar", "Tea", "Tobacco")) 
# Group 
wide_df <- wide_df %>% rowwise() %>% mutate(Biomass = max(c(Sorghumbiomass/10, Miscanthus, Reedcanarygrass, Switchgrass)), # sorghum biomass is expressed in kg/ha, and not 10kg/ha as it is the case for the three other crops
                                        #  don't include Jatropha because it is expressed in seeds and not above ground biomass in GAEZ. 
                                        Fodder = max(c(Alfalfa, Napiergrass)),   #,  # it's almost only napiergrass that varies (indeed, some say that it's the highest yielding tropical forage crop https://www.sciencedirect.com/science/article/pii/S2468227619307756#bib0001
                                        Rice = max(c(Drylandrice, Wetlandrice)),
                                        Roots = max(c(Cassava, Sweetpotato, Whitepotato,Yam)), 
                                        Sugar = max(c(Sugarbeet, Sugarcane)) # Especially necessary to match the international price of sugar
                                       ) %>% as.data.frame()


# keep only crops of interest 
wide_df <- wide_df[, names(wide_df) %in% c("grid_id", "lon", "lat", mapmat[,"Crops"])]

conv_fac <- c(Wheat = 0.87, 
              Rice = 0.87, 
              Maizegrain = 0.86, 
              Sorghum = 0.87, 
              Barley = 0.87, 
              Oat = 0.87, 
              Soybean = 0.90, 
              Rapeseed = 0.90, 
              Sunflower = 0.92, 
              Groundnut = 0.65, # this is applied to go from shelled groundnuts to GAEZ dry weight 
              # Cotton = 0.33, we won't use this, see below paragraph on cotton
              Banana = 0.25, # from banana to dry weight banana
              Tobacco = 0.75) # from traded tobacco dry leaves to GAEZ dry weight (i.e. traded leaves have water content of 25%)  
# get prices for 2000
prices$Roots <- 1 # this is just for convenience, as we do not compare commodity values (multiplied by prices) with each others now. 
price_avg <- prices %>% 
  filter(year>=1995 & year <= 2004) %>% 
  #filter(year==2000) %>% 
  summarise(across(.cols = any_of(mapmat[,"Prices"]), 
                   .fns = mean, na.rm = TRUE))

df_cs <- wide_df
df_cs <- dplyr::mutate(df_cs, Banana = Banana / conv_fac["Banana"])
df_cs <- dplyr::mutate(df_cs, Barley = Barley / conv_fac["Barley"])
df_cs <- dplyr::mutate(df_cs, Fodder = Fodder * 0.05 / 0.271)
df_cs <- dplyr::mutate(df_cs, Groundnut = Groundnut / conv_fac["Groundnut"])
df_cs <- dplyr::mutate(df_cs, Maizegrain = Maizegrain / conv_fac["Maizegrain"])
df_cs <- dplyr::mutate(df_cs, Oat = Oat / conv_fac["Oat"])
df_cs <- dplyr::mutate(df_cs, Rapeseed = Rapeseed * 0.41 / conv_fac["Rapeseed"])
df_cs <- dplyr::mutate(df_cs, Rubber = Rubber / 0.35)
df_cs <- dplyr::mutate(df_cs, Sorghum = Sorghum / conv_fac["Sorghum"])
df_cs <- dplyr::mutate(df_cs, Sunflower = Sunflower * 0.42 / conv_fac["Sunflower"])
df_cs <- dplyr::mutate(df_cs, Soybean = Soybean / conv_fac["Soybean"])
df_cs <- dplyr::mutate(df_cs, Soybean_oil = Soybean * 0.18 / conv_fac["Soybean"])
df_cs <- dplyr::mutate(df_cs, Soybean_meal = Soybean * 0.82 / conv_fac["Soybean"])
df_cs <- dplyr::mutate(df_cs, Rice = Rice / conv_fac["Rice"])
df_cs <- dplyr::mutate(df_cs, Tobacco = Tobacco / conv_fac["Tobacco"])
df_cs <- dplyr::mutate(df_cs, Wheat = Wheat / conv_fac["Wheat"])
df_cs <- dplyr::mutate(df_cs, across(.cols = all_of(mapmat[,"Crops"]),
                                     .fns = ~./1000)) 
df_cs <- dplyr::mutate(df_cs, across(.cols = all_of(c("Fodder", "Biomass")),
                                     .fns = ~.*10)) 
for(aeay_i in mapmat[,"Crops"]){
  price_i <- price_avg[mapmat[mapmat[,"Crops"]==aeay_i,"Prices"]]%>%as.numeric()
  eaear_i <- paste0("eaear_", aeay_i)
  df_cs <- dplyr::mutate(df_cs, 
                         !!as.symbol(eaear_i) := !!as.symbol(aeay_i) * price_i)
}
df_cs <- df_cs %>% rowwise() %>% mutate(eaear_Cereals = max(c(eaear_Barley, eaear_Sorghum, eaear_Wheat)), 
                                        eaear_Oilfeed_crops = max(c(eaear_Groundnut, eaear_Rapeseed, eaear_Sunflower)), 
                                        eaear_Cocoa_Coffee = max(c(eaear_Cocoa, eaear_Coffee))) %>% as.data.frame()

df_cs <- dplyr::mutate(df_cs, eaear_Soy_compo =  eaear_Soybean_meal + eaear_Soybean_oil)
df_cs <- dplyr::select(df_cs, -eaear_Soybean, -eaear_Soybean_meal, -eaear_Soybean_oil)


# for(others in all_crops[all_crops!="eaear_Maizegrain"]){
#   
# }
# jnk <- layerStats(gaez_global[[c("Maizegrain", "Tobacco", "Oilpalm")]], 'pearson', na.rm=T)

# Restrict to places with positive AEAY for maize
# which is distributed like this
gaez_maize <- gaez_global[["Maizegrain"]]
plot(gaez_maize)

# For the global correlation exercise, we use Sugar (beet and cane) and not only sugarcane as in the regressions. 
names_4corr <- names(predictors_dict)
names_4corr[names_4corr=="eaear_Sugarcane"] <- "eaear_Sugar"

df_cs_posmaiz <- 
  df_cs %>% 
  filter(eaear_Maizegrain>0)
cor_mat_abs <-  cor(dplyr::select(df_cs_posmaiz, all_of(names_4corr) ))#starts_with("eaear_")
# for comparison: 
# wide_df_posmaiz <- 
#   wide_df %>% 
#   filter(Maizegrain>0)
# cor_mat_abs_raw <-  cor(dplyr::select(wide_df_posmaiz, all_of(c("Maizegrain", "Cotton", "Tobacco", "Citrus", "Coconut", "Banana", "Tea", "Rubber")) ))#starts_with("eaear_")

cortests <- cor_mat_abs
j_exposures <- list()
length(j_exposures) <- length(predictors_dict)
names(j_exposures) <- names_4corr

for(serie1 in colnames(cortests)){
  
  for(serie2 in row.names(cortests)){
    
    cortests[serie1, serie2] <- cor.test(df_cs_posmaiz[,serie1], df_cs_posmaiz[,serie2])$p.value
  }
  j_exposures[[serie1]] <- cor_mat_abs[cortests[serie1, ] < 0.05, serie1]
}

j_exposures[["eaear_Maizegrain"]] 

totable <- cor_mat_abs["eaear_Maizegrain", ] %>% as.data.frame() # %>% matrix(nrow = 1)
names(totable) <- "Corr"

for(coln in rownames(totable)){
  rownames(totable)[rownames(totable) == coln] <- predictors_dict[names_4corr==coln] #%>% unname()
}
# it spuriously gave back sugarcane name, so change it again
rownames(totable)[rownames(totable) == "Sugarcane"] <- "Sugar crops" 

totable <- 
  totable %>% 
  filter(rownames(totable)!="Maize") %>% 
  arrange(desc(Corr)) %>%
  mutate(Corr = round(Corr, 2)) 

names(totable) <- NULL
options(knitr.table.format = "latex") 
kable(totable, booktabs = T, align = "c", 
      caption = "Coefficients of correlation of crop AEAYs with maize AEAY") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) 


# store, for each crop, which other crop it is most correlated with
# corr_mapmat <- cbind(names(predictors_dict),NA)
# colnames(corr_mapmat) <- c("Crops", "fst_corr")
# for(crop in names(predictors_dict)){
#   x <- j_exposures[[crop]]
#   corr_mapmat[corr_mapmat[,"Crops"]==crop, "fst_corr"] <- names(x)[x == max(x[x<max(x)])] 
# }


#### MAP OUTCOMES  -------------------------------------------------
# it makes sense to plot outcomes (deforestation for cropland and for oil palm) as:
# maps, cumulative over time, in percentage of grid area (more explicit than absolute number of hectares)

# first, get the sample actually used for estimation 
d_clean_list <- list()
for(LT in c("loss_cropland", "loss_oilpalm_indus")){
  d_clean_out <- make_main_reg(continent = "all",
                              outcome_variable = LT,
                              
                              start_year = est_parameters[["start_year"]],
                              end_year = est_parameters[["end_year"]],
                              
                              # regress on one crop at a time, and control for other crops
                              exposure_rfs = "eaear_Maizegrain", # whatever
                              control_all_absolute_rfs = TRUE,
                              annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                              all_exposures_rfs = c(), # whatever
                              
                              trade_exposure = est_parameters[["trade_exposure"]],
                              trade_exposure_period = est_parameters[["trade_exposure_period"]],
                              
                              rfs_lead = est_parameters[["leads"]],
                              rfs_lag = est_parameters[["lags"]],
                              aggr_lead = est_parameters[["aggr_lead"]],
                              aggr_lag = est_parameters[["aggr_lag"]],
                              sjpos = est_parameters[["sjpos"]],
                              exposure_quantiles = est_parameters[["exposure_quantiles"]],
                              s_trend = est_parameters[["s_trend"]],
                              s_trend_loga = est_parameters[["s_trend_loga"]],
                              s_trend_sq = est_parameters[["s_trend_sq"]],
                              
                              glm_iter = est_parameters[["glm_iter"]],
                              
                              output = "everything")
  
  d_clean_out <- d_clean_out[[2]] # data is in 2nd element 
  
  d_clean_list[[LT]] <- ddply(d_clean_out, "grid_id", summarise, 
                              accu := sum(!!as.symbol(LT)), 
                              lon = unique(lon), 
                              lat = unique(lat), 
                              continent_name = unique(continent_name)) 
  
  d_clean_list[[LT]]$loss_type <- LT
}

d_clean_accu <- bind_rows(d_clean_list)

# spatialize
d_clean_sf <- st_as_sf(d_clean_accu, coords = c("lon", "lat"), crs = 4326)
# make grid shapes
d_clean_sf <- st_transform(d_clean_sf, crs = mercator_world_crs) 
d_clean_sf <- st_buffer(d_clean_sf, dist = 4500) # half the size of a grid (approximately)
st_geometry(d_clean_sf) <- sapply(st_geometry(d_clean_sf), FUN = function(x){st_as_sfc(st_bbox(x))}) %>% st_sfc(crs = mercator_world_crs)
# go back to unprojected crs
d_clean_sf <- st_transform(d_clean_sf, crs = 4326)

# accumulated deforestation in hectares (i.e. 10000m2) to percentage of grid cell area 
# but don't take the actual area of the shape, because it's made up by the construction above. as.numeric(st_area(geometry))
# So just take 9km2 for every cell, it's just an approximation for visualization 
d_clean_sf <- mutate(d_clean_sf, accu_pct = 100 * accu * 10000 / 9000^2) # first multiply by 100 to get pct point
summary(d_clean_sf$accu_pct)

# make three different datasets 
am <- d_clean_sf[d_clean_sf$continent_name == "America", ]
af <- d_clean_sf[d_clean_sf$continent_name == "Africa", ]
as <- d_clean_sf[d_clean_sf$continent_name == "Asia", ]

land <- st_read(here("input_data", "ne_50m_land"))
unique(land$scalerank)
land <- land[land$scalerank==0, c("geometry")]
#plot(land)
# spLand <- as(land, "Spatial")

d_clean_sf <- mutate(d_clean_sf, loss_type = if_else(loss_type=="loss_cropland", 
                                                     true = "Commodity-driven cropland deforestation",
                                                     false = "Commodity-driven oil palm deforestation"))
# d_clean_sf <- mutate(d_clean_sf, accu_pct_d = cut_width(x = accu_pct, width = 0., label = FALSE))
d_clean_sf <- mutate(d_clean_sf, accu_pct_d = cut(x = accu_pct, 
                                                  breaks = c(0, 0.01, 0.1, 1, 10, 79), 
                                                  right = T, 
                                                  include.lowest = F, 
                                                  label = c("< 0.01%", 
                                                            "0.01-0.1%",
                                                            "0.1-1%",
                                                            "1-10%",
                                                            "> 10%")))

summary(d_clean_sf$accu_pct_d)
d_clean_sf$accu_pct_d <- as.character(d_clean_sf$accu_pct_d)
unique(d_clean_sf$accu_pct_d)
min(d_clean_sf$accu)==0
sum(d_clean_sf$accu)

base_map <- ggplot(data = land) +
  geom_sf(alpha = 0) + theme_bw() +# alpha = 0.8fill = "lightgrey", 
  geom_sf(data = d_clean_sf, aes(fill = accu_pct), lwd = NA) + # 
  coord_sf(xlim = c(-95, 150), ylim = c(-30, 30), expand = FALSE) 

base_map + scale_fill_viridis(name = "Cumul. 2010-2016\nin % of cell area", 
                              option="viridis",
                              direction = -1) + 
  facet_wrap(facets = ~loss_type, ncol = 1, nrow = 2, 
             strip.position = "top") +
  
  scale_y_discrete(breaks = c(0))

### WIP ### 

ggplot(data = land) +
  geom_sf(alpha = 0) + theme_bw() +# alpha = 0.8fill = "lightgrey", 
  geom_sf(data = d_clean_sf, aes(fill = accu_pct_d), lwd = NA) + # 
  coord_sf(xlim = c(-95, 150), ylim = c(-30, 30), expand = FALSE) +
  
  scale_fill_viridis(discrete = TRUE# name = "Cumul. 2010-2016\nin % of cell area",
                    #option = "viridis", 
                    #direction = -1, 
                    #breaks = c(0.1, 1, 10, 79), # unique(d_clean_sf$accu_pct_d), 
                    ) +

                     # breaks = c(0.0001, 0.001, 0.01, 0.78)
  labels = c("< 0.01%", 
                               "0.01-0.1%",
                               "0.1-1%",
                               "1-10%",
                               "> 10%")
  facet_wrap(facets = ~loss_type, ncol = 1, nrow = 2, 
             strip.position = "top") +
  scale_y_discrete(breaks = c(0))


ggplot(data = land) +
  geom_sf(fill = "lightgrey") + theme_minimal() +# alpha = 0.8, alpha = 0.5
  geom_sf(data = d_clean_sf, aes(fill = accu_pct), lwd = NA) + # 
  #coord_sf(xlim = c(-95, -30), ylim = c(-33, 25), expand = FALSE) +
  scale_fill_viridis(name = '% cropland deforestation\n in cell area', 
                     option="viridis",
                     direction = -1) +
  facet_zoom(xlim = c(-95, -30)) +
  
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())

  matrix(c(-95, 25, -33, 25,
           -33, -30, -95, -30, 
           -95, 25), ncol = 2, byrow = TRUE)

  rm(d_clean, d_clean_accu, d_clean_list, d_clean_out, d_clean_sf, land, base_map, af, am, as)

  

  #### GROUP 6 SEPARATELY #### 
  rm(all_tests_est, d_clean_list)
  all_tests_est <- list()
  d_clean_list <- list()
  for(CNT in continents){
    ### GROUP 6 CROPS 
    loss_type <- "loss_pasture"
    
    for(CROP in c("eaear_Fodder", "eaear_Citrus", "eaear_Coconut", g6_crops)){
      ## Determine the crops to control for
      #crops_ctrl <- crops_groups[sapply(crops_groups, match, x = CROP, nomatch = FALSE)==0] %>% unlist() %>% unname()  # this removes from the set of controls, the crops that are in the same group as CROP
      crops_ctrl <- all_crops
      
      # the order of the strata in the list are determined by how we want to plot results
      temp_list_g6 <- make_main_reg(continent = CNT,
                                    outcome_variable = loss_type,
                                    
                                    start_year = est_parameters[["start_year"]],
                                    end_year = est_parameters[["end_year"]],
                                    
                                    # regress on one crop at a time, and control for other crops
                                    exposure_rfs = CROP, # "eaear_Maizegrain",
                                    control_all_absolute_rfs = TRUE,
                                    annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                    all_exposures_rfs = crops_ctrl,
                                    
                                    trade_exposure = est_parameters[["trade_exposure"]],
                                    trade_exposure_period = est_parameters[["trade_exposure_period"]],
                                    
                                    rfs_lead = est_parameters[["leads"]],
                                    rfs_lag = est_parameters[["lags"]],
                                    aggr_lead = est_parameters[["aggr_lead"]],
                                    aggr_lag = est_parameters[["aggr_lag"]],
                                    sjpos = est_parameters[["sjpos"]],
                                    exposure_quantiles = est_parameters[["exposure_quantiles"]],
                                    s_trend = est_parameters[["s_trend"]],
                                    s_trend_loga = est_parameters[["s_trend_loga"]],
                                    s_trend_sq = est_parameters[["s_trend_sq"]],
                                    
                                    glm_iter = est_parameters[["glm_iter"]],
                                    
                                    output = "everything")
      
      all_tests_est[[CNT]][[CROP]] <- temp_list_g6[[1]] # this is the fixest object
      d_clean_list[[CNT]][[loss_type]] <- temp_list_g6[[2]] # this is the data (always the same at every CROP iteration)
      rm(temp_list_g6)
      rm(CROP)
    }
  }
  for(CNT in "all"){# "all", , "Africa", "Asia"
    
    for(CROP in c(c("eaear_Citrus", "eaear_Coconut", g6_crops))){
      
      loss_type <- "loss_pasture"  
      d_clean <- d_clean_list[[CNT]][[loss_type]]
      
      # store post-estimation outputs of all clustering levels
      postest_list <- list()
      # THREE LEVELS OF SPATIAL CLUSTERING
      # Main 
      postest_list[["grid_id_10"]] <- post_est_fnc(est_obj = all_tests_est[[CNT]][[CROP]], # loss_type
                                                   base_exposure = CROP,
                                                   trade_exposure = NULL,
                                                   aggr_lag = est_parameters[["aggr_lag"]], 
                                                   aggr_lead = est_parameters[["aggr_lead"]], 
                                                   clustering = "oneway",
                                                   cluster_var1 = "grid_id_10",
                                                   output = "tstat")
      # 5 times larger
      postest_list[["grid_id_5"]] <- post_est_fnc(est_obj = all_tests_est[[CNT]][[CROP]], # loss_type
                                                  base_exposure = CROP,
                                                  trade_exposure = NULL,
                                                  aggr_lag = est_parameters[["aggr_lag"]], 
                                                  aggr_lead = est_parameters[["aggr_lead"]], 
                                                  clustering = "oneway",
                                                  cluster_var1 = "grid_id_5",
                                                  output = "tstat")
      # 20 times larger
      postest_list[["grid_id_20"]] <- post_est_fnc(est_obj = all_tests_est[[CNT]][[CROP]], # loss_type
                                                   base_exposure = CROP,
                                                   trade_exposure = NULL,
                                                   aggr_lag = est_parameters[["aggr_lag"]], 
                                                   aggr_lead = est_parameters[["aggr_lead"]], 
                                                   clustering = "oneway",
                                                   cluster_var1 = "grid_id_20",
                                                   output = "tstat")
      
      ## TWO-WAY  CLUSTERING
      # need to do this here and not in post_est_fnc for fixest to find the data to construct the new vcov 
      d_clean <- left_join(d_clean, main_data[,c("grid_id", "year", CROP)], by = c("grid_id", "year"))
      if(is.numeric(est_parameters[["cluster_var1"]])){
        # make sure that the exposure var is not there already
        d_clean <- dplyr::select(d_clean, -!!as.symbol(CROP[1]))
        d_clean <- mutate(d_clean, !!as.symbol(paste0(CROP[1],"_",est_parameters[["cluster_var1"]],"tiles")) := cut(!!as.symbol(CROP[1]), breaks = est_parameters[["cluster_var1"]], labels = paste0(CROP[1],"_Q", 1:est_parameters[["cluster_var1"]])))
      }
      # extract t stats of coef of interest 
      postest_list[["twoway_percentile"]] <- post_est_fnc(est_obj = all_tests_est[[CNT]][[CROP]], # loss_type
                                                          base_exposure = CROP,
                                                          trade_exposure = NULL,
                                                          aggr_lag = est_parameters[["aggr_lag"]], 
                                                          aggr_lead = est_parameters[["aggr_lead"]],
                                                          clustering = "twoway",
                                                          cluster_var1 = est_parameters[["cluster_var1"]],
                                                          cluster_var2 = est_parameters[["cluster_var2"]],
                                                          output = "tstat")
      rm(d_clean)
      
      # save row names 
      estimate_names <- row.names(postest_list[[1]])
      
      x <- bind_cols(postest_list, .name_repair = "minimal")
      names(x) <- names(postest_list)
      x <- dplyr::mutate(x, across(.cols = everything(), 
                                   ~abs(.)>1.645, 
                                   .names = paste0("{.col}","_isprecise")))
      x <- dplyr::mutate(x, rob_preci = rowMeans(across(.cols = contains("_isprecise")) ))
      # keep only the main t-stat
      x <- dplyr::select(x, -grid_id_5, -grid_id_20, -twoway_percentile)
      names(x)[names(x)=="grid_id_10"] <- "tstat"
      
      x$model <- CNT
      x$via_trade <- grepl("_expo_", estimate_names)    
      x$crop <- CROP
      x$Dynamics <- ""
      x$Dynamics[grepl("aggrall", estimate_names)] <- "Cumulative"
      x$Dynamics[grepl("lead1", estimate_names)] <- "t+1"
      x$Dynamics[grepl("lead2", estimate_names)] <- "t+2"
      x$Dynamics[estimate_names==paste0(CROP,"_X_statute_conv")] <- "t"
      x$Dynamics[grepl("lag1", estimate_names)] <- "t-1"
      x$Dynamics[grepl("lag2", estimate_names)] <- "t-2"
      x$Dynamics[grepl("lag3", estimate_names)] <- "t-3"
      
      # order dynamics 
      # estimate_names <- x$Dynamics  
      # x <- x[c("Cumulative", "t+1", "t", "t-1", "t-2", "t-3", "t-4"),] 
      # # remove NA rows generated if some dynamics are absent for this crop 
      # x <- x[!is.na(x$tstat), ]
      row.names(x) <- NULL # otherwise they duplicate while stacking
      cnt_crop_list[[CNT]][[CROP]] <- x
      
    }
    cnt_list[[CNT]] <- cnt_crop_list[[CNT]] %>% bind_rows()
  }
  df <- cnt_list %>% bind_rows()
  
#### LOSS CATEGORIES ------------------------------------------

# what's different here is that regressions are run with output = coef_table, 
# CROP by CROP in order to have crop specirfic trends and/or interactions with trade exposures. 
# but not with output = est_object otherwise eats all memory. 
# and with data pre processing to spare time 

### DATA 
# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "loss_cropcommo_op_aeaycompo_long_final.Rdata"))
main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "loss_crop_opindus_aeaycompo_long_final.Rdata"))
# this is with cropland outside commo, and without pasture
#main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "loss_aeaycompo_long_final.Rdata"))

#main_data <- dplyr::mutate(main_data, loss_commodity = loss_cropland + loss_oilpalm_both + loss_pasture)
# release some memory upfront
main_data <- dplyr::filter(main_data, year >= 2008, year <= 2019)

# this may be used to remove some controls in loss_cropland regressions
woody_perrenials <- c("eaear_Citrus", "eaear_Coconut", "eaear_Oilpalm", 
                      "eaear_Banana", "eaear_Cocoa_Coffee", "eaear_Tea", "eaear_Rubber")

# this may be used to remove some controls in loss residual regressions (called loss_pasture)
cropland_crops <- c(all_crops[!(all_crops %in% c(woody_perrenials))])

# loss_type <- "loss_cropland"
# CNT <- "America"
# CROP <- "eaear_Soy_compo"

all_tests_est <- list()
for(CNT in continents){
  
  ### PRE PROCESS ### 
  d <- main_data
  
  # - are in study period
  d <- dplyr::filter(d, year >= est_parameters[["start_year"]])
  d <- dplyr::filter(d, year <= est_parameters[["end_year"]])
  
  # - are in study area
  if(CNT != "all"){
    d <- dplyr::filter(d, continent_name == CNT)
  }
  
  # Prepare interactions with trade, if specified
  all_eaear_trade_exposure_rfs <- c()
  if(length(est_parameters[["trade_exposure"]])>0){
    
    trade_expo_dat <- readRDS(here("temp_data", "processed_trade_exposures", paste0("trade_exposures_", est_parameters[["trade_exposure_period"]],".Rdata")))
    # if trade_exposure is not _imp, the following line will keep both expo and expo_imp variables, but it already divides by 2 the number of cols that are joined to the bigger d dataframe
    trade_expo_dat <- trade_expo_dat[, (grepl(x = names(trade_expo_dat), pattern = est_parameters[["trade_exposure"]]) | names(trade_expo_dat)=="country_name") ]
    
    d <- left_join(d, trade_expo_dat, by = "country_name")
    
    # make ALL eaear-trade exposures
    for(exp_rfs in est_parameters[["all_exposures_rfs"]]){ 
      # identify the corresponding trade_expo 
      trade_expo_j <- names(d)[names(d) == paste0(est_parameters[["trade_exposure"]],"_",str_to_lower(gsub(pattern = "eaear_", replacement = "", x = exp_rfs) ) ) ]
      eaear_trade_expo_j <- paste0("eaear_",trade_expo_j)
      d <- mutate(d, !!as.symbol(eaear_trade_expo_j) := !!as.symbol(exp_rfs) * !!as.symbol(trade_expo_j))
      all_eaear_trade_exposure_rfs <- c(all_eaear_trade_exposure_rfs, eaear_trade_expo_j)
    }
    rm(trade_expo_dat)
  }
  
  # Keep only in data the useful variables 
  d <- dplyr::select(d, all_of(unique(c("grid_id", "year", "lat", "lon", "continent_name", "country_name", "country_year",
                                        # "fc_2000", "fc_2009", "remaining_fc", # "accu_defo_since2k",
                                        "grid_id_5", "grid_id_10", "grid_id_20", "grid_id_5_year", "grid_id_10_year", "grid_id_20_year",
                                        # "tmf_agri", "tmf_flood", "tmf_plantation", "tmf_deforestation",
                                        "loss_cropland", "loss_oilpalm_indus",
                                        "pasture_share_2000",
                                        all_crops, all_eaear_trade_exposure_rfs ))))
  
  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", unique(c(all_rfs_treatments)))], by = c("year"))
  
  # d[d[,est_parameters[["outcome_variable"]]] < est_parameters[["grid_event_threshold"]], est_parameters[["outcome_variable"]]] <- 0
  
  ### CROPLAND ### 
  loss_type <- "loss_cropland"

  temp_est <- feglm(fml = as.formula(paste0(loss_type, " ~ 1 | grid_id + country_year")),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    pre_d_clean_cropland <- d[unlist(temp_est$obs_selection),]
  }  else {
    pre_d_clean_cropland <- d
  }
  rm(temp_est)

  # GROUP 1-2
  for(CROP in c(gentest_crops, displtest_crops)){
    ## Determine the crops to control for
    #crops_ctrl <- crops_groups[sapply(crops_groups, match, x = CROP, nomatch = FALSE)==0] %>% unlist() %>% unname()  # this removes from the set of controls, the crops that are in the same group as CROP
    crops_ctrl <- all_crops

    # **NOPE** (control transitory LUC)
    # further remove woody perennials and pastures, if we regress cropland anyway
    crops_ctrl <- crops_ctrl[sapply(crops_ctrl, function(c){!(c %in% c(woody_perrenials))})]# , "eaear_Fodder"

    # the order of the strata in the list are determined by how we want to plot results
    all_tests_est[[CNT]][[CROP]] <- make_main_reg(pre_process = TRUE,
                                                  pre_processed_data = pre_d_clean_cropland,
                                                  continent = CNT,
                                                  outcome_variable = loss_type,

                                                  start_year = est_parameters[["start_year"]],
                                                  end_year = est_parameters[["end_year"]],

                                                  # regress on one crop at a time, and control for other crops
                                                  exposure_rfs = CROP, # "eaear_Maizegrain",
                                                  control_all_absolute_rfs = TRUE,
                                                  annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                  all_exposures_rfs = crops_ctrl,

                                                  trade_exposure = est_parameters[["trade_exposure"]],
                                                  trade_exposure_period = est_parameters[["trade_exposure_period"]],

                                                  rfs_lead = est_parameters[["leads"]],
                                                  rfs_lag = est_parameters[["lags"]],
                                                  aggr_lead = est_parameters[["aggr_lead"]],
                                                  aggr_lag = est_parameters[["aggr_lag"]],
                                                  sjpos = est_parameters[["sjpos"]],
                                                  s_trend = est_parameters[["s_trend"]],
                                                  s_trend_loga = est_parameters[["s_trend_loga"]],
                                                  s_trend_sq = est_parameters[["s_trend_sq"]],

                                                  glm_iter = est_parameters[["glm_iter"]],

                                                  output = "coef_table")
    rm(CROP)
  }
  
  ### OIL PALM (GROUP 3) ###
  oilpalm_loss <- "loss_oilpalm_indus"

  temp_est <- feglm(fml = as.formula(paste0(oilpalm_loss, " ~ 1 | grid_id + country_year")),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    pre_d_clean_op <- d[unlist(temp_est$obs_selection),]
  }  else {
    pre_d_clean_op <- d
  }
  rm(temp_est)

  oilpalm_ctrl <- cropland_crops
  # oilpalm_ctrl <- all_crops # [all_crops != "eaear_Coconut"]
  # oilpalm_ctrl <- c()
  #for(oilpalm_loss in c("loss_oilpalm_both", "loss_oilpalm_indus")){
  all_tests_est[[CNT]][["eaear_Oilpalm"]] <- make_main_reg(pre_process = TRUE,
                                                           pre_processed_data = pre_d_clean_op,
                                                           continent = CNT,
                                                           outcome_variable = oilpalm_loss,
                                                           start_year = est_parameters[["start_year"]],
                                                           end_year = est_parameters[["end_year"]],

                                                           # regress on crop at a time, but control for annual interactions with other crops
                                                           exposure_rfs = "eaear_Oilpalm",
                                                           control_all_absolute_rfs = TRUE,
                                                           annual_rfs_controls = TRUE,
                                                           all_exposures_rfs = oilpalm_ctrl,

                                                           trade_exposure = est_parameters[["trade_exposure"]],
                                                           trade_exposure_period = est_parameters[["trade_exposure_period"]],

                                                           rfs_lead = est_parameters[["leads"]],
                                                           rfs_lag = est_parameters[["lags"]],
                                                           aggr_lead = est_parameters[["aggr_lead"]],
                                                           aggr_lag = est_parameters[["aggr_lag"]],
                                                           sjpos = est_parameters[["sjpos"]],
                                                           s_trend = est_parameters[["s_trend"]],
                                                           s_trend_loga = est_parameters[["s_trend_loga"]],
                                                           s_trend_sq = est_parameters[["s_trend_sq"]],

                                                           glm_iter = est_parameters[["glm_iter"]],

                                                           output = "coef_table")
  #}
  
}

# est_obj <- gentest_crops_est[[loss_type]][[CNT]][[CROP]][[1]]
# est_dat <- gentest_crops_est[[loss_type]][[CNT]][[CROP]][[2]]
# est_df <- gentest_crops_est[[loss_type]][[CNT]][[CROP]][[3]]
# unique(est_dat$year)
for(CNT in continents){
  sapply(all_tests_est[[CNT]], FUN = function(eo){eo$convStatus}) %>% print()
}


all_tests_est[["all"]][["loss_oilpalm_both"]] %>% post_est_fnc(base_exposure = "eaear_Oilpalm")
all_tests_est[["America"]][["eaear_Tobacco"]] %>% summary()
all_tests_est[["America"]][["eaear_Soy_compo"]] %>% summary()
all_tests_est[["Asia"]][["eaear_Oilpalm"]]$convStatus
all_tests_est[["Asia"]][["eaear_Oilpalm"]] %>% summary()

### PLOT COEFFICIENTS ### 
# prepare regression outputs in a tidy data frame readable by dwplot
# this works annual_controls is TRUE and we add (as controls) all crops that we are interested in. 
rm(df)
cnt_crop_list <- list()
cnt_list <- list()
for(CNT in continents){# "all", , "Africa", "Asia"
  for(CROP in c(gentest_crops, displtest_crops, "eaear_Oilpalm")){
    # extract t stats of coef of interest 
    x <- all_tests_est[[CNT]][[CROP]]
    x <- data.frame(tstat = x[grepl("_aggrall",row.names(x)),"t value"])
    
    x$model <- CNT
    x$via_trade <- grepl("_expo_", row.names(x))    
    x$crop <- CROP
    
    cnt_crop_list[[CNT]][[CROP]] <- x
  }
  cnt_list[[CNT]] <- cnt_crop_list[[CNT]] %>% bind_rows()
}
df <- cnt_list %>% bind_rows() # "loss_cropland"

df$significant01 <- ""
#df$significant01 <- abs(df[,"tstat"])
df[abs(df[,"tstat"]) > 1.645, "significant01"] <- "p-value < .1"

# attribute crop Group 
# horrible code mais j'ai pas réussi à faire mieux ... 
df_gn <- sapply(names(crops_groups), function(gn){if_else(df$crop %in% crops_groups[[gn]], true = gn, false = NULL)}) 
df$Group <- ""
for(i in 1:nrow(df_gn)){
  row <- df_gn[i,]
  df$Group[i] <- row[!is.na(row)] %>% unname()
}

# make some variable name changes necessary for dotwhisker
# df <- df[df$crop != "eaear_Maizegrain",]
df <- mutate(df, sizes = if_else(model == "all", true = 4, false = 2.5))
names(df)[names(df)=="tstat"] <- "estimate"
# names(df)[names(df)=="model"] <- "continent"
# names(df)[names(df)=="Dynamics"] <- "model"
names(df)[names(df)=="crop"] <- "term"
# df$cumulative <- factor(df$Dynamics=="Cumulative")
head(df)


dwplot(df,
       dot_args = list(aes(color = model, shape = model, size = model, alpha = significant01)) ) %>%  #size = 2.5,  shape = Dynamics, size = cumulative
  relabel_predictors(predictors_dict) + 
  # # critical values  : 1.645   1.960  2.576
  geom_vline(xintercept = 1.96, colour = "grey60", linetype = "dotdash", alpha = 1) +
  geom_vline(xintercept = -1.96, colour = "grey60", linetype = "dotdash", alpha = 1) +
  geom_vline(xintercept = 2.576, colour = "grey60", linetype = "dotted", alpha = 1) +
  geom_vline(xintercept = -2.576, colour = "grey60", linetype = "dotted", alpha = 1) +
  
  facet_wrap(facets = ~Group, scales = "free_y") + 
  
  scale_color_brewer(type = "qual",palette="Set1", #  Accent  
                     breaks=c("all", "America", "Africa", "Asia"),
                     labels = c("Pan-tropical", "America (tropical lat.)", "Africa (tropical lat.)", "Asia (tropical lat.)"),
                     name="") +
  
  scale_shape_manual(values = c(18, 17, 15, 16), # c(21, 24, 22, 25),# c(1, 2, 0, 3), # c(16, 17, 15, 18)
                     breaks=c("all", "America", "Africa", "Asia"),
                     labels = c("Pan-tropical", "America (tropical lat.)", "Africa (tropical lat.)", "Asia (tropical lat.)"),
                     name="") +
  
  scale_size_manual(values = c(3.5, 2.5, 2.5, 2.5),
                    breaks=c("all", "America", "Africa", "Asia"),
                    labels = c("Pan-tropical", "America (tropical lat.)", "Africa (tropical lat.)", "Asia (tropical lat.)"),
                    name="") +
  
  scale_x_continuous(breaks = c(-1.96, 1.96), 
                     labels = c("-1.96","1.96")) +
  
  theme_bw() + xlab("t-statistics") + ylab("") +  
  theme(plot.title = element_text(face="bold", size=c(10)),
        legend.position = "bottom",# c(0.8, 0.05),
        legend.justification = c(0, 0), 
        legend.background = element_rect(colour="grey80"),
        legend.title = element_blank()) 















































#### ALL TESTS - LOSS COMMODITY ------------------------------------------
main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "loss_commodity_aeaycompo_long_final.Rdata"))
# release some memory upfront
main_data <- dplyr::filter(main_data, year >= 2008, year <= 2019)

loss_type <- "loss_commodity"
# loss commodity is correctly measured up to 2019. 
est_parameters[["end_year"]] <- 2019
# and the controlled crops are all available, as the outcome captures potentially deforestation for any of them
crops_ctrl <- all_crops

all_tests_losscommo_est <- list()
for(CNT in continents){
  
  # joint estimation du coup 
#  for(CROP in c(gentest_crops, displtest_crops, "eaear_Citrus", "eaear_Oilpalm", "eaear_Coconut", g6_crops)){
    ## Determine the crops to control for 
   # crops_ctrl <- crops_groups[sapply(crops_groups, match, x = CROP, nomatch = FALSE)==0] %>% unlist() %>% unname()  # this removes from the set of controls, the crops that are in the same group as CROP
    
    # the order of the strata in the list are determined by how we want to plot results
    all_tests_losscommo_est[[CNT]] <- make_main_reg(continent = CNT,
                                                  outcome_variable = loss_type,
                                                  
                                                  start_year = est_parameters[["start_year"]],
                                                  end_year = est_parameters[["end_year"]],
                                                  
                                                  # regress on one crop at a time, and control for other crops
                                                  exposure_rfs = all_crops, # CROP - THIS DICTATES JOINT ESTIMATION 
                                                  control_all_absolute_rfs = FALSE, # ANBD THIS THAT PREVENTS ANY CROP TO BE IN THE CONTROL CATEGORIE
                                                  annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                  all_exposures_rfs = crops_ctrl,
                                                  
                                                  trade_exposure = est_parameters[["trade_exposure"]],
                                                  trade_exposure_period = est_parameters[["trade_exposure_period"]],
                                                  
                                                  rfs_lead = est_parameters[["leads"]],
                                                  rfs_lag = est_parameters[["lags"]],
                                                  aggr_lead = est_parameters[["aggr_lead"]],
                                                  aggr_lag = est_parameters[["aggr_lag"]],
                                                  sjpos = est_parameters[["sjpos"]],
                                                  s_trend = est_parameters[["s_trend"]],
                                                  s_trend_loga = est_parameters[["s_trend_loga"]],
                                                  s_trend_sq = est_parameters[["s_trend_sq"]],
                                                  
                                                  output = est_parameters[["output"]])
    
  #}
}

# est_obj <- gentest_crops_est[[loss_type]][[CNT]][[CROP]][[1]]
# est_dat <- gentest_crops_est[[loss_type]][[CNT]][[CROP]][[2]]
# est_df <- gentest_crops_est[[loss_type]][[CNT]][[CROP]][[3]]
# unique(est_dat$year)


all_tests_losscommo_est[["America"]]%>% summary()

### PLOT COEFFICIENTS ### 
# prepare regression outputs in a tidy data frame readable by dwplot
rm(df)
cnt_list <- list()
for(CNT in continents){# "all", , "Africa", "Asia"
    # extract t stats of coef of interest 
    x <- post_est_fnc(est_obj = all_tests_losscommo_est[[CNT]], 
                      base_exposure = all_crops, # NOTE THIS FOR POST JOINT ESTIMATION
                      aggr_lag = est_parameters[["aggr_lag"]], 
                      aggr_lead = est_parameters[["aggr_lead"]], 
                      output = "tstat")
    
    x$model <- CNT
    x$via_trade <- grepl("_expo_", row.names(x))    
    x$crop <- gsub(pattern = "_X_.*$", x = row.names(x), replacement = "") # replace everything after and including _X_ with nothing

    x$Dynamics <- ""
    x$Dynamics[grepl("aggrall", row.names(x))] <- "Cumulative"
    x$Dynamics[grepl("lead1", row.names(x))] <- "t+1"
    x$Dynamics[row.names(x)%in%paste0(all_crops,"_X_statute_conv")] <- "t"
    x$Dynamics[grepl("lag1", row.names(x))] <- "t-1"
    x$Dynamics[grepl("lag2", row.names(x))] <- "t-2"
    x$Dynamics[grepl("lag3", row.names(x))] <- "t-3"
    
    # order dynamics 
    # row.names(x) <- x$Dynamics  
    # x <- x[c("Cumulative", "t+1", "t", "t-1", "t-2", "t-3", "t-4"),] 
    # # remove NA rows generated if some dynamics are absent for this crop 
    # x <- x[!is.na(x$tstat), ]
    row.names(x) <- NULL # otherwise they duplicate while stacking
    cnt_list[[CNT]] <- x
    
}
df <- cnt_list %>% bind_rows() # "loss_cropland"

df[df$model == "all", "model"] <- "Pan-tropical"
df[df$model == "America", "model"] <- "America (tropical lat.)"
df[df$model == "Africa", "model"] <- "Africa (tropical lat.)"
df[df$model == "Asia", "model"] <- "Asia (tropical lat.)"
df$significant01 <- ""
#df$significant01 <- abs(df[,"tstat"])
df[abs(df[,"tstat"]) > 1.645, "significant01"] <- "p-value < .1"

# attribute crop Group 
# horrible code mais j'ai pas réussi à faire mieux ... 
df_gn <- sapply(names(crops_groups), function(gn){if_else(df$crop %in% crops_groups[[gn]], true = gn, false = NULL)}) 
df$Group <- ""
for(i in 1:nrow(df_gn)){
  row <- df_gn[i,]
  df$Group[i] <- row[!is.na(row)] %>% unname()
}

# relabel terms (for y axis)
# row.names(df) <- NULL
# df <- mutate(df, term = paste0(crop,".",Dynamics))
# predictors_dict <- rep("", length(df$term)) 
# names(predictors_dict) <- df$term
# predictors_dict[grepl("Cumulative", names(predictors_dict))] <- df$crop[grepl("Cumulative", names(predictors_dict))]
# df <- relabel_predictors(df,predictors_dict) 
predictors_dict <- c(eaear_Maizegrain = "Maize",
                     eaear_Fodder = "Pasture",
                     eaear_Cereals = "Cereals",
                     eaear_Rice = "Rice",
                     eaear_Soy_compo = "Soy",
                     eaear_Cotton = "Cotton",
                     eaear_Oilfeed_crops = "Other oil crops",
                     eaear_Biomass = "Biomass crops",
                     eaear_Sugarcane = "Sugarcane",
                     eaear_Tobacco = "Tobacco",
                     eaear_Citrus = "Citrus",
                     eaear_Oilpalm = "Oil palm",
                     eaear_Coconut = "Coconut",
                     eaear_Banana = "Banana",
                     eaear_Cocoa_Coffee = "Cocoa or Coffee",
                     eaear_Tea = "Tea",
                     eaear_Rubber = "Rubber" )

# eaear_Maizegrain.Cumulative = "Maize",
# `eaear_Maizegrain.t+1` = "",
# `eaear_Maizegrain.t` = "",
# `eaear_Maizegrain.t-1` = "",
# `eaear_Maizegrain.t-2` = "",
# `eaear_Maizegrain.t-3` = "",

# df <- mutate(df, y_mapping = paste0(term, "-", Dynamics))

head(df)
# make some variable name changes necessary for dotwhisker
df_save <- df

# df <- df[df$crop != "eaear_Maizegrain",]
df <- df[df$Dynamics=="Cumulative",]
names(df)[names(df)=="tstat"] <- "estimate"
# names(df)[names(df)=="model"] <- "continent"
# names(df)[names(df)=="Dynamics"] <- "model"
names(df)[names(df)=="crop"] <- "term"
# df$cumulative <- factor(df$Dynamics=="Cumulative")

dwplot(df,
       dot_args = list(size = 2.5, aes(color = model, shape = model, alpha = significant01)) ) %>%  # shape = Dynamics, size = cumulative
  relabel_predictors(predictors_dict) + 
  # # critical values  : 1.645   1.960  2.576
  # geom_vline(xintercept = 1.645, colour = "grey60", linetype = "solid", alpha = 0.5) +
  # geom_vline(xintercept = -1.645, colour = "grey60", linetype = "solid", alpha = 0.5) +
  geom_vline(xintercept = 1.96, colour = "grey60", linetype = "dotdash", alpha = 1) +
  geom_vline(xintercept = -1.96, colour = "grey60", linetype = "dotdash", alpha = 1) +
  geom_vline(xintercept = 2.576, colour = "grey60", linetype = "dotted", alpha = 1) +
  geom_vline(xintercept = -2.576, colour = "grey60", linetype = "dotted", alpha = 1) +
  
  facet_wrap(facets = ~Group, scales = "free_y") + 
  
  scale_color_brewer(type = "qual", palette="Set1") +
  
  scale_x_continuous(breaks = c(-1.96, 1.96), 
                     labels = c("-1.96","1.96")) +
  
  theme_bw() + xlab("t-statistics") + ylab("") +  
  theme(plot.title = element_text(face="bold", size=c(10)),
        legend.position = "bottom",# c(0.8, 0.05),
        legend.justification = c(0, 0), 
        legend.background = element_rect(colour="grey80"),
        legend.title = element_blank()) 
# scale_size_manual(breaks = c(TRUE, FALSE), 
#                   values = c(3, 1.5), 
#                   labels = c("Cumulative", "Annual")) +
















































#### DES STATS RFS --------------------------------------------------------------------------------------------------
# add PSD data to the chart 
psd <- readRDS(here("temp_data", "prepared_psd.Rdata"))
# UNITS
# from https://apps.fas.usda.gov/psdonline/app/index.html#/app/about#G5
# Production, Trade, & Use: 1000 metric tons.
# So once divided by 1000, it's expressed in million tons 
psd$us_dc_maize <- psd$UnitedStates.Domestic_Consumption.Maize / 1000
psd$us_ts_maize <- psd$UnitedStates.Total_Supply.Maize / 1000
psd$us_im_maize <- psd$UnitedStates.Imports.Maize / 1000

#head(psd)
w_rfs <- left_join(rfs, psd[,c("year", "us_dc_maize")], by = "year")

# w_rfs <- dplyr::mutate(w_rfs, us_ex_maize = us_ts_maize - us_dc_maize)

# give back RFS1 values, to display them
w_rfs[w_rfs$year>=2006 & w_rfs$year<=2007, ] <- w_rfs %>% 
  dplyr::filter(year>=2006 & year<=2007) %>% 
  dplyr::mutate(statute_total = c(0, 0), # give zeros for total, because this is going to represent RFS2
                final_total = c(0, 0), 
                statute_advanced = c(0, 0), 
                final_advanced = c(0, 0))
# reconstitute conventional
w_rfs[w_rfs$year>=2006 & w_rfs$year<=2007, ] <- w_rfs %>% 
  dplyr::filter(year>=2006 & year<=2007) %>% 
  dplyr::mutate(statute_conv = statute_total - statute_advanced, 
                final_conv = final_total - final_advanced)

# continue final_conv
w_rfs[is.na(w_rfs$final_conv), "final_conv"] <- 15
w_rfs[is.na(w_rfs$us_dc_maize), "us_dc_maize"] <- 309.132

# group biodiesel and (other) advanced, not relevant to show distinctively
w_rfs <- w_rfs %>% 
  dplyr::mutate(statute_alladvanced = statute_advanced + statute_biodiesel, 
                final_alladvanced = final_advanced + final_biodiesel)


# take what's needed
w_rfs <- w_rfs %>% 
  dplyr::filter(year>=2000) %>% 
  dplyr::select(year, statute_conv, final_conv, us_dc_maize) #  statute_alladvanced, final_alladvanced,

# add RFS1 (2000-2012 only)
w_rfs$statute_rfs1 <- c(0, 0, 0, 0, 0, 0, 4, 4.7, 5.4, 6.1, 6.8, 7.4, 7.5, rep(NA, 10) ) 

# convert mandates from billion gallons to tons maize 
# According to USDA, 1 bushel maize gives 2.7 gallons ethanol https://www.ers.usda.gov/about-ers/partnerships/strengthening-statistics-through-the-icars/biofuels-data-sources/
# Given that 1 bushel maize is 0.0254 metric ton maize https://www.sagis.org.za/conversion_table.html
# 1 gallon ethanol = (1/2.7) * 0.0254 metric ton
# 1bgal = (1/2.7) * 0.0254 billion metric tons = (1/2.7) * 0.0254 * 1000 million metric tons
w_rfs <- dplyr::mutate(w_rfs, across(.cols = contains("_conv") | contains("_rfs1"), .fns = ~.*(1/2.7) * 0.0254 * 1000) )

# pile up 
l_rfs <- pivot_longer(w_rfs, cols = c("statute_conv", "final_conv", "statute_rfs1"), names_to = "mandates", values_to = "maize_milton") %>% as.data.frame() %>% arrange(mandates)


ggplot(l_rfs, aes(x = year, y = maize_milton, group = mandates)) +
  geom_line(aes(linetype=mandates)) + 
  geom_line(aes(x = year, y = us_dc_maize)) +
  geom_label(aes(label="US domestic consumption", 
                 x=2019,
                 y=290)) +
  
  geom_label(aes(label="RFS2 (statutory)", 
                 x=2020,
                 y=150)) +
  
  geom_label(aes(label="RFS2 (final rule)", 
                 x=2016,
                 y=120)) +
  
  geom_label(aes(label="RFS1 (statutory)", 
                 x=2014,
                 y=70)) +
  scale_linetype_manual(breaks=c("statute_conv", "final_conv", "statute_rfs1"),
                        values=c("solid", "dotted", "twodash"),
                        labels=c("RFS2 (statutory)", "RFS2 (final)", "RFS1 (statutory)"),
                        name="Mandates") +
  
  scale_x_continuous(breaks = c(2000, 2005, 2006, 2008, 2012, 2015, 2022)) +
  scale_y_continuous(name = "Million tonnes maize") + 
  theme_minimal() +
  theme(legend.position="none", 
        plot.title = element_text(size = 10, face = "bold"), 
        axis.title.y.left=element_text(size=10,face="bold", hjust = 1),
        axis.title.y.right=element_text(size=10,face="bold", hjust = 1),
        axis.title.x=element_blank(), #element_text(size=10,face="bold", hjust = 0.5), 
        panel.grid = element_line(inherit.blank = TRUE))  









