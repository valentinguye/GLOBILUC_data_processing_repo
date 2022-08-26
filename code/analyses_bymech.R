

### GLOBAL CRS USED ### 
mercator_world_crs <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "

### OBJECTS USED IN RFS PROCESSES ###
main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "loss_aeaycompo_long_final.Rdata"))
main_data <- dplyr::mutate(main_data, loss_commodity = loss_cropland + loss_oilpalm_both + loss_pasture)

# loss_commo <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "loss_commodity_aeaycompo_long_final.Rdata"))
# summary(loss_commo$loss_commodity)
# summary(main_data$loss_commodity)
# rm(loss_commo)

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "loss_commodity_aeaycompo_long_final.Rdata"))

# release some memory upfront
main_data <- dplyr::filter(main_data, year >= 2008, year <= 2019)

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
exposure_rfs <- "eaear_Cereals"
control_all_absolute_rfs <- TRUE # whether to control for other crops
annual_rfs_controls <- TRUE # on annual or year averaged dynamics
all_exposures_rfs = cash_crops_ctrl # which crops

trade_exposure = "trade_expo_imp"
trade_exposure_period = "20012007"

# treatment dynamics
original_rfs_treatments <- c("statute_conv")
# These 4 arg. below determine the mandates that are interacted with the exposure(s) of the crop(s) of interest.
# Not the dynamics wpecified for control crops (this is done by annual_rfs_controls), 
# nor the dynamic effects that are to be aggregated eventually. 
rfs_lead = 0
rfs_lag = 0
rfs_fya = 0 
rfs_pya = 0
# These, below, determine the dynamic effects to aggregate
aggr_lead = rfs_lead
aggr_lag = rfs_lag

# those are always FALSE 
group_exposure_rfs <- FALSE
most_correlated_only = FALSE
exposure_pasture <- FALSE 
sjpos <- FALSE # should the sample be restricted to cells where sj is positive? 

# control heterogeneity
control_pasture <- FALSE
pasture_trend <- FALSE
fc_trend <- FALSE
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
continent = CNT 
outcome_variable = loss_type 
start_year = est_parameters[["start_year"]]
end_year = est_parameters[["end_year"]]
exposure_rfs = CROP
control_all_absolute_rfs = TRUE
annual_rfs_controls = TRUE
all_exposures_rfs = cash_crops_ctrl
trade_exposure = est_parameters[["trade_exposure"]]
trade_exposure_period = est_parameters[["trade_exposure_period"]]
rfs_lead = est_parameters[["leads"]] 
rfs_lag = est_parameters[["lags"]]
rfs_fya = est_parameters[["fya"]] 
rfs_pya = est_parameters[["pya"]] 
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
                          control_pasture = FALSE,
                          pasture_trend = FALSE,
                          remaining = FALSE, # should remaining forest be controlled for STOP DOING THIS BECAUSE IT INTRODUCES NICKELL BIAS
                          s_trend = FALSE,
                          s_trend_loga = FALSE,
                          s_trend_sq = FALSE,
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
    
    all_eaear_trade_exposure_rfs <- c()
    if(length(trade_exposure)>0){
      trade_expo_dat <- readRDS(here("temp_data", "processed_trade_exposures", paste0("trade_exposures_", trade_exposure_period,".Rdata")))
      # if trade_exposure is not _imp, the following line will keep both expo and expo_imp variables, but it already divides by 2 the number of cols that are joined to the bigger d dataframe
      trade_expo_dat <- trade_expo_dat[, (grepl(x = names(trade_expo_dat), pattern = trade_exposure) | names(trade_expo_dat)=="country_name") ]
      
      d <- left_join(d, trade_expo_dat, by = "country_name")
      
      # make ALL the eaear-trade exposures
      names(d)
      for(exp_rfs in unique(c(exposure_rfs, all_exposures_rfs))){
        # identify the corresponding trade_expo 
        trade_expo_j <- names(d)[names(d) == paste0(trade_exposure,"_",str_to_lower(gsub(pattern = "eaear_", replacement = "", x = exp_rfs) ) ) ]
        eaear_trade_expo_j <- paste0("eaear_",trade_expo_j)
        d <- mutate(d, !!as.symbol(eaear_trade_expo_j) := !!as.symbol(exp_rfs) * !!as.symbol(trade_expo_j))
        all_eaear_trade_exposure_rfs <- c(all_eaear_trade_exposure_rfs, eaear_trade_expo_j)
      }
      rm(trade_expo_dat)
    }
    
    # code below should work whether d is from tmf or losscommo 
    # Keep only in data the useful variables 
    d <- dplyr::select(d, all_of(unique(c("grid_id", "year", "year_ln", "year_sq", "lat", "lon", "continent_name", "country_name", "country_year",
                                          # "fc_2000", "fc_2009", "remaining_fc", # "accu_defo_since2k",
                                          "grid_id_5", "grid_id_10", "grid_id_20", "grid_id_5_year", "grid_id_10_year", "grid_id_20_year",
                                          outcome_variable,# "tmf_agri", "tmf_flood", "tmf_plantation",
                                          "pasture_share_2000",
                                          exposure_rfs, all_eaear_trade_exposure_rfs, all_exposures_rfs )))) #sj, 
    
    # Merge only the prices needed, not the whole price dataframe
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
    for(exp_rfs in exposure_rfs){
      # identify the corresponding trade_expo 
      trade_expo_j <- paste0("eaear_",trade_exposure,"_",str_to_lower(gsub(pattern = "eaear_", replacement = "", x = exp_rfs) ) ) 
      eaear_trade_exposure_rfs <- c(eaear_trade_exposure_rfs, trade_expo_j)
    }
  }
  
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
      for(exp_rfs in c(exposure_rfs, eaear_trade_exposure_rfs)){
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
  
  
  # from here, if we are in rfs process, sj may refer to the group of exposures
  if(group_exposure_rfs){
    sj <- "grp_exp"
  }
  

  # VARYING SLOPES makes computations much longer (and may not converge)
  # We specify heterogeneous trends using fixest. 
  # as inspired by the answer from Laurent Bergé here https://stackoverflow.com/questions/34232834/fixed-effects-regression-with-state-specific-trends
  # and by the vignettes https://cran.r-project.org/web/packages/fixest/vignettes/fixest_walkthrough.html#41_Interactions_involving_fixed-effects
  # and documentation of fixest https://lrberge.github.io/fixest/reference/feglm.html#varying-slopes-1
  # the point is: we want to allow the year variable to have a different slope (coefficient) for every level of exposure variable
  
  # if(s_trend){
  #   for(eaear_exp_rfs in exposure_rfs){
  #     fe <- paste0(fe, paste0(" + ",eaear_exp_rfs,"[[year]]")) 
  #   }
  # }
  # 
  # if(s_trend_loga){
  #   for(eaear_exp_rfs in exposure_rfs){
  #     fe <- paste0(fe, paste0(" + ",eaear_exp_rfs,"[[year_ln]]")) 
  #   }
  # }
  # if(s_trend_sq){
  #   for(eaear_exp_rfs in exposure_rfs){
  #     fe <- paste0(fe, paste0(" + ",eaear_exp_rfs,"[[year_sq]]")) 
  #   }
  # }
  
  
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
  if(s_trend){
    # Suitability time trend
    # d <- mutate(d, suitability_trend := !!as.symbol(sj)*(year))#-2000 # doesn't change anything that it's multiplied by 1...19 or 2001...2019
    # controls <- c(controls, "suitability_trend")
    for(eaear_exp_rfs in exposure_rfs){
      varname <- paste0(eaear_exp_rfs, "_trend")
      controls <- c(controls, varname)
      d <- mutate(d, !!as.symbol(varname) := !!as.symbol(eaear_exp_rfs) * year) # for the linear trend it does not matter
    }
  }
  if(s_trend_loga){
    for(eaear_exp_rfs in exposure_rfs){
      varname <- paste0(eaear_exp_rfs, "_trend_loga")
      controls <- c(controls, varname)
      #make the logarithmic trend start when RFS actually starts, such that it fits best
      # so it should be log(1) in 2005 (although never occurring in the data used for analysis but the log trend starts from there)
      d <- mutate(d, !!as.symbol(varname) := !!as.symbol(eaear_exp_rfs) * year_ln)
    }
  }
  if(s_trend_sq){
    for(eaear_exp_rfs in exposure_rfs){
      varname <- paste0(eaear_exp_rfs, "_trend_sq")
      controls <- c(controls, varname)
      d <- mutate(d, !!as.symbol(varname) := !!as.symbol(eaear_exp_rfs) * year_sq)
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
                           #fixef.iter = 100000,
                           nthreads = 3,
                           fixef.rm = "perfect",
                           glm.iter = glm_iter,
                           notes = TRUE, 
                           verbose = 1) 
  
  df_res <- summary(reg_res)$coeftable#[paste0(original_sj, "_X_", original_Pk), ]
  
  fixest_df <- degrees_freedom(reg_res, type = "t")
  
  ## MAKE AGGREGATE RESULTS 
  
  ## Annual mandates
  # joint OR DISJOINT estimation
  if((rfs_lead > 0 | rfs_lag > 0) & rfs_fya == 0 & rfs_pya == 0 ){
    for(EOI in c(exposure_rfs, eaear_trade_exposure_rfs)){
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
  if(output == "everything"){
    toreturn <- list(reg_res, d_clean, df_res) 
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
  
  
  rm(d_clean, df_res)
  return(toreturn)
  rm(toreturn)
}

#### Main specification -----------------------------------------------------------------------
est_parameters <- list(start_year = 2008, 
                       end_year = 2022, 
                       trade_exposure = NULL, # "trade_expo", # "trade_expo_imp", #  "trade_expo_imp", # NULL, # "trade_expo_imp",
                       trade_exposure_period = "20012007",
                       annual_rfs_controls = TRUE,
                       lags = 3,
                       leads = 1,
                       fya = 0, 
                       pya = 0,
                       aggr_lag = 1, 
                       aggr_lead = 1,
                       s_trend = FALSE, 
                       s_trend_loga = FALSE,
                       s_trend_sq = FALSE,
                       clustering = "oneway",
                       cluster_var1 = "grid_id_10")

crops <- list(maize = "eaear_Maizegrain",
              marg_land = c("eaear_Fodder", "eaear_Biomass"), # this only serves as a control 
              group1 = c("eaear_Cereals", "eaear_Rice"), # "eaear_Fodder", 
              group2 = c("eaear_Soy_compo", "eaear_Cotton",  "eaear_Oilfeed_crops"), 
              group4 = c("eaear_Sugarcane"), # "eaear_Biomass", 
              group5 = c("eaear_Tobacco", "eaear_Citrus"), 
              group3 = c("eaear_Coconut", "eaear_Oilpalm"), 
              group6 = c("eaear_Banana", "eaear_Cocoa_Coffee", "eaear_Tea", "eaear_Rubber"))

# this is used to remove some controls in loss_cropland regressions
woody_perrenials <- c("eaear_Citrus", "eaear_Coconut", "eaear_Oilpalm", 
                      "eaear_Banana", "eaear_Cocoa_Coffee", "eaear_Rubber")


crops_outcomes <- c("loss_cropland")# "loss_cropland", 

#### SUBSTITUTION OR DISPLACEMENT TESTS  ------------------------------------------
# All crops, by common set of mechanisms predicted to apply on them

crops_pooled12 <- list(maize = "eaear_Maizegrain",
                        group1 = c("eaear_Fodder", "eaear_Cereals", "eaear_Rice", "eaear_Soy_compo", "eaear_Cotton",  "eaear_Oilfeed_crops"), 
                        group4 = c("eaear_Biomass", "eaear_Sugarcane"), 
                        group5 = c("eaear_Tobacco", "eaear_Citrus"),
                        group3 = c("eaear_Coconut", "eaear_Oilpalm"), 
                        group6 = c("eaear_Banana", "eaear_Cocoa_Coffee", "eaear_Tea", "eaear_Rubber"))

crop_groups <- crops 
# crop_groups <- crops_pooled12


# The set of crops which share the same dynamics 
gentest_crops <- c("eaear_Maizegrain", 
                "eaear_Cereals", "eaear_Rice",
                "eaear_Soy_compo", "eaear_Cotton",  "eaear_Oilfeed_crops")


gentest_crops_est <- list()

loss_type <- "loss_commodity"
CNT <- "America"
CROP <- "eaear_Soy_compo"

for(loss_type in crops_outcomes){
  for(CNT in c("America", "Africa", "Asia")){# "all", 
    for(CROP in gentest_crops){
      ## Determine the crops to control for 
      # this removes from the set of controls, the crops that are in the same group as CROP
      crops_ctrl <- crop_groups[sapply(crop_groups, match, x = CROP, nomatch = FALSE)==0] %>% unlist() %>% unname()
      # crops_ctrl <- unlist(unname(crop_groups))
      if(loss_type =="loss_cropland"){
        # further remove woody perennials and pastures, if we regress cropland anyway
       crops_ctrl <- crops_ctrl[sapply(crops_ctrl, function(c){!(c %in% c(woody_perrenials))})]# , "eaear_Fodder"
       # and loss_cropland is correctly measured only up to 2016 
       est_parameters[["end_year"]] <- 2016
      }
      # the order of the strata in the list are determiend by how we want to plot results
      gentest_crops_est[[loss_type]][[CNT]][[CROP]] <- make_main_reg(continent = CNT,
                                                     outcome_variable = loss_type,

                                                     start_year = est_parameters[["start_year"]],
                                                     end_year = est_parameters[["end_year"]],

                                                     # regress on one crop at a time, and control for other crops
                                                     exposure_rfs = CROP,
                                                     control_all_absolute_rfs = TRUE,
                                                     annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                     all_exposures_rfs = crops_ctrl,

                                                     trade_exposure = est_parameters[["trade_exposure"]],
                                                     trade_exposure_period = est_parameters[["trade_exposure_period"]],

                                                     rfs_lead = est_parameters[["leads"]],
                                                     rfs_lag = est_parameters[["lags"]],
                                                     rfs_fya = est_parameters[["fya"]],
                                                     rfs_pya = est_parameters[["pya"]],
                                                     aggr_lead = est_parameters[["aggr_lead"]],
                                                     aggr_lag = est_parameters[["aggr_lag"]],
                                                     s_trend = est_parameters[["s_trend"]],
                                                     s_trend_loga = est_parameters[["s_trend_loga"]],
                                                     s_trend_sq = est_parameters[["s_trend_sq"]],

                                                     output = "coef_table"  )
      
      # gentest_crops_est[[loss_type]][[CNT]][[CROP]] <- gentest_crops_est[[loss_type]][[CNT]][[CROP]][[3]]
    }
  }    
}

# est_obj <- gentest_crops_est[[loss_type]][[CNT]][[CROP]][[1]]
# est_dat <- gentest_crops_est[[loss_type]][[CNT]][[CROP]][[2]]
# est_df <- gentest_crops_est[[loss_type]][[CNT]][[CROP]][[3]]
# unique(est_dat$year)
### PLOT COEFFICIENTS ### 
# prepare regression outputs in a tidy data frame readable by dwplot
rm(df)
loss_cnt_crop_list <- list()
loss_cnt_list <- list()
loss_list <- list()
for(loss_type in crops_outcomes){
  for(CNT in c("America", "Africa", "Asia")){# "all", , "Africa", "Asia"
    for(CROP in gentest_crops){
      x <- gentest_crops_est[[loss_type]][[CNT]][[CROP]]%>% as.data.frame()
      x$term <- CROP
      x$model <- CNT
      x$loss_type <- loss_type
      x <- x[grepl("_X_aggrall", row.names(x)), ]
      x$via_trade <- grepl("_expo_", row.names(x))
      loss_cnt_crop_list[[loss_type]][[CNT]][[CROP]] <- x
    }
    loss_cnt_list[[loss_type]][[CNT]] <- loss_cnt_crop_list[[loss_type]][[CNT]] %>% bind_rows()
  }
  loss_list[[loss_type]] <- loss_cnt_list[[loss_type]] %>% bind_rows()
}
# Choose one outcome type here 
df <- loss_list %>% bind_rows() # "loss_cropland"
df[df$model == "all", "model"] <- "Global"
#df$term <- rep(eaear_mapmat[,"Crops"], length(dyn_df_list))
df$significant01 <- ""
df[df[,"Pr(>|t|)"] < 0.1, "significant01"] <- "p-value < .1"

names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"
# decide to remove or not sugarcane 
# x <- x[!grepl(pattern = "expo_X_sugarcane", x = x$term),]
row.names(df) <- NULL
df

{dwplot(df,
        dot_args = list(size = 2.5,aes(shape = model, alpha = significant01)),#, alpha = significant01
        whisker_args = list(size = 1, aes(alpha = significant01)),#alpha = significant01
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(
      eaear_Maizegrain = "Maize",
      eaear_Fodder = "Pasture",
      eaear_Cereals = "Cereals",
      eaear_Rice = "Rice",
      eaear_Soy_compo = "Soy", 
      eaear_Cotton = "Cotton",
      eaear_Oilfeed_crops = "Other oil crops"
    )) +
    #scale_color_manual(values=wes_palette(n=4, name="Darjeeling1", type = "discrete")) +
    scale_color_brewer(type = "qual", palette="Set1") +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    #ggtitle(title) +
    theme(plot.title = element_text(face="bold", size=c(10)),
          legend.position = c(0.8, 0.05),
          legend.justification = c(0, 0), 
          legend.background = element_rect(colour="grey80"),
          legend.title = element_blank()) }



#### DISPLACEMENT TESTS -------------------------------------------------------------
crops_pooled45 <- list(maize = "eaear_Maizegrain",
                       group1 = c("eaear_Fodder", "eaear_Cereals", "eaear_Rice"), 
                       group2 = c("eaear_Soy_compo", "eaear_Cotton",  "eaear_Oilfeed_crops"), 
                       group4 = c("eaear_Biomass", "eaear_Sugarcane", "eaear_Tobacco", "eaear_Citrus"),
                       group3 = c("eaear_Coconut", "eaear_Oilpalm"), 
                       group6 = c("eaear_Banana", "eaear_Cocoa_Coffee", "eaear_Tea", "eaear_Rubber"))

crop_groups <- crops 
# crop_groups <- crops_pooled45

# The set of crops which share the same dynamics 
displtest_crops <- c("eaear_Biomass", "eaear_Sugarcane", "eaear_Tobacco", "eaear_Citrus")

displtestcrops_est <- list()

for(loss_type in crops_outcomes){
  for(CNT in c("America", "Africa", "Asia")){# "all", 
    for(CROP in displtest_crops){
      # determine the crops to control for 
      crops_ctrl <- crop_groups[sapply(crop_groups, match, x = CROP, nomatch = FALSE)==0] %>% unlist() %>% unname()
      if(loss_type =="loss_cropland"){
        crops_ctrl <- crops_ctrl[sapply(crops_ctrl, function(c){!(c %in% c(woody_perrenials, "eaear_Fodder"))})]
        est_parameters[["end_year"]] <- 2016
      }
      # the order of the strata in the list are determiend by how we want to plot results
      displtestcrops_est[[loss_type]][[CNT]][[CROP]] <- make_main_reg(continent = CNT,
                                                                  outcome_variable = loss_type,
                                                                  
                                                                  start_year = est_parameters[["start_year"]],
                                                                  end_year = est_parameters[["end_year"]],
                                                                  
                                                                  # regress on one crop at a time, and control for other crops
                                                                  exposure_rfs = CROP,
                                                                  control_all_absolute_rfs = TRUE,
                                                                  annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                                  all_exposures_rfs = crops_ctrl,
                                                                  
                                                                  trade_exposure = est_parameters[["trade_exposure"]],
                                                                  trade_exposure_period = est_parameters[["trade_exposure_period"]],
                                                                  
                                                                  rfs_lead = est_parameters[["leads"]],
                                                                  rfs_lag = est_parameters[["lags"]],
                                                                  rfs_fya = est_parameters[["fya"]],
                                                                  rfs_pya = est_parameters[["pya"]],
                                                                  aggr_lead = est_parameters[["aggr_lead"]],
                                                                  aggr_lag = est_parameters[["aggr_lag"]],
                                                                  s_trend = est_parameters[["s_trend"]],
                                                                  s_trend_loga = est_parameters[["s_trend_loga"]],
                                                                  s_trend_sq = est_parameters[["s_trend_sq"]],
                                                                  
                                                                  output = "coef_table"  )
      
      # displtestcrops_est[[loss_type]][[CNT]][[CROP]] <- displtestcrops_est[[loss_type]][[CNT]][[CROP]][[3]]
    }
  }    
}

est_obj <- displtestcrops_est[[loss_type]][[CNT]][[CROP]][[1]]
est_dat <- displtestcrops_est[[loss_type]][[CNT]][[CROP]][[2]]
est_df <- displtestcrops_est[[loss_type]][[CNT]][[CROP]][[3]]
unique(est_dat$year)
### PLOT COEFFICIENTS ### 
# prepare regression outputs in a tidy data frame readable by dwplot
rm(df)
loss_cnt_crop_list <- list()
loss_cnt_list <- list()
loss_list <- list()
for(loss_type in crops_outcomes){
  for(CNT in c("America", "Africa", "Asia")){# "all", , "Africa", "Asia"
    for(CROP in displtest_crops){
      x <- displtestcrops_est[[loss_type]][[CNT]][[CROP]]%>% as.data.frame()
      x$term <- CROP
      x$model <- CNT
      x$loss_type <- loss_type
      x <- x[grepl("_X_aggrall", row.names(x)), ]
      x$via_trade <- grepl("_expo_", row.names(x))
      loss_cnt_crop_list[[loss_type]][[CNT]][[CROP]] <- x
    }
    loss_cnt_list[[loss_type]][[CNT]] <- loss_cnt_crop_list[[loss_type]][[CNT]] %>% bind_rows()
  }
  loss_list[[loss_type]] <- loss_cnt_list[[loss_type]] %>% bind_rows()
}
# Choose one outcome type here 
df <- loss_list %>% bind_rows() # "loss_cropland"
df[df$model == "all", "model"] <- "Global"
#df$term <- rep(eaear_mapmat[,"Crops"], length(dyn_df_list))
df$significant01 <- ""
df[df[,"Pr(>|t|)"] < 0.1, "significant01"] <- "p-value < .1"

names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"
# decide to remove or not sugarcane 
# x <- x[!grepl(pattern = "expo_X_sugarcane", x = x$term),]
row.names(df) <- NULL
df

{dwplot(df,
        dot_args = list(size = 2.5,aes(shape = model, alpha = significant01)),#, alpha = significant01
        whisker_args = list(size = 1, aes(alpha = significant01)),#alpha = significant01
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(
      eaear_Biomass = "Biomass crops",
      eaear_Sugarcane = "Sugarcane", 
      
      eaear_Tobacco = "Tobacco", 
      eaear_Citrus = "Citrus"
    )) +
    #scale_color_manual(values=wes_palette(n=4, name="Darjeeling1", type = "discrete")) +
    scale_color_brewer(type = "qual", palette="Set1") +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    #ggtitle(title) +
    theme(plot.title = element_text(face="bold", size=c(10)),
          legend.position = c(0.8, 0.05),
          legend.justification = c(0, 0), 
          legend.background = element_rect(colour="grey80"),
          legend.title = element_blank()) }




#### OIL PALM #### 
oilpalm_outcomes <- c("loss_oilpalm_both", "loss_oilpalm_indus")# "loss_commodity", 
oilpalm_est <- list()

for(loss_type in oilpalm_outcomes){
  for(CNT in c("America", "Africa", "Asia")){# "all", "America", "Africa", 
    # determine the crops to control for 
    crops_ctrl <- crops[sapply(crops, match, x = "eaear_Oilpalm", nomatch = FALSE)==0] %>% unlist() %>% unname()
    if(grepl(loss_type, pattern = "oilpalm")){
      crops_ctrl <- c()
      est_parameters[["end_year"]] <- 2016
    }
    # the order of the strata in the list are determiend by how we want to plot results
    oilpalm_est[[loss_type]][[CNT]][["eaear_Oilpalm"]] <- make_main_reg(continent = CNT, 
                                                                outcome_variable = loss_type, 
                                                                start_year = 2011,
                                                                end_year = 2019,
                                                                
                                                                # regress on crop at a time, but control for annual interactions with other crops
                                                                exposure_rfs = "eaear_Oilpalm",
                                                                control_all_absolute_rfs = TRUE,
                                                                annual_rfs_controls = TRUE,
                                                                all_exposures_rfs = crops_ctrl,
                                                                
                                                                trade_exposure = NULL,
                                                                trade_exposure_period = "20012007",
                                                                
                                                                rfs_lead = 1,
                                                                rfs_lag = 3,
                                                                rfs_fya = 0,
                                                                rfs_pya = 0,
                                                                aggr_lead = 1,
                                                                aggr_lag = 3,
                                                                s_trend = FALSE,
                                                                s_trend_loga = FALSE,
                                                                s_trend_sq = FALSE,

                                                                output = "coef_table"  )
  }
}

rm(df)
loss_cnt_crop_list <- list()
loss_cnt_list <- list()
loss_list <- list()
for(loss_type in oilpalm_outcomes){
  for(CNT in c("America", "Africa", "Asia")){# "all", , "Africa", "Asia" , 
    x <- oilpalm_est[[loss_type]][[CNT]][["eaear_Oilpalm"]]%>% as.data.frame()
    x$term <- "eaear_Oilpalm"
    x$model <- CNT
    x$loss_type <- loss_type
    x <- x[grepl("_X_aggrall", row.names(x)), ]
    x$via_trade <- grepl("_expo_", row.names(x))
    loss_cnt_crop_list[[loss_type]][[CNT]][["eaear_Oilpalm"]] <- x
  
    loss_cnt_list[[loss_type]][[CNT]] <- loss_cnt_crop_list[[loss_type]][[CNT]] %>% bind_rows()
  }
  loss_list[[loss_type]] <- loss_cnt_list[[loss_type]] %>% bind_rows()
}
# Choose one outcome type here 
df <- loss_list %>% bind_rows() # "loss_cropland"
df[df$model == "all", "model"] <- "Global"
#df$term <- rep(eaear_mapmat[,"Crops"], length(dyn_df_list))
df$significant01 <- ""
df[df[,"Pr(>|t|)"] < 0.1, "significant01"] <- "p-value < .1"

names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"
# decide to remove or not sugarcane 
# x <- x[!grepl(pattern = "expo_X_sugarcane", x = x$term),]
row.names(df) <- NULL
df


#### GROUP6 TESTS -----------------------------------------------------------------------------------------------------

perennials_outcomes <- c("loss_commodity") # , "loss_pasture"
crop_groups <- crops 

# The set of crops which share the same dynamics 
g6_crops <- c("eaear_Banana", "eaear_Cocoa_Coffee", "eaear_Tea", "eaear_Rubber")

g6_crops_est <- list()

for(loss_type in perennials_outcomes){
  for(CNT in c("America", "Africa", "Asia")){# "all", 
    for(CROP in g6_crops){
      # determine the crops to control for 
      crops_ctrl <- crop_groups[sapply(crop_groups, match, x = CROP, nomatch = FALSE)==0] %>% unlist() %>% unname()
      if(loss_type =="loss_pasture"){
        crops_ctrl <- c("eaear_Fodder", "eaear_Citrus", "eaear_Coconut")
        est_parameters[["end_year"]] <- 2016
      }
      # the order of the strata in the list are determiend by how we want to plot results
      g6_crops_est[[loss_type]][[CNT]][[CROP]] <- make_main_reg(continent = CNT,
                                                                      outcome_variable = loss_type,
                                                                      
                                                                      start_year = est_parameters[["start_year"]],
                                                                      end_year = est_parameters[["end_year"]],
                                                                      
                                                                      # regress on one crop at a time, and control for other crops
                                                                      exposure_rfs = CROP,
                                                                      control_all_absolute_rfs = TRUE,
                                                                      annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                                      all_exposures_rfs = crops_ctrl,
                                                                      
                                                                      trade_exposure = est_parameters[["trade_exposure"]],
                                                                      trade_exposure_period = est_parameters[["trade_exposure_period"]],
                                                                      
                                                                      rfs_lead = est_parameters[["leads"]],
                                                                      rfs_lag = est_parameters[["lags"]],
                                                                      rfs_fya = est_parameters[["fya"]],
                                                                      rfs_pya = est_parameters[["pya"]],
                                                                      aggr_lead = est_parameters[["aggr_lead"]],
                                                                      aggr_lag = est_parameters[["aggr_lag"]],
                                                                      s_trend = est_parameters[["s_trend"]],
                                                                      s_trend_loga = est_parameters[["s_trend_loga"]],
                                                                      s_trend_sq = est_parameters[["s_trend_sq"]],
                                                                      
                                                                      output = "coef_table"  )
      
      # g6_crops_est[[loss_type]][[CNT]][[CROP]] <- g6_crops_est[[loss_type]][[CNT]][[CROP]][[3]]
    }
  }    
}

est_obj <- g6_crops_est[[loss_type]][[CNT]][[CROP]][[1]]
est_dat <- g6_crops_est[[loss_type]][[CNT]][[CROP]][[2]]
est_df <- g6_crops_est[[loss_type]][[CNT]][[CROP]][[3]]
unique(est_dat$year)
### PLOT COEFFICIENTS ### 
# prepare regression outputs in a tidy data frame readable by dwplot
rm(df)
loss_cnt_crop_list <- list()
loss_cnt_list <- list()
loss_list <- list()
for(loss_type in perennials_outcomes){
  for(CNT in c("America", "Africa", "Asia")){# "all", , "Africa", "Asia"
    for(CROP in g6_crops){
      x <- g6_crops_est[[loss_type]][[CNT]][[CROP]]%>% as.data.frame()
      x$term <- CROP
      x$model <- CNT
      x$loss_type <- loss_type
      x <- x[grepl("_X_aggrall", row.names(x)), ]
      x$via_trade <- grepl("_expo_", row.names(x))
      loss_cnt_crop_list[[loss_type]][[CNT]][[CROP]] <- x
    }
    loss_cnt_list[[loss_type]][[CNT]] <- loss_cnt_crop_list[[loss_type]][[CNT]] %>% bind_rows()
  }
  loss_list[[loss_type]] <- loss_cnt_list[[loss_type]] %>% bind_rows()
}
# Choose one outcome type here 
df <- loss_list %>% bind_rows() # "loss_cropland"
df[df$model == "all", "model"] <- "Global"
#df$term <- rep(eaear_mapmat[,"Crops"], length(dyn_df_list))
df$significant01 <- ""
df[df[,"Pr(>|t|)"] < 0.1, "significant01"] <- "p-value < .1"

names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"
# decide to remove or not sugarcane 
# x <- x[!grepl(pattern = "expo_X_sugarcane", x = x$term),]
row.names(df) <- NULL
df

{dwplot(df,
        dot_args = list(size = 2.5,aes(shape = model, alpha = significant01)),#, alpha = significant01
        whisker_args = list(size = 1, aes(alpha = significant01)),#alpha = significant01
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(
      eaear_Banana = "Banana", 
      eaear_Cocoa_Coffee = "Cocoa or Coffee",
      eaear_Rubber = "Rubber", 
      eaear_Tea = "Tea"
    )) +
    #scale_color_manual(values=wes_palette(n=4, name="Darjeeling1", type = "discrete")) +
    scale_color_brewer(type = "qual", palette="Set1") +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    #ggtitle(title) +
    theme(plot.title = element_text(face="bold", size=c(10)),
          legend.position = c(0.8, 0.05),
          legend.justification = c(0, 0), 
          legend.background = element_rect(colour="grey80"),
          legend.title = element_blank()) }


