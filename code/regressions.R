
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


### READ ALL POSSIBLE DATASETS HERE 
# but not all together because of memory issues. 
# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_gaez_long_final.Rdata"))

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_gaez_long_final.Rdata"))
# 
# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_gaez_long_final.Rdata"))

prices <- readRDS(here("temp_data", "prepared_prices.Rdata"))

rm(dataset, start_year, end_year, crop_j, j_soy, price_k, extra_price_k, standardized_si, price_lag, SjPj, SkPk, fe, distribution, output, se, cluster, 
   controls, regressors, outcome_variable)

dataset = "glass"
start_year = 1983
end_year = 2020
further_lu_evidence = "none"
crop_j = "Oilpalm"
j_soy = "Soybean_oil"
# for the k variables hypothesized in overleaf for palm oil, feglm quasipoisson converges within 25 iter. 
# but maybe not with skPk controls. 
price_k <- c( "Soybean_oil", "Rapeseed_oil", "Sunflower_oil", 
              "Coconut_oil", "Soybeans","Sugar", "Maize") 
# "Barley",  "Chicken", "Sheep", "Banana", "Beef", "Olive_oil",
#               "Orange",  "Cotton",  "Groundnuts",  "Rubber", "Sorghum","Cocoa",  "Coffee",
#                "Rice",   "Wheat",  "Palm_oil", ), 
#                "Tea", "Tobacco",  "Oat",  
#               , "Pork")
extra_price_k = c("Crude_oil") # , 
standardized_si = TRUE
price_lag = 1
SjPj = TRUE
SkPk = FALSE
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

make_base_reg <- function(dataset = "glass", # one of "glass", "fl8320", "phtfl"
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
  
  
  #### SPECIFICATIONS  
  
  ## Outcome variable
  if(dataset=="glass"){
    outcome_variable <- "first_loss"} #   "sbqt_direct_lu"      "sbqt_mode_lu" 
  if(dataset=="fl8320"){
    outcome_variable <- "firstloss_glassgfc"}
  if(dataset=="phtfl"){
    outcome_variable <- "phtf_loss"}
  
  # restrict the outcome variable to evidence of further LU
  if(dataset == "glass" & further_lu_evidence != "none"){
    if(crop_j == "Pasture"){
      d <- dplyr::mutate(d, lu_evidence = (!!as.symbol(further_lu_evidence) == 30)) # either direct or mode subsequent lu are grassland
    }
    if(crop_j == "Soybean"){
      d <- dplyr::mutate(d, lu_evidence = (!!as.symbol(further_lu_evidence) == 10)) # either direct or mode subsequent lu are cropland
    }
    if(crop_j %in% c("Oilpalm", "Cocoa", "Coffee")){
      d <- dplyr::mutate(d, lu_evidence = (sbqt_direct_lu != 20 & sbqt_mode_lu == 20)) # is not forest in the year directly after, but is forest again in the mode subsequent lu. 
    }
  }
  
  
  
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


#### RUN AND PLOT ####
DS <- "glass"
SY <- 1983
EY <- 2020
if(DS == "glass"){main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_gaez_long_final.Rdata"))}
if(DS == "fl8320"){main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "firstloss8320_gaez_long_final.Rdata"))}
if(DS == "phtfl"){main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_gaez_long_final.Rdata"))}


### PALM OIL
K_palmoil <- c("Soybean_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil",
               "Sugar", "Maize") # in prices spelling "Olive_oil", "Soybeans", "Soybean_meal",

oilpalm_res <- make_base_reg(dataset = DS,
                              start_year = SY, end_year = EY,
                              crop_j = "Oilpalm",
                              price_k = K_palmoil,
                              extra_price_k = "Crude_oil")


### SOY 
K_soy <- c("Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", 
           "Sugar", "Maize") # in prices spelling "Olive_oil",


soy_res <- make_base_reg(dataset = DS,
                        start_year = SY, end_year = EY, 
                        crop_j = "Soybean", 
                        price_k = K_soy, 
                        extra_price_k = "Crude_oil")


### COCOA 
K_cocoa <- c("Coffee", "Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", "Sugar") # in prices spelling


cocoa_res <- make_base_reg(dataset = DS,
                          start_year = SY, end_year = EY, 
                          crop_j = "Cocoa", 
                          price_k = K_cocoa)


### COFFEE 
K_coffee <- c("Tea", "Cocoa", "Sugar", "Tobacco") # in prices spelling


coffee_res <- make_base_reg(dataset = DS,
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
main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_gaez_long_final.Rdata"))

### PALM OIL
K_palmoil <- c("Soybean_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil",
               "Sugar", "Maize") # in prices spelling "Olive_oil", "Soybeans", "Soybean_meal",

oilpalm_res <- make_base_reg(#start_year = 1983, end_year = 2015,
  crop_j = "Oilpalm",
  price_k = K_palmoil,
  extra_price_k = "Crude_oil")


### SOY 
K_soy <- c("Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", 
           "Sugar", "Maize") # in prices spelling "Olive_oil",


soy_res <- make_base_reg(#start_year = 1983, end_year = 2015, 
  crop_j = "Soybean", 
  price_k = K_soy, 
  extra_price_k = "Crude_oil")


### COCOA 
K_cocoa <- c("Coffee", "Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", "Sugar") # in prices spelling


cocoa_res <- make_base_reg(#start_year = 1983, end_year = 2015, 
  crop_j = "Cocoa", 
  price_k = K_cocoa)


### COFFEE 
K_coffee <- c("Tea", "Cocoa", "Sugar", "Tobacco") # in prices spelling


coffee_res <- make_base_reg( #start_year = 1983, end_year = 2015, 
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
main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "phtfloss_gaez_long_final.Rdata"))

### PALM OIL
K_palmoil <- c("Soybean_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil",
               "Sugar", "Maize") # in prices spelling "Olive_oil", "Soybeans", "Soybean_meal",

oilpalm_res <- make_base_reg(dataset = "phtfl",
                             start_year = 2002, end_year = 2020,
                             crop_j = "Oilpalm",
                             price_k = K_palmoil,
                             extra_price_k = "Crude_oil")


### SOY 
K_soy <- c("Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", 
           "Sugar", "Maize") # in prices spelling "Olive_oil",


soy_res <- make_base_reg(dataset = "phtfl",
                         start_year = 2002, end_year = 2020, 
                         crop_j = "Soybean", 
                         price_k = K_soy, 
                         extra_price_k = "Crude_oil")


### COCOA 
K_cocoa <- c("Coffee", "Palm_oil", "Rapeseed_oil", "Sunflower_oil", "Coconut_oil", "Sugar") # in prices spelling


cocoa_res <- make_base_reg(dataset = "phtfl",
                           start_year = 2002, end_year = 2020, 
                           crop_j = "Cocoa", 
                           price_k = K_cocoa)


### COFFEE 
K_coffee <- c("Tea", "Cocoa", "Sugar", "Tobacco") # in prices spelling


coffee_res <- make_base_reg(dataset = "phtfl",
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



