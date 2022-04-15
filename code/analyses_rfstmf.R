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
                   "ggplot2", "dotwhisker", "viridis", "hrbrthemes", "rasterVis", #"tmap",# "leaflet", "htmltools"
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
dir.create(here("temp_data","reg_results", "rfs"))



### GLOBAL CRS USED ### 
mercator_world_crs <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "

### OBJECTS USED IN RFS PROCESSES ###


# this is mapmat for eaear 
eaear_mapmat_data <- c(
  "Banana","eaear_Banana",
  "Crude_oil", "eaear_Biomass",
  #"Barley", "eaear_Barley",
  "Cereals", "eaear_Cereals",
  "Orange", "eaear_Citrus",
  #"Cocoa", "eaear_Cocoa",
  "Cocoa_Coffee", "eaear_Cocoa_Coffee",
  "Coconut_oil", "eaear_Coconut",
  #"Coffee", "eaear_Coffee",
  "Cotton", "eaear_Cotton",
  #"Groundnuts", "eaear_Groundnut",
  "Beef", "eaear_Fodder", 
  "Maize", "eaear_Maizegrain",
  "Oilfeed_crops", "eaear_Oilfeed_crops",
  # "Oat", "eaear_Oat",
  # "Olive_oil", "eaear_Olive",  
  "Palm_oil", "eaear_Oilpalm",
  #"Rapeseed_oil", "eaear_Rapeseed",
  "Rice", "eaear_Rice",
  "Rubber", "eaear_Rubber",
  #"Sorghum", "eaear_Sorghum2", 
  "Soy_index", "eaear_Soy_compo",
  "Sugar", "eaear_Sugar", 
  #"Sunflower_oil", "eaear_Sunflower",
  "Tea", "eaear_Tea",
  "Tobacco", "eaear_Tobacco" 
  #"Wheat", "eaear_Wheat"
)

eaear_mapmat <- matrix(data = eaear_mapmat_data, 
                       nrow = length(eaear_mapmat_data)/2,
                       ncol = 2, 
                       byrow = TRUE)

colnames(eaear_mapmat) <- c("Prices", "Crops")


### MAIN DATA SET ### 
main_data <- readRDS(here("temp_data", "merged_datasets", "tmf_aoi", "tmf_aeay_pantrop_long_final_1990_2020.Rdata"))
# release some memory upfront
main_data <- dplyr::filter(main_data, year >= 2011, year <= 2019)

# main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "driverloss_aesi_long_final.Rdata"))

### PRICE DATA ###
prices <- readRDS(here("temp_data", "prepared_international_prices.Rdata"))



# renewable fuel standards, as from https://www.epa.gov/renewable-fuel-standard-program/renewable-fuel-annual-standards
# remember: biodiesel are nested within advanced biofuels, which are themeselves nested within total
# conventional biofuels are the part of the total that is not "advanced"
rfs <- data.frame(year = min(prices$year):2022, 
                  statute_total = 0, final_total = 0, 
                  statute_advanced = 0, final_advanced = 0, 
                  statute_biodiesel = 0, final_biodiesel = 0, 
                  statute_conv_earliest = 0)

rfs <- dplyr::arrange(rfs, year)
# this is if RFS2 takes over as soon as 2008 (9bgal)
rfs[rfs$year >= 2006 & rfs$year <= 2022, c("statute_total")] <- c(4, 4.7, 9, 11.1, 12.95,	13.95,	15.2,	16.55,	18.15,	20.5,	22.25,	24.0,	26.0,	28.0, 30.0, 33.0, 36.0)
rfs[rfs$year >= 2006 & rfs$year <= 2022, c("final_total")] <- c(4, 4, 9, 11.1, 12.95,	13.95,	15.2,	16.55,	16.28,	16.93,	18.11,	19.28,	19.29,	19.92, 20.09, NA, NA)

rfs[rfs$year >= 2006 & rfs$year <= 2022, c("statute_advanced")] <- c(0, 0, 0, 0.6, 0.95,	1.35,	2.0,	2.75,	3.75,	5.5,	7.25,	9.0,	11.0,	13.0, 15.0, 18.0, 21.0)
rfs[rfs$year >= 2006 & rfs$year <= 2022, c("final_advanced")] <- c(0, 0, 0, 0.6, 0.95,	1.35,	2.0,	2.75,	2.67,	2.88,	3.61,	4.28,	4.29,	4.92, 5.09, NA, NA)

rfs[rfs$year >= 2009 & rfs$year <= 2022, c("final_biodiesel")] <- c(0.5, 1.15, 0.8, 1.0, 1.28, 1.63, 1.73, 1.9, 2.0, 2.1, 2.1, 2.43, 2.43, NA)
rfs[rfs$year >= 2009 & rfs$year <= 2022, c("statute_biodiesel")] <- c(0.5, 0.65, 0.8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

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
  for(lead in c(1:6)){
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
  
  for(lag in c(1:6)){
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

#### DES STATS RFS #### 
# add PSD data to the chart 
psd <- readRDS(here("temp_data", "prepared_psd.Rdata"))
psd$us_ah_maize_mha <- psd$UnitedStates.Area_Harvested.Maize / 1000

#head(psd)
w_rfs <- left_join(rfs, psd[,c("year", "us_ah_maize_mha")], by = "year")

w_rfs <- w_rfs %>% 
  dplyr::filter(year<=2020) %>% 
  dplyr::select(year, statute_advanced, statute_conv, final_advanced, final_conv, 
                us_ah_maize_mha) %>% 
  dplyr::mutate(statute_conv = statute_conv - final_conv, 
                statute_advanced = statute_advanced - final_advanced) %>%
  tidyr::gather(RFS, BOG, c("statute_advanced", "statute_conv", "final_advanced", "final_conv")) 


# rename RFS
w_rfs[w_rfs$RFS == "statute_advanced", "RFS"] <- "Advanced biofuels, statutory"
w_rfs[w_rfs$RFS == "statute_conv", "RFS"] <- "Conventional biofuels, statutory"
w_rfs[w_rfs$RFS == "final_advanced", "RFS"] <- "Advanced biofuels, final"
w_rfs[w_rfs$RFS == "final_conv", "RFS"] <- "Conventional biofuels, final"

# factor to give a specific order
w_rfs$RFS <- factor(w_rfs$RFS, 
                    levels = c("Advanced biofuels, statutory", 
                               "Advanced biofuels, final", 
                               "Conventional biofuels, statutory", 
                               "Conventional biofuels, final"))
# rename axis
names(w_rfs)[names(w_rfs)=="RFS"] <- "RFS mandates (right axis)"
names(w_rfs)[names(w_rfs)=="BOG"] <- "Billions of gallons (ethanol-equivalent)"

coeff <- 0.4

ggplot(w_rfs, aes(x=year, y=`Billions of gallons (ethanol-equivalent)`/coeff, fill=`RFS mandates (right axis)`)) + 
  geom_area(alpha=0.6 , size=.5, colour="white") +
  geom_line(aes(x = year, y = us_ah_maize_mha)) +
  scale_y_continuous(name = "Maize area harvested in US, Mha",
                     sec.axis = sec_axis(~.*coeff, name="Billions of gallons (ethanol-equivalent)")) +# 
  scale_fill_viridis(discrete = T) +
  theme_ipsum() + 
  #ggtitle("RFS mandates, statutory and final, ramping up for advanced and conventinal biofuels") + 
  theme(plot.title = element_text(size = 10, face = "bold"), 
        axis.title.y.left=element_text(size=10,face="bold", hjust = 1),
        axis.title.y.right=element_text(size=10,face="bold", hjust = 1),
        axis.title.x=element_text(size=10,face="bold", hjust = 0.5))


#### DES STATS EAEAR #### 

### MAP OF TROPICAL COMMODITY DRIVEN DEFORESTATION ### 
prepared_maps <- list(America = list(), 
                      Africa = list(), 
                      Asia = list())
for(CNT in c("America", "Africa", "Asia")){
  fullstack <- brick(here("temp_data", "merged_datasets", "tmf_aoi", paste0("anytype_masked_stack_",CNT,"_1990_2020.tif")) )
  agri <- fullstack[[1:31]]
  plan <- fullstack[[32:62]]
  floo <- fullstack[[63:93]]
  tmfx <- fullstack[[94:124]]
  gaez <- fullstack[[125:174]]
  pst2k <- fullstack[[175]]
  
  names(agri) <- paste0("tmf_agri.",seq(1990, 2020, 1)) 
  names(plan) <- paste0("tmf_plantation.",seq(1990, 2020, 1)) 
  names(floo) <- paste0("tmf_flood.",seq(1990, 2020, 1)) 
  names(tmfx) <- paste0("tmf_ext.",seq(1990, 2020, 1)) 
  
  agri_accu <- sum(agri[[22:30]], na.rm = TRUE)
  plan_accu <- sum(plan[[22:30]], na.rm = TRUE)
  floo_accu <- sum(floo[[22:30]], na.rm = TRUE)
  
  # vcdl <- values(agri_accu)
  # summary(vcdl)
  # rm(vcdl)
  # turn zeros into NAs
  agri_accu <- reclassify(agri_accu, cbind(0, NA))
  #plot(agri_accu)
  a <- area(agri_accu)
  agri_accu <- stack(agri_accu, a)
  agri_accu_pct <- overlay(agri_accu, fun = function(x, y){pct <- x/(100*y)
                                                           pct <- if_else(pct>1, 1, pct)
                                                           return(pct)})
  # vpct <- values(agri_accu_pct)
  # summary(vpct)
  
  plan_accu <- reclassify(plan_accu, cbind(0, NA))
  #plot(plan_accu)
  a <- area(plan_accu)
  plan_accu <- stack(plan_accu, a)
  plan_accu_pct <- overlay(plan_accu, fun = function(x, y){pct <- x/(100*y)
                                                           pct <- if_else(pct>1, 1, pct)
                                                           return(pct)})

  floo_accu <- reclassify(floo_accu, cbind(0, NA))
  #plot(floo_accu)
  a <- area(floo_accu)
  floo_accu <- stack(floo_accu, a)
  floo_accu_pct <- overlay(floo_accu, fun = function(x, y){pct <- x/(100*y)
                                                           pct <- if_else(pct>1, 1, pct)
                                                           return(pct)})
  
  prepared_maps[[CNT]][["agri"]] <- agri_accu_pct
  prepared_maps[[CNT]][["plan"]] <- plan_accu_pct
  prepared_maps[[CNT]][["floo"]] <- floo_accu_pct
  
  rm(a, agri_accu, plan_accu, floo_accu, agri_accu_pct, plan_accu_pct, floo_accu_pct)
}

### MAKE CONTINENT PANEL MAPs 
land <- st_read(here("input_data", "ne_50m_land"))
unique(land$scalerank)
land <- land[land$scalerank==0, c("geometry")]
#plot(land)
spLand <- as(land, "Spatial")

am <- prepared_maps[["America"]][["floo"]]
af <- prepared_maps[["Africa"]][["floo"]]
as <- prepared_maps[["Asia"]][["floo"]]

#library(gridExtra)
library(rasterVis)

pam <- levelplot(am, margin = FALSE, 
                 xlab = "", ylab = "",
                 #colorkey=FALSE,
                 #scales=list(draw = FALSE),            
                 col.regions=magma(n = 10000, direction = -1),                   
                 #at=seq(0, 10000),
                 names.attr=rep('', 1)) + 
  layer(sp.polygons(spLand, lwd=1)) 


paf <-levelplot(af, margin = FALSE,
                xlab = "", ylab = "",
                colorkey=FALSE,
                # scales=list(draw = FALSE), x=list(draw=FALSE)           
                col.regions=magma(n = 10000, direction = -1),                   
                #at=seq(0, 10000),
                names.attr=rep('', 1)) + 
  layer(sp.polygons(spLand, lwd=1)) 

pas <-levelplot(as, margin = FALSE, 
                xlab = "", ylab = "",
                colorkey = FALSE,
                # colorkey=list(
                #   space='bottom',                   
                #   labels=list(at=-5:5, font=4),
                #   axis.line=list(col='black'),
                #   width=0.75
                # ),    
                # par.settings=list(
                #   strip.border=list(col='transparent'),
                #   strip.background=list(col='transparent'),
                #   axis.line=list(col='transparent')
                # ),names.attr=rep('', 1),
                #scales=list(draw = FALSE),            
                col.regions=magma(n = 10000, direction = -1),                   
                #at=seq(0, 10000),
) + 
  layer(sp.polygons(spLand, lwd=1)) 

# It's necessary to call the objects to display the maps

pam
paf
pas
# make some comparisons with driven loss 
drivenloss  <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "driverloss_all_aeay_long_final.Rdata"))

drivenloss <- dplyr::filter(drivenloss, year >=2011 & year <= 2019)

drivenloss_den <- density(drivenloss$driven_loss_commodity)

main_data <- dplyr::mutate(main_data, tmf_deforestation = tmf_agri + tmf_plantation)

tmf_agri_den <- density(main_data$tmf_deforestation)




#### REGRESSION FUNCTION #### 

### TEMPORARY OBJECTS 
outcome_variable = "tmf_agri" # "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
start_year = 2011
end_year = 2019
continent = "America"
further_lu_evidence = "none"
original_sj = "Fodder"# ,"Fodder",  "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber" in GAEZ spelling, one, part, or all of the 6 main drivers of deforestation: 

pasture_shares <- FALSE

estimated_effect = "alpha"
rfs_reg <- TRUE
rfs_rando <- ""
original_rfs_treatments <- c("statute_conv")
rfs_lead <- 3
rfs_lag <- 3 
rfs_fya <-  0
rfs_pya <- 0
aggr_dyn <- TRUE
original_exposure_rfs <- "eaear_Soy_compo"
group_exposure_rfs <- FALSE
control_all_absolute_rfs <- TRUE
annual_rfs_controls <- FALSE

control_pasture <- FALSE
pasture_trend <- FALSE

fc_trend <- FALSE
s_trend <- TRUE
fc_s_trend <- FALSE

sjpos <- FALSE # should the sample be restricted to cells where sj is positive? 

fe = "grid_id + country_year" #   
distribution <- "quasipoisson"
offset <- FALSE
invhypsin = TRUE
conley_cutoff <- 100
clustering = "twoway" # either "oneway" or "twoway". If oneway, it clusters on cluster_var1. 
cluster_var1 = "grid_id" 
cluster_var2 = "grid_id_5_year"

output = "est_object"
glm_iter <- 25 




rm(outcome_variable, start_year, end_year, continent, further_lu_evidence, original_sj, original_Pk, pasture_shares, 
   standardization, price_dyn, estimated_effect, group_prices, biofuel_focus, aggregate_K, control_interact_sj, control_interact_Pk, reference_crop, control_direct, 
   rfs_reg, original_rfs_treatments, rfs_lead, rfs_lag, original_exposure_rfs, group_exposure_rfs, control_absolute_rfs, control_all_absolute_rfs, remaining, sjpos, fe, distribution, invhypsin, conley_cutoff, se, boot_cluster, coefs_to_aggregate, 
   output, glm_iter)

make_main_reg <- function(outcome_variable = "tmf_agri", # one of "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
                          start_year = 2011, 
                          end_year = 2019, 
                          continent = "all", # one of "Africa", "America", "Asia", or "all"
                          further_lu_evidence = "none", # either "none", "sbqt_direct_lu", or "sbqt"mode_lu"
                          original_sj = c("Fodder"),  # in GAEZ spelling, one, part, or all of the 6 main drivers of deforestation: 
                          all_but_k = FALSE,
                          # , "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber"
                          # note that this is the vector that defines what is removed from full_control
                          original_Pk = c("Beef"), # in price spelling. One, part, or all of the full set of commodities having a price-AESI match.
                          # , "Palm_oil", "Cocoa", "Coffee", "Rubber", 
                          # "Rapeseed_oil", "Sunflower_oil","Rice", "Wheat", "Maize", "Sugar", "Sorghum"
                          # focal_j_extra_price = "Soybean", # in GAEZ spelling, the crop of interest if the treatment is price that has no match in GAEZ (an "extra price" in previous saying)
                          # extra_price_k = c(), # "Chicken", "Pork", "Sheep", "Crude_oil" in price spelling. One, part, or all of the full set of commodities NOT having a price-AESI match.
                          
                          pasture_shares = FALSE, # if TRUE, and crop_j = "Fodder", then qj is proxied with the share of pasture area in 2000. 
                          standardization = "_std2", # one of "", "_std", or "_std2"
                          price_dyn = "_lag1", # one of main - and then prices are lagged 1 year if sj is fodder or soybean and if it's tree plantations it is 3pya or "alt" and it is the other way round.
                          # available price dynamics are "_lag1", "_2pya", "_3pya", "_4pya", "_5pya",
                          estimated_effect = "alpha",# if "alpha", estimates the (aggregated or not) cross-elasticity effect. If different from "alpha" (e.g. "gamma") it estimates the ILUC effect from price covariation. 
                          rfs_reg = TRUE,
                          rfs_rando = "", # either "between", "within", or any other string. If one of the former two, randomization inference of the type is performed
                          original_rfs_treatments = c("statute_conv"),
                          rfs_lead = 0,
                          rfs_lag = 0,
                          rfs_fya = 0, 
                          rfs_pya = 0,
                          aggr_dyn = TRUE, # whether to report aggregate coefficients of all leads and lags ("all") or leads and lags separately ("leadlag"), or no aggregate (any other string)
                          original_exposure_rfs = "eaear_Soy_compo",
                          group_exposure_rfs = FALSE, # if original_exposure_rfs is length 1, it does noe matter whether this option is TRUE or FALSE
                          control_all_absolute_rfs = TRUE,
                          annual_rfs_controls = FALSE,
                          
                          control_pasture = FALSE,
                          pasture_trend = FALSE,
                          
                          aggregate_K = "", # if not "", either "meats", "cereals_feeds", "vegetable_oils", or "biofuel_feedstocks"
                          control_interact_sj = FALSE, 
                          group_prices = FALSE,
                          biofuel_focus = FALSE, # how to group prices? With a focus on biofuel feedstocks or rather on cereals/feeds and vegetable oils? 
                          
                          control_interact_Pk = FALSE,
                          reference_crop = "Olive",
                          
                          control_direct = FALSE,
                          control_sjPj = FALSE,
                          incl_sjPj = FALSE,
                          control_skPk = FALSE,
                          incl_skPk = FALSE,
                          # sjpj_lag = "_lag1", # either "" or "_lag1" or "_lag2"
                          # skpk_lag = "_lag1", # either "" or "_lag1" or "_lag2"
                          remaining = FALSE, # should remaining forest be controlled for STOP DOING THIS BECAUSE IT INTRODUCES NICKELL BIAS
                          s_trend = TRUE,
                          fc_trend = FALSE,
                          fc_s_trend = FALSE,
                          
                          sjpos = FALSE, # should the sample be restricted to cells where sj is positive? 
                          # open_path = FALSE,
                          # commoXcommo = "Fodder_X_Beef",
                          #commo_m = c(""), comment coder ça pour compatibilité avec loops over K_commo ? 
                          fe = "grid_id + country_year", 
                          distribution = "quasipoisson",#  "quasipoisson", 
                          invhypsin = TRUE, # if distribution is gaussian, should the dep. var. be transformed to inverse hyperbolic sine?
                          offset = FALSE,
                          se = "twoway",# # passed to vcov argument. Currently, one of "cluster", "twoway", or an object of the form: 
                          # - "exposure2ways" for the two-way clustering to be customed as original_exposure_rfs + country_year, and not along FE.   
                          # - vcov_conley(lat = "lat", lon = "lon", cutoff = 100, distance = "spherical")
                          # with cutoff the distance, in km, passed to fixest::vcov_conley, if se = "conley"  
                          # boot_cluster ="grid_id",
                          # old argument: cluster ="grid_id", # the cluster level if se = "cluster" (i.e. one way)
                          # coefstat = "confint", # one of "se", "tstat", "confint"
                          glm_iter = 25,
                          output = "coef_table" # one of "data", est_object, or "coef_table" 
){
  
  # Define the outcome_variable based on the crop under study (if we are not in the placebo case)
  if(original_exposure_rfs %in% c("eaear_Oilpalm", "eaear_Rubber") & outcome_variable != "tmf_flood"){
    outcome_variable <- "tmf_plantation"
  }
  
  # and this matrix is used to select exposure type, either eaear or aesi 
  if(grepl("eaear_", original_exposure_rfs) & rfs_reg){
    exp_matmap <- eaear_mapmat
  }else{
    exp_matmap <- mapmat_si
  }
  #### PREPARE NEEDED VARIABLE NAMES
  # this does not involve data, just arguments of the make_reg function
  # original_ names are used to get generic covariate names in coefficient tables (irrespective of modelling choices)
  
  ## Standardized suitability index to find in the main data 
  # this just those for which focus is set in current specification
  sj <- paste0(original_sj, standardization)  
  # this is all possible exposures, necessary in every specificaton
  original_all_exposures <- exp_matmap[,"Crops"]
  
  # (maybe for alpha only currently) identify the L - 1 exposures that will be interacted with the treatment.
  # this is necessary to avoid perfect colinearity with time FE
  # if(original_treatments %in% mapmat_si[,"Prices"]){
  #   original_all_exposures_but_k <- original_all_exposures[original_all_exposures != mapmat_si[mapmat_si[,"Prices"]==original_treatments, "Crops"]]
  # }else{
  #   # handles cases when treatments is an extra price with no match in GAEZ 
  #   original_all_exposures_but_k <- original_all_exposures[original_all_exposures != focal_j_extra_price]
  # }
  
  # standardize them
  all_exposures <- paste0(original_all_exposures, standardization)
  
  exposure_rfs <- paste0(original_exposure_rfs, standardization)
  
  ## Price variable names to find in the price data
  price_info <- price_dyn 
  # if(price_dyn == "main"){
  #   if(original_sj %in% c("Fodder", "Soybean")){
  #     price_info <- "_lag1"
  #   }else{
  #     price_info <- "_3pya"
  #   }
  # }
  # if(price_dyn == "alt"){
  #   if(original_sj %in% c("Fodder", "Soybean")){
  #     price_info <- "_3pya"
  #   }else{
  #     price_info <- "_lag1"
  #   }
  # }
  
  # actual name to find in data
  Pk <- paste0("ln_", original_Pk, price_info)
  
  # it's corresponding exposure
  original_sk <- mapmat_si[mapmat_si[,"Prices"]==original_Pk, "Crops"]
  sk <- paste0(original_sk, standardization)
  
  # individual prices mobilized in this specific sj estimation
  original_Pj <- mapmat_si[mapmat_si[,"Crops"]==original_sj, "Prices"]
  Pj <- paste0("ln_", original_Pj, price_info)
  
  original_treatments <- unlist(maplist[[original_sj]])
  treatments <- paste0("ln_", original_treatments, price_info)
  
  # this is all possible treatments, necessary in every specification 
  original_all_treatments <- mapmat_si[,"Prices"]
  all_treatments <- paste0("ln_", original_all_treatments, price_info)
  
  all_rfs_treatments <- grep(pattern = original_rfs_treatments, names(prices), value = TRUE)
  
  # weighted prices, necessary if control_interact_sj (we don"t need an original version)
  w_treatments <- c(paste0("w_kca_ln_", vegetable_oils, price_info), paste0("w_kca_ln_", cereals_feeds, price_info))
  # note that "Soybean" is not one of them. Only Soybean_meal and Soybean_oil. 
  
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
  
  
  
  # this is used at some points 
  sjPj <- paste0(original_sj, "_X_", original_Pj)
  skPk <- paste0(original_sk, "_X_", original_Pk)
  sjPk <- paste0(original_sj, "_X_", original_Pk)
  
  
  #### MAKE THE VARIABLES NEEDED IN THE DATA
  # manipulate a different data set so that original one can be provided to all functions and not read again every time. 
  d <- main_data
  
  # Keep only in data the useful variables 
  d <- dplyr::select(d, all_of(c("grid_id", "year", "lat", "lon", "country_name", "country_year", "continent_name", 
                                 # "fc_2000", "fc_2009", "remaining_fc", # "accu_defo_since2k",
                                 outcome_variable, "pasture_share_2000",
                                 unique(c(original_all_exposures, all_exposures))))) #sj, 
  
  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", unique(c(treatments, all_treatments, all_rfs_treatments, w_treatments, rfs_treatments)))], by = c("year"))#, all_treatments
  
  
  # If we want the exposure to deforestation for pasture to be proxied with the share of pasture area in 2000, rather than suitability index, 
  if(pasture_shares){
    d[,grepl("Fodder", names(d))] <- d$pasture_share_2000
  }
  
  if((distribution == "gaussian" | estimated_effect == "beta") & invhypsin){
    # transform dependent variable, if gaussian GLM 
    d <- dplyr::mutate(d, !!as.symbol(outcome_variable) := asinh(!!as.symbol(outcome_variable)))
  }
  
  # /!\ important! : if this is TRUE, then sj and original_sj are still CALLED the same way, so referring to crop j. 
  # Yet, the value now is (1-sk)  
  if(all_but_k){
    
    d <- mutate(d, !!as.symbol(sj) := (1 - !!as.symbol(sk)))
    
    # in this case, we don't want to group prices and to interact them with sj 
    control_interact_sj <- FALSE
    
    # we don't want skPk to be controlled for, even in the dir_ctrl_grp
    control_skPk <- FALSE
    
    # and we don't want to aggregate over K
    aggregate_K <- ""
  }
  
  
  ### DEFINE SPECIFICATION ### 
  regressors <- sjPk
  
  # d <- mutate(d,
  #             !!as.symbol(sjPk) := (!!as.symbol(sj)) * (!!as.symbol(Pk)) )
  
  # if(original_sj == "Soybean"){
  #   sjPsoy <- c()
  #   for(soy_commo in unlist(maplist[[original_sj]][["soy"]])){
  #     varname <- paste0(original_sj, "_X_", soy_commo)
  #     sjPsoy <- c(sjPsoy, varname)
  #     Pks <- paste0("ln_", soy_commo, price_info)
  #     d <- mutate(d,
  #                 !!as.symbol(varname) := (!!as.symbol(sj)) * (!!as.symbol(Pks)) )
  #   }
  #   rm(soy_commo, Pks, varname)
  # 
  #   d <- mutate(d,
  #               Soybean_X_soy := rowMeans(across(.cols = (any_of(sjPsoy)))))
  # 
  #   regressors <- c(regressors, "Soybean_X_soy")
  # }    
  
  ## RFS 
  
  if(rfs_reg){
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
        for(exp_rfs in exposure_rfs){
          # make regressors of interest
          varname <- paste0(exp_rfs, "_X_", rfs_var)
          regressors <- c(regressors, varname)
          d <- mutate(d, !!as.symbol(varname) := !!as.symbol(exp_rfs) * !!as.symbol(rfs_var))
        }
      }
    }
  }
  
  
  if(control_sjPj & !incl_sjPj & (original_Pj != original_Pk)){
    regressors <- c(regressors, sjPj)
    d <- mutate(d,
                !!as.symbol(sjPj) := (!!as.symbol(sj)) * (!!as.symbol(Pj)) )
  }
  
  if(control_skPk & !incl_skPk & (original_Pj != original_Pk) & original_Pk %in% mapmat_si[,"Prices"]){
    regressors <- c(regressors, skPk)
    d <- mutate(d,
                !!as.symbol(skPk) := (!!as.symbol(sk)) * (!!as.symbol(Pk)) )
  }
  
  # if Pk is one of the two soy commodities, 
  if((original_Pk == "Soybean_meal" | original_Pk == "Soybean_oil")){ 
    
    # we don't control for Soybean_X_Soy_index (which is in there only if sj is soybean)
    if(original_sj == "Soybean"){
      regressors <- regressors[regressors != "Soybean_X_Soy_index"]
      # if moreover, we control for sjPj, add the interaction between suitability for Soybean and the other soy commodity than that of interest. 
      if(control_sjPj & !incl_sjPj){
        original_Psoy <- soy[soy != original_Pk]
        Psoy <- paste0("ln_", original_Psoy, price_info)
        varname <- paste0("Soybean_X_", original_Psoy)
        regressors <- c(regressors, varname)
        d <- mutate(d,
                    !!as.symbol(varname) := (!!as.symbol(sj)) * (!!as.symbol(Psoy)) )
      }
    }else{ # if sj is not soybean and we control for skPk, then simply add its interaction with Soybean suitability  
      if(control_skPk & !incl_skPk){
        varname <- paste0("Soybean_X_", original_Pk)
        regressors <- c(regressors, varname)
        s_soy <- paste0("Soybean", standardization)
        d <- mutate(d,
                    !!as.symbol(varname) := (!!as.symbol("s_soy")) * (!!as.symbol(Pk)) )#s_soy Soybean
      }
    }
  }
  
  if(control_interact_sj){
    # useful for later
    protected <- regressors
    
    # remove Pk from crops to group (bc sjPk is the regressor of interest)
    original_treatments_but_k <- original_treatments[original_treatments != original_Pk]
    treatments_but_k <- paste0("ln_", original_treatments_but_k, price_info)
    
    # remove Pj from crops to group, because in any case we want sjPj to be either as a stand-alone control or absorbed in direct controls
    original_treatments_but_jk <- original_treatments_but_k[original_treatments_but_k != original_Pj]
    treatments_but_jk <- paste0("ln_", original_treatments_but_jk, price_info)
    
    # if sj is Soybean, we want to remove both soy commodities from treatments to interact with sj: 
    if(original_sj == "Soybean"){
      original_treatments_but_jk <- original_treatments_but_jk[!grepl("Soybean", original_treatments_but_jk)]
      treatments_but_jk <- paste0("ln_", original_treatments_but_jk, price_info)
    }
    
    # in this case too, we want to remove both soy prices from their groups - but currently this case does not arise, as Soy_index is not in maplist thus not in treatmments
    if(original_Pk=="Soy_index"){
      original_treatments_but_jk <- original_treatments_but_jk[!grepl("Soybean", original_treatments_but_jk)]
      treatments_but_jk <- paste0("ln_", original_treatments_but_jk, price_info)
    }
    # unnecessary currently, because Pj is Soy_index if sj is Soybean, and soy_index is never in maplist 
    
    
    # make groups specific to the present sj
    sj_group_list <- list()
    for(sj_grp in names(maplist[[original_sj]])){
      sj_group_list[[sj_grp]] <- treatments_but_jk[original_treatments_but_jk %in% maplist[[original_sj]][[sj_grp]]]
    }
    
    # sj_group_list <- list(meats = treatments_but_j[original_treatments_but_j %in% maplist[[original_sj]][["meats"]]], 
    #                       cereals_feeds = treatments_but_j[original_treatments_but_j %in% maplist[[original_sj]][["cereals_feeds"]]], 
    #                       vegetable_oils = treatments_but_j[original_treatments_but_j %in% maplist[[original_sj]][["vegetable_oils"]]], 
    #                       biofuel_feedstocks = treatments_but_j[original_treatments_but_j %in% maplist[[original_sj]][["biofuel_feedstocks"]]],
    
    used_group_names <- names(sj_group_list[lengths(sj_group_list)>0]) # currently useless
    for(grp in used_group_names){
      # average the weighted logarithms of prices (preserving the price info in the name, and thus replacing non-grouped crops)
      grp_commodities <- sj_group_list[[grp]]
      
      # if we are combining crops that are weighted, then SUM them  (prices for these crops are expressed for $/kca)
      if(grp %in% c("vegetable_oils", "cereals_feeds")){
        w_kca_grp_commodities <- paste0("w_kca_", grp_commodities)
        
        d <- mutate(d, 
                    !!as.symbol(paste0("ln_", grp, price_info)) := rowSums(across(.cols = (any_of(w_kca_grp_commodities)))))
      }else{
        # if we are not to use weighted prices, then AVERAGE them
        d <- mutate(d, 
                    !!as.symbol(paste0("ln_", grp, price_info)) := rowMeans(across(.cols = (any_of(grp_commodities)))))
      }
    }
    # log Pj and Pk 
    # if(Pj != Pk){
    #   d <- mutate(d, !!as.symbol(Pk) := log(!!as.symbol(Pk)), 
    #               !!as.symbol(Pj) := log(!!as.symbol(Pj)))
    # }else{
    #   d <- mutate(d, !!as.symbol(Pk) := log(!!as.symbol(Pk)))
    # }    
    if(group_prices){  
      # interact all the combined treatments with the exposure 
      for(original_gn in used_group_names){
        varname <- paste0(original_sj, "_X_", original_gn)
        regressors <- c(regressors, varname)
        
        gn <- paste0("ln_", original_gn, price_info)
        d <- mutate(d,
                    !!as.symbol(varname) := (!!as.symbol(sj)) * (!!as.symbol(gn)) )
      }
      rm(varname, original_gn, gn)
      
      # finally, remove interactions with prices that are not correlated with Pk
      protected <- regressors[regressors %in% protected |
                                regressors == paste0(original_sj, "_X_cereals_feeds") |
                                regressors == paste0(original_sj, "_X_vegetable_oils")]
      
      unprotected <- regressors[!(regressors %in% protected)]
      regressors <- c(protected, unprotected[(unprotected %in% paste0(original_sj, "_X_", k_treatments[[k_price]]))])
      
    }else{ # i.e. if we do not group terms of interaction between sj and other prices 
      
      for(tr in treatments_but_jk){
        varname <- paste0(original_sj, "_X_", original_treatments_but_jk[match(tr, treatments_but_jk)])
        regressors <- c(regressors, varname)
        
        # if we are to weight covariates by the size of the respective market
        # w_tr <- paste0("w_", tr)
        # if(w_tr %in% names(d)){
        #   tr <- w_tr
        # }
        
        d <- mutate(d,
                    !!as.symbol(varname) := (!!as.symbol(sj)) * ( (!!as.symbol(tr)) ))
      }
      rm(varname, tr)
      
    }  
  }
  
  
  ### HANDLE SPECIFICATION TO AGGREGATE EFFECTS ####
  if(aggregate_K!=""){
    original_Pks <- unlist(maplist[[original_sj]][[aggregate_K]])
    
    # remove Pj from the Pks (necessary when the group to aggregate over includes Pj)
    # if original_sj == "Soybean" this is currently not needed because handled in maplist directly. But this is necessary for other sj.   
    original_Pks <- original_Pks[original_Pks != original_Pj]
    
    # coefficients to aggregate (necessary in post estimation below)
    sjPks <- paste0(original_sj, "_X_", original_Pks)
    
    # In any case, we may need these objects
    sjPgrp <- c()
    for(commo in original_Pks){
      varname <- paste0(original_sj, "_X_", commo)
      sjPgrp <- c(sjPgrp, varname)
      Pks <- paste0("ln_", commo, price_info)
      d <- mutate(d,
                  !!as.symbol(varname) := (!!as.symbol(sj)) * (!!as.symbol(Pks)) )# don't log Pk here, it as been pre logged above, if soy is specified as a group of interest for sj
    }
    rm(commo, Pks, varname)
    
    # remove any individual sjPk that may already be in the specification
    regressors <- regressors[!(regressors %in% sjPgrp)] 
    
    # We don't do this currently, as we prefer specifying sjPj as a single variable, and the interaction with the group it belongs to excluding it. 
    # UNLESS THE GROUP INCLUDES SOY AND sj IS SOYBEAN, hence the condition below
    # # if Pj is NOT in aggregate_K, sjPgrp does not include sjPj 
    # if(!(original_Pk %in% original_Pks)){
    #   # hence, in order to be able to aggregate properly over groups that do not include Pj, we need to remove sjPj from regressors
    #   regressors <- regressors[regressors != paste0(original_sj, "_X_", original_Pj)]
    #   # and recompute the interaction with the group that includes Pj, and that is used as a control in this specification 
    #   # this group is 
    #   group_Pj <- names(maplist[[original_sj]])[sapply(maplist[[original_sj]], FUN = function(x){return(original_Pj %in% x)})]
    #   sjPgrp_j <- paste0(original_sj, "_X_", unlist(maplist[[original_sj]][[group_Pj]]))
    #   
    #   d <- mutate(d, 
    #               !!as.symbol(paste0(original_sj, "_X_", group_Pj)) := rowMeans(across(.cols = (any_of(sjPgrp_j)))))
    # } 
    
    # add the individual sjPgrp variables in the specification 
    regressors <- c(regressors, sjPgrp)
    
    # remove the group variable
    regressors <- regressors[regressors != paste0(original_sj, "_X_", aggregate_K)] 
    
    
    ## This is a particular case that is not handled by the code above. 
    # in this case, we recompute the interaction with soy, including all 3 soy commodities
    if(original_sj == "Soybean" & aggregate_K != "soy"){
      sjPsoy <- c()
      for(soy_commo in unlist(maplist[[original_sj]][["soy"]])){
        varname <- paste0(original_sj, "_X_", soy_commo)
        sjPsoy <- c(sjPsoy, varname)
        Pks <- paste0("ln_", soy_commo, price_info)
        d <- mutate(d,
                    !!as.symbol(varname) := (!!as.symbol(sj)) * (!!as.symbol(Pks)) )# don't log Pk here, it as been pre logged above, if soy is specified as a group of interest for sj
      }
      rm(soy_commo, Pks, varname)
      
      d <- mutate(d, 
                  Soybean_X_soy := rowMeans(across(.cols = (any_of(sjPsoy)))))
      
      # remove the single variable sjPj from specification
      regressors <- regressors[regressors != "Soybean_X_Soybean"]
    }
  }  
  ### OLD, NOT WORKING WAY TO HANDLE AGGREGATION ####  
  #   # In any case, we may need these objects
  #   sjPsoy <- c()
  #   for(soy_commo in unlist(maplist[[original_sj]][["soy"]])){
  #     varname <- paste0(original_sj, "_X_", soy_commo)
  #     sjPsoy <- c(sjPsoy, varname)
  #     Pks <- paste0("ln_", soy_commo, price_info)
  #     d <- mutate(d,
  #                 !!as.symbol(varname) := (!!as.symbol(sj)) * (!!as.symbol(Pks)) )# don't log Pk here, it as been pre logged above, if soy is specified as a group of interest for sj
  #   }
  #   rm(soy_commo)
  #     
  #   if(aggregate_K=="soy"){
  #     # if we aggregate of soy commodities, and Pj is Soybean, then we don"t want it as a stand alone regressor, we want it in soy, just like Soybean_meal and Soybean_oil
  #     # if sj is Soybean, we remove sjPj because we want it to be interacted with soy now
  #     # if we aggregate over soy commodities, but not for sj = soybean, then we still need the sjPsoy variables, the only difference is we do not need to remove sjPj
  #     # it is not necessary to condition on sj == Soybean
  #     # BUT IT IS IMPORTANT that is is run before the following line of code, otherwise Soybean_X_Soybean does not appear at all. 
  #     regressors <- regressors[regressors != "Soybean_X_Soybean"]
  #     
  #     # the sjPsoy variables are added as regressors 
  #     regressors <- c(regressors, sjPsoy)
  #     
  #     # sj_X_soy is removed below
  #   }
  #   if(aggregate_K != "soy"){
  #     # sjPsoy are averaged within obs. 
  #     # (equivalent to first averaging the log(prices) and then multiplying by sj)
  #     # necessary to redo it here if sj == soybean because then Soybean_X_soy computed above did not include Soybean_X_Soybean
  #     if(original_sj == "Soybean"){
  #       d <- mutate(d, 
  #                   Soybean_X_soy := rowMeans(across(.cols = (any_of(sjPsoy)))))
  #     }      
  #     # and remove any individual sjPsoy term (should not change anything that this is within the condition above or not, as long as aggregate procedure specifies Pk = Pj, because then there is no Oilpalm_X_Soybean in regressors)
  #     regressors <- regressors[!(regressors %in% sjPsoy)]
  #     
  #     # as above, remove the sjPj term so that it does not appear twice 
  #     regressors <- regressors[regressors != sjPk]
  #     
  #     # and add to regressors all the individual terms we want to aggregate
  #     regressors <- c(regressors, paste0(original_sj, "_X_", original_Pks))
  #     
  #   }
  # 
  #   # remove the group variable
  #   regressors <- regressors[regressors != paste0(original_sj, "_X_", aggregate_K)]  
  
  
  
  
  
  
  ### Controls - mechanisms ####
  # it's important that this is not conditioned on anything so these objects exist
  controls <- c()
  # gamma_controls <- c() 
  # delta_controls <- c() 
  
  # IT IS NORMAL THAT ABS. EXPOSURE CONTROLS ARE INTERACTED WITH 2fya AND 1pya: BOTH ARE AVERAGES OF TWO YEARS. 
  # 1pya is the avg of current and lag1 values, which are grouped together because they reflect the same kind of mechanism. 
  if(control_all_absolute_rfs){
    all_abs <- exp_matmap[, "Crops"][!(exp_matmap[, "Crops"] %in% original_exposure_rfs)]
    
    # add the share of pasture in grid cell in 2000 as an exposure to pastures (if the exposure of interest is not pasture)
    if(control_pasture & !grepl("Fodder", original_exposure_rfs)){
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
  if((rfs_reg) & group_exposure_rfs){
    sj <- "grp_exp"
  }
  
  # d <- mutate(d, sj_year := paste0(!!as.symbol(sj),"_",year))
  
  # if(fe == "country_year" | fe == "year"){
  #   controls <- sj
  # }
  
  # add remainging forest cover as a control
  if(remaining & !offset){
    controls <- c(controls, "remaining_fc")
    # gamma_controls <- c(gamma_controls, "remaining_fc")
    # delta_controls <- c(delta_controls, "remaining_fc")
  }
  
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
      d <- mutate(d, !!as.symbol(varname) := !!as.symbol(eaear_exp_rfs) * (year-2000))
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
  
  # if we don't control for indirect interactions with sj, then we still add the interactions with the extra prices in the specification
  # if(!(control_interact_sj) & control_direct){
  #   meats_ctrl <- c()
  #   for(original_mc in meats[meats != original_Pj & meats != original_Pk]){ # remove original_Pj for cases when sj is fodder, and Fodder_X_Beef is already in direct_controls
  #     varname <- paste0(original_sj, "_X_", original_mc)
  #     meats_ctrl <- c(meats_ctrl, varname)
  #     mc <- paste0("ln_", original_mc, price_info)
  #     d <- mutate(d,
  #                 !!as.symbol(varname) := !!as.symbol(sj) * (!!as.symbol(mc)) )
  #   }
  #   
  #   d <- mutate(d, sj_X_meats = rowMeans(across(.cols = (any_of(meats_ctrl)))))
  #   controls <- c(controls, "sj_X_meats")
  #   
  #   if(original_Pk != "Crude_oil"){
  #     co <- paste0("ln_", "Crude_oil", price_info)
  #     d <- mutate(d,
  #                 sj_X_Crude_oil := !!as.symbol(sj) * (!!as.symbol(co)) )
  #     controls <- c(controls, "sj_X_Crude_oil")
  #   }
  #   
  #   if(original_Pk != "Fertilizer"){
  #     co <- paste0("ln_", "Fertilizer", price_info)
  #     d <- mutate(d,
  #                 sj_X_Fertilizer := !!as.symbol(sj) * (!!as.symbol(co)) )
  #     controls <- c(controls, "sj_X_Fertilizer")
  #   }
  # }
  
  
  
  
  # if we are to group these controls in one single variable
  # if(define_an_option_name){
  #   sj_interactions <- c()
  #   for(at in all_treatments[all_treatments!=treatments]){
  #     varname <- paste0(original_sj, "_X_", original_all_treatments[match(at, all_treatments)])
  #     sj_interactions <- c(sj_interactions, varname)
  #     
  #   }
  #   rm(varname, at)
  #   
  #   d <- mutate(d, one_ctrl = rowSums(across(.cols = (any_of(sj_interactions)))))
  #   controls <- c(controls, "one_ctrl")
  # }
  
  
  # add indirect interactions with sj
  if(control_interact_Pk){
    
    # code like this allows to control which controls to drop
    s_ref <- paste0(reference_crop, standardization)
    # handle cases where s_ref has several elements
    d <- dplyr::mutate(d, s_ref = rowSums(across(.cols = (any_of(s_ref)))))
    
    # for(si in all_exposures[!(all_exposures %in% c(sj, s_ref))]){
    #   varname <- paste0(original_all_exposures[match(si, all_exposures)], "_X_", original_treatments)
    #   controls <- c(controls, varname)
    #   # Log prices so their variations are comparable
    #   d <- mutate(d,
    #               !!as.symbol(varname) := (!!as.symbol(si)) * log( (!!as.symbol(treatments)) ))#+1
    # }
    
    # for(si in all_exposures_but_k[all_exposures_but_k != sj]){ # do not put the regressor a second time
    #   varname <- paste0(original_all_exposures_but_k[match(si, all_exposures_but_k)], "_X_", original_treatments)
    #   controls <- c(controls, varname)
    #   # Log prices so their variations are comparable
    #   d <- mutate(d,
    #               !!as.symbol(varname) := (!!as.symbol(si)) * log( (!!as.symbol(treatments)) ))#+1
    # }
    # rm(varname, si)
    
    # Code below to control for other indirect interactions with Pk through controlling only for (1-sj-s_ref)Pk  
    # (which is supposed to be equivalent to main approach)
    
    varname <- "one_ctrl"
    controls <- c(controls, varname)
    d <- mutate(d,
                !!as.symbol(varname) := (1 - !!as.symbol(sj) - s_ref ) * (!!as.symbol(treatments)) )#
    rm(varname)
    
    # Code below to control for other indirect interactions with Pk through controlling only for skPk 
    # (which is bad, because it is equivalent to (1-sk)Pk and thus throws sjPk once again)
    # original_sk <- mapmat_si[mapmat_si[,"Prices"]==original_treatments, "Crops"]
    # sk <- paste0(original_sk, standardization)
    # varname <- paste0(original_sk, "_X_", original_treatments)
    # controls <- c(controls, varname)
    # # Log prices so their variations are comparable
    # d <- mutate(d,
    #             !!as.symbol(varname) := (!!as.symbol(sk)) * log( (!!as.symbol(treatments)) ))#+1
    # rm(original_sk, sk, varname)
    
    
  }
  
  
  # add direct interactions
  if(control_direct){
    # this is to control for direct interactions for main drivers of tropical deforestation
    # main_drivers <- c("Fodder", "Soybean", "Oilpalm", "Cocoa", "Coffee", "Rubber")
    # for(original_sm in main_drivers[main_drivers != original_sj]){
    #   original_Pm <- mapmat_si[mapmat_si[,"Crops"]==original_sm, "Prices"]
    #   Pm <- paste0("ln_", original_Pm, price_info)
    #   sm <- paste0(original_sm, standardization)
    #   varname <- paste0(original_sm, "_X_", original_Pm)
    #   controls <- c(controls, varname)
    #   d <- mutate(d,
    #               !!as.symbol(varname) := (!!as.symbol(sm)) * (!!as.symbol(Pm)) )# if group_price (the only way currently), prices are already logged. 
    # }
    
    # original_sk <- mapmat_si[mapmat_si[,"Prices"]==original_Pk, "Crops"]
    # sk <- paste0(original_sk, standardization)
    # skPk <- paste0(original_sk, "_X_", original_Pk)
    # d <- mutate(d,
    #             !!as.symbol(skPk) := (!!as.symbol(sk) * !!as.symbol(Pk)))
    # regressors <- c(regressors, skPk)
    
    # this is to control for direct interactions from all treatments, grouped in one single variable 
    # ||| original_treatments != original_Pk &  original_treatments != original_Pj &   |||
    # currently: remove Pj because we control for sjPj individually. Par contre, skPk does not make much sense given sj>0 condition, so we just include it in dir_ctrl_grp
    
    Pms <- all_treatments[original_all_treatments != original_Pj & original_all_treatments != original_Pk]
    # add skPk in direct controls if we want to control for it, include it in combined direct controls, and if it has a suitability index to be matched with
    if(control_skPk & incl_skPk  & original_Pk %in% mapmat_si[,"Prices"]){
      Pms <- c(Pk, Pms)
    }
    if(control_sjPj & incl_sjPj & (Pj != Pk)){
      Pms <- c(Pj, Pms)        
    }
    
    # if Pk is one of the soy commodities, then we don't want soy_index to be interacted with suitability for soy in direct controls 
    # (be we in sj = Soybean or not)
    if(original_Pk == "Soybean_meal" | original_Pk == "Soybean_oil"){
      Pms <- Pms[Pms != paste0("ln_", "Soy_index", price_info)]
    }   
    
    #   if(all_but_k){# in this case, let sjPJ be in the grouped direct controls
    #     Pms <- all_treatments[original_all_treatments != original_Pk]# & original_treatments %in% mapmat_si[,"Prices"]
    #     # but it's then also necessary to rename it before it is produced and called the same way in the loop below
    #     noskPk <- paste0("no",original_sk, "_X_", original_Pk)
    #     names(d)[names(d)==paste0(original_sj, "_X_", original_Pk)] <- noskPk
    #     regressors <- noskPk
    #     
    #   }else{ # but otherwise, it is just a test with sjPj being specified individually so remove it from the group of direct controls
    #     Pms <- all_treatments[original_all_treatments != original_Pj & original_all_treatments != original_Pk]# & original_treatments %in% mapmat_si[,"Prices"]
    #   }
    
    direct_controls <- c()
    for(Pm in Pms){ 
      # find the suitability index that matches Pm
      original_Pm <- original_all_treatments[match(Pm, all_treatments)]
      original_sm <- mapmat_si[mapmat_si[,"Prices"]==original_Pm, "Crops"]
      sm <- paste0(original_sm, standardization)
      varname <- paste0(original_sm, "_X_", original_Pm)
      direct_controls <- c(direct_controls, varname)
      d <- mutate(d,
                  !!as.symbol(varname) := (!!as.symbol(sm)) * (!!as.symbol(Pm)) )
    }
    # if Pk is one of the soy commodities, then we didn't include Soybean_X_Soy_index in direct controls, 
    # but we want to include the other soy commodity (neither has been included in control_interact_sj) - IF skPk IS INCLUDED IN DIRECT CONTROLS 
    # AND IF sj is not Soybean, bc if it is, then Soybean_X_othersoycommodity is sjPj so it is already in the regressors  
    # need to do this separately to handle matching between price and suitability
    if((original_Pk == "Soybean_meal" | original_Pk == "Soybean_oil") & control_skPk & incl_skPk & original_sj != "Soybean"){
      varname <- paste0("Soybean_X_", original_Pk)
      direct_controls <- c(direct_controls, varname)
      d <- mutate(d,
                  !!as.symbol(varname) := (!!as.symbol(sm)) * (!!as.symbol(Pk)) )
    }  
    
    # make sure there is no variable in direct_controls that is already in the regressors or the controls
    direct_controls <- direct_controls[!(direct_controls %in% c(regressors, controls))]
    
    d <- mutate(d, dir_ctrl_grp = rowSums(across(.cols = (any_of(direct_controls)))))
    controls <- c(controls, "dir_ctrl_grp")
    
    rm(varname, original_Pm, original_sm, sm, Pm)
    
    # We don't do this anymore because sjPj is not included in direct_controls, but add as an individual covariate.  
    # # just so that the objects exist, but they won"t be used
    # gamma_regressors <- regressors
    # gamma_controls <- controls    
    # # if we are in beta procedure, estimate a specification without direct control sjPj 
    # if(estimated_effect == "beta"){
    #   direct_controls_but_sjPj <- direct_controls[direct_controls != paste0(original_sj, "_X_", original_Pj)]
    #   d <- mutate(d, dir_ctrl_grp_but_sjPj = rowSums(across(.cols = (any_of(direct_controls_but_sjPj)))))
    #   gamma_controls <- c(controls[controls != "dir_ctrl_grp"], "dir_ctrl_grp_but_sjPj")
    # }
    
  }
  
  #### MODEL SPECIFICATION FORMULAE - alpha model ####
  
  if(estimated_effect == "beta"){
    if(original_sj == "Soybean"){
      gamma_regressors <- regressors[regressors != "Soybean_X_soy"]
    }else{
      gamma_regressors <- regressors[regressors != paste0(original_sj, "_X_", original_Pj)]
    }
    gamma_controls <- controls
    
    # MODEL SPECIFICATION FORMULAE - gamma and delta models
    gamma_model <- as.formula(paste0(outcome_variable,
                                     " ~ ",
                                     paste0(gamma_regressors, collapse = "+"),
                                     " + ",
                                     paste0(gamma_controls, collapse = "+"),
                                     " | ",
                                     fe)) 
  }
  
  
  # do NOT condition this to estimated_effect, as we will need it for beta estimation too
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
  
  # - are suitable to crop j 
  if(sjpos){
    d <- dplyr::filter(d, !!as.symbol(sj) > 0)  
  }
  
  # now that this is done, log rfs regressors: 
  # if(rfs_reg){
  #   d <- dplyr::mutate(d, !!as.symbol(regressors) := log(!!as.symbol(regressors)))  
  # }
  
  # # - are in study period 
  # if(start_year != 2011 | end_year != 2019){
  #   d <- dplyr::filter(d, year >= start_year)
  #   d <- dplyr::filter(d, year <= end_year)
  # }
  
  # - are in study area
  if(continent != "all"){
    d <- dplyr::filter(d, continent_name == continent)
  }
  
  # have remaining forest
  # d <- dplyr::filter(d, remaining_fc > 0)
  
  # remove units with no country name (in the case of all_drivers data set currently, because nearest_feature function has not been used in this case, see add_variables.R)
  d <- dplyr::filter(d, !is.na(country_name))
  
  if(estimated_effect == "alpha"){
    used_vars <- unique(c("grid_id", "year", "lat", "lon", "country_name", "country_year", "continent_name", #"remaining_fc", "accu_defo_since2k", # "sj_year",
                          original_exposure_rfs, original_rfs_treatments, rfs_treatments, outcome_variable, regressors, controls))
  }else{
    used_vars <- unique(c("grid_id", "year", "lat", "lon", "country_name", "country_year", "continent_name", # "remaining_fc", "accu_defo_since2k", # "sj_year",
                          original_exposure_rfs, original_rfs_treatments, rfs_treatments, outcome_variable, regressors, controls, gamma_controls))
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
  
  # experience deforestation at least once (automatically removed if distribution is quasipoisson, but not if it's gaussian. 
  # Though we want identical samples to compare distributional assumptions)
  if(distribution == "gaussian"){
    # # - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove
    obstormv <- obs2remove(fml = as.formula(paste0(outcome_variable, " ~ ", fe)),
                           data = d,
                           family = "poisson")
    # this is necessary to handle cases when there is no obs. to remove
    if(length(obstormv)>0){
      d <- d[-obstormv,]
    } else {
      d <- d
    }
    rm(obstormv)
  }
  
  if(fe != "grid_id + country_year"){
    # # - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove
    obstormv <- obs2remove(fml = as.formula(paste0(outcome_variable, " ~ grid_id + country_year")),
                           data = d,
                           family = "poisson")
    # this is necessary to handle cases when there is no obs. to remove
    if(length(obstormv)>0){
      d <- d[-obstormv,]
    } else {
      d <- d
    }
    rm(obstormv)
  }
  
  # if we don't remove obs that have no remaining forest, and if we don't remove always 0 outcome obs. in previous step, d IS A BALANCED PANEL !!
  # this is a matter only if we do beta process, with cluster bootstrap
  # if(!nrow(d) == length(unique(d$grid_id))*length(unique(d$year))){
  #   warning("data is not balanced")
  # }
  d_clean <- d
  rm(d)
  
  
  ### REGRESSIONS
  
  # Store only information necessary, in a dataframe. otherwise the output of fixest estimation is large and we can't collect too many at the same time (over loops)  
  # either there are several elemnts in regressors, and then we want to aggregate them, or there is only one. 
  # In both cases, we are interested in a one-line output
  
  # handle SE computation flexibly within feglm now, through argument vcov
  
  # in case of unit fe, adjust also the vcov argument, because two way wont work, but we might still want to cluster on country year
  if(fe == "grid_id" & se == "twoway"){
    se <- as.formula("~ grid_id + country_year")
  }
  if(fe == "country_year" & se == "twoway"){
    se <- as.formula("~ grid_id + country_year")
  }
  if(se == "exposure2ways"){
    se <- as.formula(paste0("~ ", original_exposure_rfs, " + country_year"))
  }
  
  
  if(estimated_effect == "alpha"){
    if(offset){
      offset_fml <- ~log(remaining_fc)
      alpha_reg_res <- fixest::feglm(alpha_model,
                                     data = d_clean, 
                                     family = distribution,# "gaussian",#  # "poisson" ,
                                     offset = offset_fml,
                                     vcov = se,
                                     # this is just to get the same p value by recomputing by hand below. 
                                     # see https://cran.r-project.org/web/packages/fixest/vignettes/standard_errors.html
                                     # ssc = ssc(cluster.df = "conventional", t.df = "conventional"),
                                     # glm.iter = 25,
                                     #fixef.iter = 100000,
                                     nthreads = 3,
                                     glm.iter = glm_iter,
                                     notes = TRUE, 
                                     verbose = 4)  
      
    }else{
      alpha_reg_res <- fixest::feglm(alpha_model,
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
                                     verbose = 4) 
    }    
    # alpha_reg_res <- summary(alpha_reg_res, vcov = as.formula("twoway ~ grid_id + country_year"))
    # keep only the coeff estimate for the sjPk term of interest
    if(Pj != Pk){
      df_res <- summary(alpha_reg_res)$coeftable#[paste0(original_sj, "_X_", original_Pk), ]
    }
    if(Pj==Pk){
      df_res <- summary(alpha_reg_res)$coeftable
    }
    
    fixest_df <- degrees_freedom(alpha_reg_res, type = "t")
    ## MAKE AGGREGATE RESULTS 
    # In this case, we are interested in LEADS only
    if(aggr_dyn == "leads"){ 
      # deprecated: in this case, there ARE lags in the specification, but we are not interested in them
      
      df_res <- rbind(rep(NA, ncol(df_res)), df_res)
      
      # ORDER MATTERS
      aggr_names <- paste0(original_exposure_rfs, c("_X_aggrleads"))
      row.names(df_res)[1] <- aggr_names
      # regressors of interest
      # base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
      # Contemporaneous value is NOT in leads
      lead_roi <- c(grep(pattern = paste0(base_reg_name,"_lead"), 
                         regressors, value = TRUE))
      
      df_res[aggr_names[1],"Estimate"] <- alpha_reg_res$coefficients[lead_roi] %>% sum()
      
      # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
      # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
      df_res[aggr_names[1],"Std. Error"] <- alpha_reg_res$cov.scaled[lead_roi, lead_roi] %>% as.matrix() %>% sum() %>% sqrt()
      
      df_res[aggr_names[1],"t value"]  <- (df_res[aggr_names[1],"Estimate"] - 0)/(df_res[aggr_names[1],"Std. Error"])
      
      # use t distribution with degrees of freedom equal to that used by fixest, i.e. after two way cluster adjustment.
      # does not make a significant difference given sample size
      df_res[aggr_names[1],"Pr(>|t|)"]  <- (2*pt(abs(df_res[aggr_names[1],"t value"]), 
                                                 lower.tail = FALSE, 
                                                 df = fixest_df)) 
    }
    # In this case, we are interested in LAGS only
    if(aggr_dyn == "lags"){
      df_res <- rbind(rep(NA, ncol(df_res)), df_res)
      
      # ORDER MATTERS
      aggr_names <- paste0(original_exposure_rfs, c("_X_aggrlags"))
      row.names(df_res)[1] <- aggr_names
      # regressors of interest
      base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
      # Contemporaneous value is in lags
      lag_roi <- c(base_reg_name, 
                   grep(pattern = paste0(base_reg_name,"_lag"), 
                        regressors, value = TRUE))
      
      df_res[aggr_names[1],"Estimate"] <- alpha_reg_res$coefficients[lag_roi] %>% sum()
      
      # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
      # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
      df_res[aggr_names[1],"Std. Error"] <- alpha_reg_res$cov.scaled[lag_roi, lag_roi] %>% as.matrix() %>% sum() %>% sqrt()
      
      df_res[aggr_names[1],"t value"]  <- (df_res[aggr_names[1],"Estimate"] - 0)/(df_res[aggr_names[1],"Std. Error"])
      
      # use t distribution with degrees of freedom equal to that used by fixest, i.e. after two way cluster adjustment.
      # does not make a significant difference given sample size
      df_res[aggr_names[1],"Pr(>|t|)"]  <- (2*pt(abs(df_res[aggr_names[1],"t value"]), 
                                                 lower.tail = FALSE, 
                                                 df = fixest_df)) 
    }
    
    
    if(aggr_dyn & rfs_lead > 0 & rfs_lag > 0 & rfs_fya == 0 & rfs_pya == 0){
      # In this case, we are interested in LEAD AND LAG effects, aggregated separately, and all together.
      df_res <- rbind(rep(NA, ncol(df_res)), rep(NA, ncol(df_res)), df_res)
      
      # ORDER MATTERS
      aggr_names <- paste0(original_exposure_rfs, c("_X_aggrleads", "_X_aggrlags"))
      row.names(df_res)[1:2] <- aggr_names
      # regressors of interest
      base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
      # Contemporaneous value is NOT in leads
      lead_roi <- c(grep(pattern = paste0(base_reg_name,"_lead"), 
                         regressors, value = TRUE))
      # Contemporaneous value is in lags
      lag_roi <- c(base_reg_name, 
                   grep(pattern = paste0(base_reg_name,"_lag"), 
                        regressors, value = TRUE))
      
      df_res[aggr_names[1],"Estimate"] <- alpha_reg_res$coefficients[lead_roi] %>% sum()
      df_res[aggr_names[2],"Estimate"] <- alpha_reg_res$coefficients[lag_roi] %>% sum()
      
      # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
      # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
      df_res[aggr_names[1],"Std. Error"] <- alpha_reg_res$cov.scaled[lead_roi, lead_roi] %>% as.matrix() %>% sum() %>% sqrt()
      df_res[aggr_names[2],"Std. Error"] <- alpha_reg_res$cov.scaled[lag_roi, lag_roi] %>% as.matrix() %>% sum() %>% sqrt()
      
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
      aggr_names <- paste0(original_exposure_rfs, c("_X_aggrall"))
      row.names(df_res)[1] <- aggr_names
      # regressors of interest
      base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
      
      # Contemporaneous, lead, and lag values
      all_roi <- c(base_reg_name, 
                   grep(pattern = paste0(base_reg_name,"_lead"), 
                        regressors, value = TRUE),
                   grep(pattern = paste0(base_reg_name,"_lag"), 
                        regressors, value = TRUE))
      
      df_res[aggr_names[1],"Estimate"] <- alpha_reg_res$coefficients[all_roi] %>% sum()
      
      # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
      # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
      df_res[aggr_names[1],"Std. Error"] <- alpha_reg_res$cov.scaled[all_roi, all_roi] %>% as.matrix() %>% sum() %>% sqrt()
      
      df_res[aggr_names[1],"t value"]  <- (df_res[aggr_names[1],"Estimate"] - 0)/(df_res[aggr_names[1],"Std. Error"])
      
      # use t distribution with degrees of freedom equal to that used by fixest, i.e. after two way cluster adjustment.
      # does not make a significant difference given sample size
      df_res[aggr_names[1],"Pr(>|t|)"]  <- (2*pt(abs(df_res[aggr_names[1],"t value"]), 
                                                 lower.tail = FALSE, 
                                                 df = fixest_df)) 
    }
    if(aggr_dyn & rfs_lead == 0 & rfs_lag == 0 & rfs_fya > 0 & rfs_pya > 0){
      # In this case, we are interested in AVERAGE LEAD AND LAG effects, separately and aggregated
      df_res <- rbind(rep(NA, ncol(df_res)), df_res)
      
      # ORDER MATTERS
      aggr_names <- paste0(original_exposure_rfs, c("_X_aggrall"))
      row.names(df_res)[1] <- aggr_names
      # regressors of interest
      base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
      
      # Contemporaneous, lead, and lag values
      all_roi <- c(grep(pattern = paste0(base_reg_name,"_",rfs_fya,"fya"), 
                        regressors, value = TRUE),
                   grep(pattern = paste0(base_reg_name,"_",rfs_pya,"pya"), 
                        regressors, value = TRUE))
      
      df_res[aggr_names[1],"Estimate"] <- alpha_reg_res$coefficients[all_roi] %>% sum()
      
      # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
      # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
      df_res[aggr_names[1],"Std. Error"] <- alpha_reg_res$cov.scaled[all_roi, all_roi] %>% as.matrix() %>% sum() %>% sqrt()
      
      df_res[aggr_names[1],"t value"]  <- (df_res[aggr_names[1],"Estimate"] - 0)/(df_res[aggr_names[1],"Std. Error"])
      
      # use t distribution with degrees of freedom equal to that used by fixest, i.e. after two way cluster adjustment.
      # does not make a significant difference given sample size
      df_res[aggr_names[1],"Pr(>|t|)"]  <- (2*pt(abs(df_res[aggr_names[1],"t value"]), 
                                                 lower.tail = FALSE, 
                                                 df = fixest_df)) 
    }
    if(aggr_dyn & rfs_lead > 0 & rfs_lag == 0 & rfs_fya == 0 & rfs_pya > 0){
      # In this case, we are interested in aggregated LEAD effects, and all together.
      
      # first aggregate leads
      df_res <- rbind(rep(NA, ncol(df_res)), df_res)
      
      # ORDER MATTERS
      aggr_names <- paste0(original_exposure_rfs, c("_X_aggrleads"))
      row.names(df_res)[1] <- aggr_names
      # regressors of interest
      base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
      # Contemporaneous value is NOT in leads
      lead_roi <- c(grep(pattern = paste0(base_reg_name,"_lead"), 
                         regressors, value = TRUE))
      
      df_res[aggr_names[1],"Estimate"] <- alpha_reg_res$coefficients[lead_roi] %>% sum()
      
      # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
      # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
      df_res[aggr_names[1],"Std. Error"] <- alpha_reg_res$cov.scaled[lead_roi, lead_roi] %>% as.matrix() %>% sum() %>% sqrt()
      
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
      aggr_names <- paste0(original_exposure_rfs, c("_X_aggrall"))
      row.names(df_res)[1] <- aggr_names
      # regressors of interest
      base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
      
      # Contemporaneous, lead, and lag values
      all_roi <- c(grep(pattern = paste0(base_reg_name,"_lead"), 
                        regressors, value = TRUE),
                   grep(pattern = paste0(base_reg_name,"_",rfs_pya,"pya"), 
                        regressors, value = TRUE))
      
      df_res[aggr_names[1],"Estimate"] <- alpha_reg_res$coefficients[all_roi] %>% sum()
      
      # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
      # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
      df_res[aggr_names[1],"Std. Error"] <- alpha_reg_res$cov.scaled[all_roi, all_roi] %>% as.matrix() %>% sum() %>% sqrt()
      
      df_res[aggr_names[1],"t value"]  <- (df_res[aggr_names[1],"Estimate"] - 0)/(df_res[aggr_names[1],"Std. Error"])
      
      # use t distribution with degrees of freedom equal to that used by fixest, i.e. after two way cluster adjustment.
      # does not make a significant difference given sample size
      df_res[aggr_names[1],"Pr(>|t|)"]  <- (2*pt(abs(df_res[aggr_names[1],"t value"]), 
                                                 lower.tail = FALSE, 
                                                 df = fixest_df)) 
    }
    
    if(aggregate_K != ""){
      df_res <- rbind(df_res, rep(NA, ncol(df_res)))
      
      row.names(df_res)[nrow(df_res)] <- "aggr_K"
      
      df_res["aggr_K","Estimate"] <- alpha_reg_res$coefficients[sjPks] %>% sum()
      # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
      # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
      df_res["aggr_K","Std. Error"] <- alpha_reg_res$cov.scaled[sjPks, sjPks] %>% as.matrix() %>% sum() %>% sqrt()
      
      df_res["aggr_K","t value"]  <- (df_res["aggr_K","Estimate"] - 0)/(df_res["aggr_K","Std. Error"])
      
      df_res["aggr_K","Pr(>|t|)"]  <- (2*pt(abs(df_res["aggr_K","t value"]), 
                                            lower.tail = FALSE, 
                                            df = fixest_df)) 
      
    }
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
    # select_coefs <- select_coefs[select_coefs != paste0(mapmat_si[mapmat_si[,"Prices"]==original_treatments, "Crops"], "_X_", original_treatments)]
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
    
    df_res <- data.frame(estimate = NA, std.error = NA, t.statistic = NA, p.value = NA, observations = NA, inference = NA)  
    
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
                                     verbose = 0)  
      
      alpha_reg_res <- fixest::feols(alpha_model,
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
                                     verbose = 0)  
      
      if(aggregate_K == ""){
        gamma_coeff <- gamma_reg_res$coefficients[sjPk]
        #gamma_aggr_se <- gamma_reg_res$cov.scaled[regressors,regressors] %>% as.matrix() %>% sum() %>% sqrt()
        
        alpha_coeff <- alpha_reg_res$coefficients[sjPk]
        #delta_aggr_se <- delta_reg_res$cov.scaled[regressors,regressors] %>% as.matrix() %>% sum() %>% sqrt()
        
        # # this is our coefficient of interest
        beta <- gamma_coeff - alpha_coeff
      }else{
        gamma_aggr_K <- gamma_reg_res$coefficients[sjPks] %>% sum()
        
        alpha_aggr_K <- alpha_reg_res$coefficients[sjPks] %>% sum()
        
        beta <- gamma_aggr_K - alpha_aggr_K
        
      }     
      
      return(beta)
      
      
      # # its SE is conservatively estimated as if delta and gamma coefficients had a null covariance.  
      # df_res[,"std.error"] <- gamma_aggr_se + delta_aggr_se 
      # 
      # df_res[,"observations"] <- gamma_reg_res$nobs
    }
    
    # methodology comes from Cameron and Miller (2015), and:
    # https://stats.stackexchange.com/questions/202916/cluster-boostrap-with-unequally-sized-clusters/202924#202924
    start <- Sys.time()
    # for two way clustered SE, we need to do 3 bootstrap estimations: by first-way clustering, 2nd way clustering, and without clustering (i.e. by 1st way x 2nd way)
    
    # store these 3 SE
    SEs <- c(grid_id = NA, country_year = NA, grid_id_X_country_year = NA)
    
    for(boot_cluster in c("grid_id", "country_year")){
      # list of parameters related to the dataset, and the clustering variable
      # names and numbers of clusters of size s
      par_list <- list(cluster_variable = boot_cluster, 
                       cluster_names = unique(d_clean[,boot_cluster]),
                       number_clusters = length(unique(d_clean[,boot_cluster])))
      
      
      set.seed(145)
      bootstraped <- boot(data = d_clean, 
                          statistic = make_estimate, 
                          ran.gen = ran.gen_cluster_blc,
                          mle = par_list,
                          sim = "parametric",
                          # parallel = "multicore",
                          # ncpus = detectCores() - 1,
                          R = 500)
      
      
      
      SEs[boot_cluster] <- sd(bootstraped$t)
      
    }
    
    set.seed(145)
    bootstraped <- boot(data = d_clean, 
                        statistic = make_estimate, 
                        ran.gen = ran.gen_blc,
                        mle = list(size = nrow(d_clean)),
                        sim = "parametric",
                        # parallel = "multicore",
                        # ncpus = detectCores() - 1,
                        R = 500)
    SEs["grid_id_X_country_year"] <- sd(bootstraped$t) 
    
    df_res[,"estimate"] <- bootstraped$t0 
    df_res[,"std.error"] <- SEs["grid_id"] + SEs["country_year"] - SEs["grid_id_X_country_year"]
    df_res[,"observations"] <- nrow(d_clean)
    
    end <- Sys.time()
    end - start
    
    # compute t.stat and p-value of the aggregated estimate (yields similar quantities as in alpha_reg_res if only one coeff was aggregated)
    df_res <- dplyr::mutate(df_res, t.statistic = (estimate - 0)/(std.error))
    
    df_res <- dplyr::mutate(df_res, p.value =  (2*pt(abs(t.statistic), 
                                                     lower.tail = FALSE, 
                                                     df = fixest_df)) )
    
  }
  
  # take data set as exactly used in estimation
  if(length(alpha_reg_res$obs_selection) > 0){
    d_clean <- d_clean[alpha_reg_res$obs_selection[[1]], ]
  }
  
  if(rfs_rando=="between"){
    # in this case, we resample EXPOSURES (without replacement), and keep rfs history intact. i.e. we randomize the between units variation
    exp_dfcs <- d_clean[!duplicated(d_clean$grid_id), c("grid_id", exposure_rfs)]
    nrow_dfcs <- nrow(exp_dfcs)
    
    # remove the exposure column in the original order from data (and its interactions with treatments)
    d_clean <- dplyr::select(d_clean, !matches(exposure_rfs))
    
    # just to gain some time
    sym_exposure_rfs <- as.symbol(exposure_rfs)
    
    # store permutation results there
    perm_res_list <- list()
    
    set.seed(145)
    #start_time <- Sys.time()
    for(i in 1:2000){
      # sample and merge back the permuted values  
      exp_dfcs[,exposure_rfs] <- sample(exp_dfcs[,exposure_rfs], size = nrow_dfcs, replace = FALSE)
      
      d_clean_i <- left_join(d_clean, exp_dfcs, by = "grid_id")
      
      # build regressors and exposure trend
      for(rfs_var in rfs_treatments){
        d_clean_i <- mutate(d_clean_i, 
                            !!as.symbol(paste0(exposure_rfs, "_X_", rfs_var)) := !!sym_exposure_rfs * !!as.symbol(rfs_var) ) 
      }
      d_clean_i <- mutate(d_clean_i, 
                          !!as.symbol(paste0(exposure_rfs, "_trend")) := !!sym_exposure_rfs * year)
      
      # rerun regression
      perm_i_res <- fixest::feglm(alpha_model,
                                  data = d_clean_i,
                                  family = distribution,
                                  vcov = as.formula("~ grid_id"), # do not cluster on country here, as exposures have been randomized within countries too  
                                  nthreads = 0,
                                  fixef.rm = "perfect")
      
      # aggregate coefficients on leads and lags 
      
      df_i_est <- data.frame(Estimate = c(sum(perm_i_res$coefficients[lead_roi]), 
                                          sum(perm_i_res$coefficients[lag_roi]) ),
                             row.names = aggr_names)
      
      # compute their t statistics
      # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
      # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
      df_i_est[aggr_names[1],"Std. Error"] <- perm_i_res$cov.scaled[lead_roi, lead_roi] %>% as.matrix() %>% sum() %>% sqrt()
      df_i_est[aggr_names[2],"Std. Error"] <- perm_i_res$cov.scaled[lag_roi, lag_roi] %>% as.matrix() %>% sum() %>% sqrt()
      
      df_i_est[aggr_names[1],"t value"]  <- (df_i_est[aggr_names[1],"Estimate"] - 0)/(df_i_est[aggr_names[1],"Std. Error"])
      df_i_est[aggr_names[2],"t value"]  <- (df_i_est[aggr_names[2],"Estimate"] - 0)/(df_i_est[aggr_names[2],"Std. Error"])
      
      perm_res_list[[i]] <- df_i_est
    }
    #end_time <- Sys.time()
    #end_time - start_time
    
    # 10 replications take 5 minutes 
    # so for 1000 replications, it would take 500/60 = 8.3 hours
    
    # without the regression, the loop takes 3s. which, replicated 1000 times mores, makes less than an hour (0.83). 
    # to apply to all 20 crops in three continents this would take between 50 and 60h 
    
    # on server 10 replications is 2.16 minutes with all available threads used. 
    # 216/60 # 4 hours for 1000 replications, 40h for 10000
    
    # Let's return both the coefficients and the t-statistics (and randomization inference p values of them)
    # so this returns estimates for every permutation in columns, and for each regressor in row 
    all_coeffs <- bind_cols(lapply(perm_res_list, FUN = function(x){dplyr::select(x, Estimate)}), 
                            .name_repair = "minimal")
    colnames(all_coeffs) <- NULL
    
    # t-statistics
    all_tstats <- bind_cols(lapply(perm_res_list, FUN = function(x){dplyr::select(x, `t value`)}), 
                            .name_repair = "minimal")
    colnames(all_tstats) <- NULL
    
    # randomization inference p values are equal to the share of a parameter estimated over the randomized distribution, that is higher than the original estimate 
    raninf_pval_coeffs <- c()
    raninf_pval_coeffs[1] <- sum( abs(all_coeffs[aggr_names[1], ]) > abs(df_res[aggr_names[1], "Estimate"]) ) / (i)
    raninf_pval_coeffs[2] <- sum( abs(all_coeffs[aggr_names[2], ]) > abs(df_res[aggr_names[2], "Estimate"]) ) / (i)
    names(raninf_pval_coeffs) <- aggr_names
    
    # this returns the share of absolute t ratios that are larger than the one we obtain with the actual data
    raninf_pval_tstats <- c()
    raninf_pval_tstats[1] <- sum( abs(all_tstats[aggr_names[1], ]) > abs(df_res[aggr_names[1], "t value"]) ) / (i)
    raninf_pval_tstats[2] <- sum( abs(all_tstats[aggr_names[2], ]) > abs(df_res[aggr_names[2], "t value"]) ) / (i)
    names(raninf_pval_tstats) <- aggr_names
    
  }
  
  if(rfs_rando=="within"){
    # in this case, we resample RFS TIME SERIES (without replacement), and keep spatial exposures intact. 
    # i.e. we randomize the within unit variation
    # it differs from the "between" procedure above because we have to resample the original RFS history, and then 
    # build the leads and lags
    
    rfs_ts <- prices[prices$year >=2008, c("year", original_rfs_treatments)] 
    
    # permute only unique values, 
    # and it differs because we sample only within years of the program, i.e. just permute the agenda
    rfs_agd <- rfs_ts[rfs_ts$year >=2008, ] # 2008 is the ealiest year needed in case rfs_lag is 3
    nrow_dfcs <- nrow(rfs_agd)
    avg_rfs <- mean(rfs_agd[,original_rfs_treatments])
    sd_rfs <- sd(rfs_agd[,original_rfs_treatments])
    
    # just to be sure that the datacombine operations work on an ordered dataset
    rfs_agd <- dplyr::arrange(rfs_agd, year)
    
    # remove the exposure columnS in the original order from data (and not its interactions with treatments, but those are replaced by mutate below)
    d_clean <- dplyr::select(d_clean, !any_of(rfs_treatments))
    
    # store permutation results there
    perm_res_list <- list()
    voi <- original_rfs_treatments
    
    set.seed(145)
    #start_time <- Sys.time()
    for(i in 1:5000){
      # resample RFS history
      rfs_ts[,voi] <- sample(unique(rfs_agd[,voi]), size = nrow_dfcs, replace = TRUE) # this allows other values than 15 to be repeatedly sampled
      #rfs_ts[,voi] <- sample(rfs_agd[,voi], size = nrow_dfcs, replace = FALSE) #this forces 15 to be repeatedly sampled
      #rfs_agd[,voi] <- rpois(n = nrow_dfcs, lambda = avg_rfs) 
      # normal distribution preferred because yields real values, not integers
      # round to 2 decimals to match degree of precision (and variation) of RFS
      # merge it with the full time series       
      #rfs_ts[,voi] <- rnorm(n = nrow_dfcs, mean = avg_rfs, sd = sd_rfs) %>% round(2)
      # %>% base::sort() don't sort, because the 
      
      # and rebuild leads and lags from it
      for(lead in c(1:(rfs_lead+1))){# +1 is just a trick to handle case when rfs_lead is 0
        rfs_ts <- DataCombine::slide(rfs_ts,
                                     Var = voi, 
                                     TimeVar = "year",
                                     NewVar = paste0(voi,"_lead",lead),
                                     slideBy = lead, 
                                     keepInvalid = FALSE)
        rfs_ts <- dplyr::arrange(rfs_ts, year)
      }  
      for(lag in c(1:rfs_lag)){
        rfs_ts <- DataCombine::slide(rfs_ts,
                                     Var = voi, 
                                     TimeVar = "year",
                                     NewVar = paste0(voi,"_lag",lag),
                                     slideBy = -lag, 
                                     keepInvalid = FALSE)
        rfs_ts <- dplyr::arrange(rfs_ts, year)
      }  
      
      # merge the new history with estimation dataset 
      d_clean_i <- left_join(d_clean, rfs_ts, by = "year")
      
      # build regressors
      for(rfs_var in rfs_treatments){
        d_clean_i <- mutate(d_clean_i, 
                            !!as.symbol(paste0(exposure_rfs, "_X_", rfs_var)) := !!as.symbol(rfs_var) * !!as.symbol(exposure_rfs) ) 
      }
      
      # rerun regression
      perm_i_res <- fixest::feglm(alpha_model,
                                  data = d_clean_i,
                                  family = distribution,
                                  vcov = se,
                                  # vcov =  as.formula("~ country_year"),
                                  nthreads = 0,
                                  fixef.rm = "perfect")
      
      if(length(perm_i_res$collin.var)==0){ # this condition forces passage to the next iteration if by chance the permutation causes perfect collinearity.
        # aggregate coefficients on leads and lags 
        df_i_est <- data.frame(Estimate = c(sum(perm_i_res$coefficients[all_roi])),#sum(perm_i_res$coefficients[lag_roi])
                               row.names = aggr_names)
        
        # compute their t statistics
        # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
        # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
        df_i_est[aggr_names[1],"Std. Error"] <- perm_i_res$cov.scaled[all_roi, all_roi] %>% as.matrix() %>% sum() %>% sqrt()
        #df_i_est[aggr_names[2],"Std. Error"] <- perm_i_res$cov.scaled[lag_roi, lag_roi] %>% as.matrix() %>% sum() %>% sqrt()
        
        df_i_est[aggr_names[1],"t value"]  <- (df_i_est[aggr_names[1],"Estimate"] - 0)/(df_i_est[aggr_names[1],"Std. Error"])
        #df_i_est[aggr_names[2],"t value"]  <- (df_i_est[aggr_names[2],"Estimate"] - 0)/(df_i_est[aggr_names[2],"Std. Error"])
        
        perm_res_list[[i]] <- df_i_est
      }  
    }
    #end_time <- Sys.time()
    #end_time - start_time
    
    # 10 replications take 5 minutes 
    # so for 1000 replications, it would take 500/60 = 8.3 hours
    
    # without the regression, the loop takes 3s. which, replicated 1000 times mores, makes less than an hour (0.83). 
    # to apply to all 20 crops in three continents this would take between 50 and 60h 
    
    # on server 10 replications is 2.16 minutes with all available threads used. 
    # 216/60 # 4 hours for 1000 replications, 40h for 10000
    
    # Let's return both the coefficients and the t-statistics (and randomization inference p values of them)
    # so this returns estimates for every permutation in columns, and for each regressor in row 
    all_coeffs <- bind_cols(lapply(perm_res_list, FUN = function(x){dplyr::select(x, Estimate)}), 
                            .name_repair = "minimal")
    colnames(all_coeffs) <- NULL
    
    # t-statistics
    all_tstats <- bind_cols(lapply(perm_res_list, FUN = function(x){dplyr::select(x, `t value`)}), 
                            .name_repair = "minimal")
    colnames(all_tstats) <- NULL
    
    # randomization inference p values are equal to the share of a parameter estimated over the randomized distribution, that is higher than the original estimate 
    
    ## Two-sided t tests
    raninf_pval_coeffs <- c()
    raninf_pval_coeffs[1] <- sum( abs(as.numeric(all_coeffs[aggr_names[1], ])) > abs(df_res[aggr_names[1], "Estimate"]) ) / (i)
    #raninf_pval_coeffs[2] <- sum( abs(as.numeric(all_coeffs[aggr_names[2], ])) > abs(df_res[aggr_names[2], "Estimate"]) ) / (i)
    names(raninf_pval_coeffs) <- aggr_names
    
    # this returns the share of absolute t ratios that are larger than the one we obtain with the actual data
    raninf_pval_tstats <- c()
    raninf_pval_tstats[1] <- sum( abs(as.numeric(all_tstats[aggr_names[1], ])) > abs(df_res[aggr_names[1], "t value"]) ) / (i)
    #raninf_pval_tstats[2] <- sum( abs(as.numeric(all_tstats[aggr_names[2], ])) > abs(df_res[aggr_names[2], "t value"]) ) / (i)
    names(raninf_pval_tstats) <- aggr_names
    
    ## One-sided t tests
    raninf_possided_pval_tstats <- c()
    raninf_possided_pval_tstats[1] <- sum( as.numeric(all_tstats[aggr_names[1], ]) > df_res[aggr_names[1], "t value"] ) / (i)
    #raninf_pval_tstats[2] <- sum( abs(as.numeric(all_tstats[aggr_names[2], ])) > abs(df_res[aggr_names[2], "t value"]) ) / (i)
    names(raninf_possided_pval_tstats) <- aggr_names
    
    raninf_negsided_pval_tstats <- c()
    raninf_negsided_pval_tstats[1] <- sum( as.numeric(all_tstats[aggr_names[1], ]) < df_res[aggr_names[1], "t value"] ) / (i)
    #raninf_pval_tstats[2] <- sum( abs(as.numeric(all_tstats[aggr_names[2], ])) > abs(df_res[aggr_names[2], "t value"]) ) / (i)
    names(raninf_negsided_pval_tstats) <- aggr_names
    
  }
  
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
  if(output == "est_object"){
    toreturn <- alpha_reg_res
  }
  
  if(output == "data"){
    toreturn <- list(df_res, d_clean)
  }
  # if(output == "est_obj"){
  #   toreturn <- sum_res
  # }
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



