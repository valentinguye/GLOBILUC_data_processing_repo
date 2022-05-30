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
                   "fixest", "boot", "urca", "plm",   #"msm", "car",  "sandwich", "lmtest",  "multcomp",
                   "ggplot2", "dotwhisker", "viridis", "hrbrthemes", #"tmap",# "leaflet", "htmltools"
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
dir.create(here("temp_data","reg_results", "rfs", "result_lists"))



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
  "Sugar", "eaear_Sugarcane", 
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

agri_crops <- eaear_mapmat[,"Crops"][!(eaear_mapmat[,"Crops"] %in% c("eaear_Oilpalm", "eaear_Rubber"))]
plantation_crops <- c("eaear_Oilpalm", "eaear_Rubber")

limited_crops1 <- c("eaear_Maizegrain", "eaear_Fodder", "eaear_Cereals", "eaear_Soy_compo", 
                    "eaear_Oilfeed_crops", "eaear_Sugarcane", "eaear_Oilpalm", "eaear_Tobacco",
                    "eaear_Cocoa_Coffee")

### MAIN DATA SETS #### 
main_data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "loss_commodity_aeay_long_final.Rdata"))
# release some memory upfront
main_data <- dplyr::filter(main_data, year >= 2008, year <= 2019)

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
all_rfs_treatments <- grep(pattern = "statute_conv", names(prices), value = TRUE)


#### EAEAR CORRELATION MATRIX ####
# work on the cross section
mdcs <- main_data[!duplicated(main_data$grid_id),]

cor_mat_abs <-  cor(dplyr::select(mdcs, all_of(eaear_mapmat[,"Crops"])))
#cor_mat_std2 <-  cor(dplyr::select(si, all_of(paste0(mapmat_si[,"Crops"], "_std2"))), use = "complete.obs")

cortests <- cor_mat_abs
j_exposures <- list()
length(j_exposures) <- length(eaear_mapmat[,"Crops"])
names(j_exposures) <- eaear_mapmat[,"Crops"]

for(serie1 in colnames(cortests)){
  
  for(serie2 in row.names(cortests)){
    
    cortests[serie1, serie2] <- cor.test(mdcs[,serie1], mdcs[,serie2], )$p.value
  }
  j_exposures[[serie1]] <- cor_mat_abs[cortests[serie1, ] < 0.05, serie1]
}
# store, for each crop, which other crop it is most correlated with
corr_mapmat <- cbind(eaear_mapmat[,"Crops"],NA)
colnames(corr_mapmat) <- c("Crops", "fst_corr")
for(crop in eaear_mapmat[,"Crops"]){
  x <- j_exposures[[crop]]
  corr_mapmat[corr_mapmat[,"Crops"]==crop, "fst_corr"] <- names(x)[x == max(x[x<max(x)])] 
}



#### DES STATS RFS #### 
# add PSD data to the chart 
psd <- readRDS(here("temp_data", "prepared_psd.Rdata"))
psd$us_ah_maize_mha <- psd$UnitedStates.Area_Harvested.Maize / 1000

#head(psd)
w_rfs <- left_join(rfs, psd[,c("year", "us_ah_maize_mha")], by = "year")

w_rfs <- w_rfs %>% 
  dplyr::filter(year>=2001 & year<=2020) %>% 
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

coeff <- 1

w_rfs <- mutate(w_rfs, usedvar = rep((w_rfs[w_rfs$`RFS mandates (right axis)`== "Conventional biofuels, final", 4] + 
                                        w_rfs[w_rfs$`RFS mandates (right axis)`== "Conventional biofuels, statutory", 4])/coeff, 
                                     4))
w_rfs[w_rfs$year<2008, "usedvar"] <- NA
# usedvar <- rfs %>% dplyr::filter(year>=2001 & year<=2020) %>% dplyr::select(year, statute_conv)
# w_rfs <- inner_join(w_rfs, usedvar, "year")

ggplot(w_rfs, aes(x=year, y=`Billions of gallons (ethanol-equivalent)`/coeff, fill=`RFS mandates (right axis)`)) + 
  geom_area(alpha=0.6 , size=.5, colour="white") +
  geom_line(aes(x = year, y = us_ah_maize_mha)) +
  geom_line(aes(x = year, y = usedvar), col = "red", size = 1) +
  #geom_line(aes(x = year, y = statute_conv, col = "red")) +
  scale_y_continuous(name = "Maize area harvested in US, Mha",
                     sec.axis = sec_axis(~.*coeff, name="Billions of gallons (ethanol-equivalent)")) +# 
  scale_fill_viridis(discrete = T) +
  theme_ipsum() + 
  scale_x_continuous(n.breaks = 5) +
  #ggtitle("RFS mandates, statutory and final, ramping up for advanced and conventinal biofuels") + 
  theme(plot.title = element_text(size = 10, face = "bold"), 
        axis.title.y.left=element_text(size=10,face="bold", hjust = 1),
        axis.title.y.right=element_text(size=10,face="bold", hjust = 1),
        axis.title.x=element_text(size=10,face="bold", hjust = 0.5), 
        panel.grid = element_line(inherit.blank = TRUE))


#### DES STATS EAEAR #### 


### MAP OF TROPICAL COMMODITY DRIVEN DEFORESTATION ### 
cdl <- brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", "loss_commo_resampledgaez_0119.tif")) 

cdl_accu <- sum(cdl[[09:19]], na.rm = TRUE)
vcdl <- values(cdl_accu)
summary(vcdl)
rm(vcdl)
cdl_accu <- reclassify(cdl_accu, cbind(0, NA))
#plot(cdl_accu)
a <- area(cdl_accu)
cdl_accu <- stack(cdl_accu, a)
accu_pct <- overlay(cdl_accu, fun = function(x, y){pct <- x/(100*y)
pct <- if_else(pct>1, 1, pct)
return(pct)})
vpct <- values(accu_pct)
summary(vpct)

america_coords <- matrix(c(-95, 25, -33, 25,
                           -33, -30, -95, -30, 
                           -95, 25), ncol = 2, byrow = TRUE)

america_ext <- st_polygon(list(america_coords)) %>% st_sfc(crs = 4326)

africa_coords <- matrix(c(-19, 20, 25, 20,
                          25, -10, -19, -10, 
                          -19, 20), ncol = 2, byrow = TRUE)

africa_ext <- st_polygon(list(africa_coords)) %>% st_sfc(crs = 4326)

asia_coords <- matrix(c(80, 30, 150, 30,
                        150, -10, 80, -10, 
                        80, 30), ncol = 2, byrow = TRUE)

asia_ext <- st_polygon(list(asia_coords)) %>% st_sfc(crs = 4326)



sfc <- c(asia_ext, america_ext, africa_ext)
continents <- st_sf(data.frame(continent_name = c("Asia", "America", "Africa"), geom = sfc))


land <- st_read(here("input_data", "ne_50m_land"))
unique(land$scalerank)
land <- land[land$scalerank==0, c("geometry")]
#plot(land)
spLand <- as(land, "Spatial")


am <- crop(accu_pct, continents[continents$continent_name=="America",])
af <- crop(accu_pct, continents[continents$continent_name=="Africa",])
as <- crop(accu_pct, continents[continents$continent_name=="Asia",])

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

pam
paf
pas

# make some comparisons with tmf
tmf  <- readRDS(here("temp_data", "merged_datasets", "tmf_aoi", "tmf_aeay_pantrop_long_final_1990_2020.Rdata"))

tmf <- dplyr::filter(tmf, year >=2011 & year <= 2019)
losscommo <- dplyr::filter(main_data, year >=2011 & year <= 2019)

tmf <- dplyr::mutate(tmf, tmf_deforestation = tmf_agri + tmf_plantation)

tmf_den <- density(tmf$tmf_deforestation)

losscommo_den <- density(losscommo$loss_commodity)

tmf_as <- tmf[tmf$continent == "Asia",]
losscommo_as <- losscommo_as[losscommo_as$continent == "Asia",]

summary(tmf_as$tmf_deforestation)
summary(losscommo_as$loss_commodity)

gaez_as <- crop(gaez, continents[continents$continent_name=="Asia",])



#### REGRESSION FUNCTION #### 

### TEMPORARY OBJECTS 
outcome_variable = "loss_commodity" 
start_year = 2010
end_year = 2019
continent = "all"

pre_process <- FALSE
pre_processed_data <- pre_d_clean_agri

rfs_rando <- ""
original_rfs_treatments <- c("statute_conv")
rfs_lag <- 2
rfs_lead <- 3
rfs_fya <-  0
rfs_pya <- 0
lag_controls = NULL
aggr_dyn <- TRUE
exposure_rfs <- eaear_mapmat[,"Crops"]#"eaear_Oilpalm"
group_exposure_rfs <- FALSE
control_all_absolute_rfs <- FALSE
most_correlated_only = FALSE
annual_rfs_controls <- FALSE

control_pasture <- FALSE
pasture_trend <- FALSE

fc_trend <- FALSE
s_trend <- FALSE
s_trend_loga <- FALSE
fc_s_trend <- FALSE

sjpos <- FALSE # should the sample be restricted to cells where sj is positive? 

fe = "grid_id + country_year" #   
preclean_level = "FE"
distribution <- "quasipoisson"
offset <- FALSE
invhypsin = TRUE
conley_cutoff <- 100
clustering = "oneway" # either "oneway" or "twoway". If oneway, it clusters on cluster_var1. 
cluster_var1 = "grid_id_10" 
cluster_var2 = "grid_id_5_year"
# dyn_tests = TRUE
output = "est_object"
glm_iter <- 25 




rm(outcome_variable, start_year, end_year, continent, 
   original_rfs_treatments, rfs_lead, rfs_lag, exposure_rfs, group_exposure_rfs, control_absolute_rfs, control_all_absolute_rfs, remaining, sjpos, fe, distribution, invhypsin, conley_cutoff, se, boot_cluster, coefs_to_aggregate, 
   output, glm_iter)

make_main_reg <- function(pre_process = FALSE, 
                          pre_processed_data = NULL,
                          outcome_variable = "loss_commodity", # one of "nd_first_loss", "first_loss", "firstloss_glassgfc", "phtf_loss"
                          start_year = 2011, 
                          end_year = 2019, 
                          continent = "all", # one of "Africa", "America", "Asia", or "all"
                          
                          rfs_rando = "", # either "between", "within", or any other string. If one of the former two, randomization inference of the type is performed
                          original_rfs_treatments = c("statute_conv"),
                          rfs_lead = 0,
                          rfs_lag = 0,
                          rfs_fya = 0, 
                          rfs_pya = 0,
                          lag_controls = NULL, # which of the lags specified above should be CONTROLLED for rather than counted in the cumulative effect
                          aggr_dyn = TRUE, # whether to report aggregate coefficients of all leads and lags ("all") or leads and lags separately ("leadlag"), or no aggregate (any other string)
                          exposure_rfs = "eaear_Soy_compo",
                          all_exposures_rfs = eaear_mapmat[,"Crops"],
                          group_exposure_rfs = FALSE, # This is deprecated, as it was relevant only when rfs exposures where standardized and could thus be added. If exposure_rfs is length 1, it does not matter whether this option is TRUE or FALSE.
                          control_all_absolute_rfs = TRUE, # there is not really an alternative where this could be FALSE and have sense. 
                          most_correlated_only = FALSE, # but this restricts the controls to only interactions with the most correlated crop. 
                          annual_rfs_controls = FALSE,
                          
                          control_pasture = FALSE,
                          pasture_trend = FALSE,
                          remaining = FALSE, # should remaining forest be controlled for STOP DOING THIS BECAUSE IT INTRODUCES NICKELL BIAS
                          s_trend = TRUE,
                          s_trend_loga = FALSE,
                          fc_trend = FALSE,
                          fc_s_trend = FALSE,
                          
                          sjpos = FALSE, # should the sample be restricted to cells where sj is positive? 
                          
                          fe = "grid_id + country_year", 
                          preclean_level = "FE", 
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
                          glm_iter = 25,
                          # dyn_tests = FALSE, # should the Fisher-type panel unit root test be returned, instead of the regressions, for the outcome_variable and the first regressor
                          output = "coef_table" # one of "data", est_object, or "coef_table" 
){
  
  # Define the outcome_variable based on the crop under study (if we are not in the placebo case)
  if(outcome_variable == "tmf_agri" & ("eaear_Oilpalm" %in% exposure_rfs | "eaear_Rubber" %in% exposure_rfs)){
    outcome_variable <- "tmf_plantation"
  }
  
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
    
    if(grepl("tmf_", outcome_variable)){
      d  <- readRDS(here("temp_data", "merged_datasets", "tmf_aoi", "tmf_aeay_pantrop_long_final_1990_2020.Rdata"))
      # release some memory upfront
      d <- dplyr::filter(d, year >= 2008, year <= 2019)
      
      d <- dplyr::mutate(d, tmf_deforestation = tmf_agri + tmf_plantation)
    }
    
    # code below should work whether d is from tmf or losscommo 
    # Keep only in data the useful variables 
    d <- dplyr::select(d, all_of(unique(c("grid_id", "year", "lat", "lon", "continent_name", "country_name", "country_year",
                                          # "fc_2000", "fc_2009", "remaining_fc", # "accu_defo_since2k",
                                          "grid_id_5", "grid_id_10", "grid_id_20", "grid_id_5_year", "grid_id_10_year", "grid_id_20_year",
                                          outcome_variable,# "tmf_agri", "tmf_flood", "tmf_plantation",
                                          "pasture_share_2000",
                                          eaear_mapmat[,"Crops"], all_exposures_rfs )))) #sj, 
    
    # Merge only the prices needed, not the whole price dataframe
    d <- left_join(d, prices[,c("year", unique(c(rfs_treatments, all_rfs_treatments)))], by = c("year"))#, all_treatments
  } 
  
  if((distribution == "gaussian") & invhypsin){
    # transform dependent variable, if gaussian GLM 
    d <- dplyr::mutate(d, !!as.symbol(outcome_variable) := asinh(!!as.symbol(outcome_variable)))
  }
  
  
  ### SPECIFICATION #### 
  
  ### Regressors 
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
  
  ### Controls
  # it's important that this is not conditioned on anything so these objects exist
  controls <- c()
  
  # transfer lags of the treatment of interest from regressors to controls
  if(length(lag_controls)>0){
    lags_to_transfer <- c()
    for(rfs_var in original_rfs_treatments){
      lags_to_transfer <- c(lags_to_transfer, paste0(exposure_rfs, "_X_", rfs_var, "_lag", lag_controls))
    }
    
    regressors <- regressors[regressors != lags_to_transfer]
    controls <- c(controls, lags_to_transfer)
    # no need to change their names, they are not captured by cumulative effects as long as they are in controls. 
    
    # old_controls <- c()
    # for(old_rc in lags_to_transfer){
    #   varname <- paste0("jexp_X_", old_rc)
    #   controls <- c(controls, varname)
    #   d <- mutate(d, !!as.symbol(varname) := !!as.symbol(exposure_rfs) * !!as.symbol(old_rc))
    # }
  }
  
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
  
  if(s_trend_loga){
    for(eaear_exp_rfs in exposure_rfs){
      varname <- paste0(eaear_exp_rfs, "_trend_expo")
      controls <- c(controls, varname)
      #make the logarithmic trend start when RFS actually starts, such that it fits best 
      # so it should be log(1) in 2005 (although never occurring in the data used for analysis but the log trend starts from there) 
      d <- mutate(d, !!as.symbol(varname) := !!as.symbol(eaear_exp_rfs) * log((year-2004))) 
    }    
  }
  
  if(s_trend){
    # Suitability time trend
    # d <- mutate(d, suitability_trend := !!as.symbol(sj)*(year))#-2000 # doesn't change anything that it's multiplied by 1...19 or 2001...2019
    # controls <- c(controls, "suitability_trend")
    for(eaear_exp_rfs in exposure_rfs){
      varname <- paste0(eaear_exp_rfs, "_trend")
      controls <- c(controls, varname)
      d <- mutate(d, !!as.symbol(varname) := !!as.symbol(eaear_exp_rfs) * (year-2005)) # for the linear trend it does not matter
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
    
    used_vars <- unique(c("grid_id", "year", "lat", "lon","continent_name",  "country_year",  #"country_name",  "remaining_fc", "accu_defo_since2k", # "sj_year",
                          "grid_id_5", "grid_id_10", "grid_id_20", "grid_id_5_year", "grid_id_10_year", "grid_id_20_year",
                          outcome_variable,# "tmf_agri", "tmf_flood", "tmf_plantation",
                          regressors, controls, 
                          exposure_rfs, original_rfs_treatments, rfs_treatments))    # this is necessary to reconstruct variables in randomization inference processes
    
    
    
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
  # In this case, we are interested in LEADS only
  if(aggr_dyn == "leads"){ 
    # deprecated: in this case, there ARE lags in the specification, but we are not interested in them
    
    df_res <- rbind(rep(NA, ncol(df_res)), df_res)
    
    # ORDER MATTERS
    aggr_names <- paste0(exposure_rfs, c("_X_aggrleads"))
    row.names(df_res)[1] <- aggr_names
    # regressors of interest
    # base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
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
  }
  # In this case, we are interested in LAGS only
  if(aggr_dyn == "lags"){
    df_res <- rbind(rep(NA, ncol(df_res)), df_res)
    
    # ORDER MATTERS
    aggr_names <- paste0(exposure_rfs, c("_X_aggrlags"))
    row.names(df_res)[1] <- aggr_names
    # regressors of interest
    base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
    # Contemporaneous value is in lags
    lag_roi <- c(base_reg_name, 
                 grep(pattern = paste0(base_reg_name,"_lag"), 
                      regressors, value = TRUE))
    
    df_res[aggr_names[1],"Estimate"] <- reg_res$coefficients[lag_roi] %>% sum()
    
    # select the part of the VCOV matrix that is to be used to compute the standard error of the sum
    # use formula for variance of sum of random variables : https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables
    df_res[aggr_names[1],"Std. Error"] <- reg_res$cov.scaled[lag_roi, lag_roi] %>% as.matrix() %>% sum() %>% sqrt()
    
    df_res[aggr_names[1],"t value"]  <- (df_res[aggr_names[1],"Estimate"] - 0)/(df_res[aggr_names[1],"Std. Error"])
    
    # use t distribution with degrees of freedom equal to that used by fixest, i.e. after two way cluster adjustment.
    # does not make a significant difference given sample size
    df_res[aggr_names[1],"Pr(>|t|)"]  <- (2*pt(abs(df_res[aggr_names[1],"t value"]), 
                                               lower.tail = FALSE, 
                                               df = fixest_df)) 
  }
  
  
  if(aggr_dyn & control_all_absolute_rfs & rfs_lead > 0 & rfs_lag > 0 & rfs_fya == 0 & rfs_pya == 0){
    # In this case, we are interested in LEAD AND LAG effects, aggregated separately, and all together.
    df_res <- rbind(rep(NA, ncol(df_res)), rep(NA, ncol(df_res)), df_res)
    
    # ORDER MATTERS
    aggr_names <- paste0(exposure_rfs, c("_X_aggrleads", "_X_aggrlags"))
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
    aggr_names <- paste0(exposure_rfs, c("_X_aggrall"))
    row.names(df_res)[1] <- aggr_names
    # regressors of interest
    base_reg_name <- paste0(exposure_rfs,"_X_",original_rfs_treatments)
    
    # Contemporaneous, lead, and lag values
    all_roi <- c(base_reg_name, 
                 grep(pattern = paste0(base_reg_name,"_lead"), 
                      regressors, value = TRUE),
                 grep(pattern = paste0(base_reg_name,"_lag"), 
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
  if(aggr_dyn & control_all_absolute_rfs & rfs_lead == 0 & rfs_lag == 0 & rfs_fya > 0 & rfs_pya > 0){
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
  if(aggr_dyn & control_all_absolute_rfs & rfs_lead > 0 & rfs_lag == 0 & rfs_fya == 0 & rfs_pya > 0){
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
  
  if(aggr_dyn & rfs_lead > 0 & rfs_lag > 0 & rfs_fya == 0 & rfs_pya == 0 & 
     length(controls) == 0 # this is the case when there is no trend AND control_all_absolute_rfs is FALSE 
  ){
    for(EOI in exposure_rfs){
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
      
      # Contemporaneous, lead, and lag values
      all_roi <- c(base_reg_name, 
                   grep(pattern = paste0(base_reg_name,"_lead"), 
                        regressors, value = TRUE),
                   grep(pattern = paste0(base_reg_name,"_lag"), 
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
  }
  # take data set as exactly used in estimation - NOT NECESSARY anymore, given the precleaning
  # if(length(reg_res$obs_selection) > 0){
  #   d_clean <- d_clean[reg_res$obs_selection[[1]], ]
  # }
  
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
  
  
  rm(d_clean, df_res)
  return(toreturn)
  rm(toreturn)
}

### MAIN SPECIFICATION ### 
est_parameters <- list(outcome_variable  = "loss_commodity",
                       continent = NULL, # this is filled within the loop 
                       start_year = 2009, 
                       end_year = 2019, 
                       most_correlated_only = FALSE,
                       annual_rfs_controls = FALSE,
                       all_exposures_rfs = eaear_mapmat[,"Crops"],
                       sjpos = FALSE,
                       lags = 1,
                       leads = 3,
                       fya = 0, 
                       pya = 0,
                       lag_controls = NULL,
                       control_pasture = FALSE,
                       s_trend = FALSE, 
                       s_trend_loga = FALSE,
                       fe = "grid_id + country_year",
                       clustering = "oneway",
                       cluster_var1 = "grid_id_10",
                       cluster_var2 = "grid_id_10",
                       distribution = "quasipoisson")


#### RFS CUMMULATIVE LEADS & LAGS 2011-2019 ####
rfs_1lag3lead_notrend_clt10 <- list(all = list(), 
                                    America = list(), 
                                    Africa = list(), 
                                    Asia = list())

# prepare data set in advance, to repeat in all regressions
CNT <- "America"

for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters[["continent"]] <- CNT
  
  d <- main_data
  
  # Keep only in data the useful variables 
  d <- dplyr::select(d, all_of(unique(c("grid_id", "year", "lat", "lon", "continent_name", "country_name", "country_year",
                                        # "fc_2000", "fc_2009", "remaining_fc", # "accu_defo_since2k",
                                        "grid_id_5", "grid_id_10", "grid_id_20", "grid_id_5_year", "grid_id_10_year", "grid_id_20_year",
                                        # "tmf_agri", "tmf_flood", "tmf_plantation", "tmf_deforestation",
                                        "loss_commodity",
                                        "pasture_share_2000",
                                        eaear_mapmat[,"Crops"] ))))
  
  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", unique(c(all_rfs_treatments)))], by = c("year"))
  
  # - are suitable to crop j 
  if(est_parameters[["sjpos"]]){
    d <- dplyr::filter(d, !!as.symbol(exposure_rfs) > 0)
  }
  
  # - are in study period
  # if(start_year != 2011 | end_year != 2019){
  d <- dplyr::filter(d, year >= est_parameters[["start_year"]])
  d <- dplyr::filter(d, year <= est_parameters[["end_year"]])
  # }
  
  # - are in study area
  if(est_parameters[["continent"]] != "all"){
    d <- dplyr::filter(d, continent_name == est_parameters[["continent"]])
  }
  
  # Have tmf deforestation to agriculture at least once
  temp_est <- feglm(fml = as.formula(paste0("loss_commodity", " ~ 1 | ", est_parameters[["fe"]])),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    pre_d_clean_agri <- d[unlist(temp_est$obs_selection),]
  }  else { 
    pre_d_clean_agri <- d
  }
  
  # # Have some tmf deforestation to plantation at least once
  # temp_est <- feglm(fml = as.formula(paste0("loss_commodity", " ~ 1 | ", est_parameters[["fe"]])),
  #                   data = d,
  #                   family = "poisson")
  # # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  # if(length(temp_est$obs_selection)>0){
  #   pre_d_clean_plantation <- d[unlist(temp_est$obs_selection),]
  # }  else { 
  #   pre_d_clean_plantation <- d
  # }
  
  rm(d)
  
  # Regressions of tmf_agri
  elm <- 1
  for(rfs_exp in est_parameters[["all_exposures_rfs"]]){#agri_crops eaear_mapmat[,"Crops"]
    rfs_1lag3lead_notrend_clt10[[CNT]][[elm]] <- make_main_reg(pre_process = TRUE,
                                                               pre_processed_data = pre_d_clean_agri, # NOTICE THIS
                                                               
                                                               outcome_variable = "loss_commodity", # AND THIS
                                                               continent = est_parameters[["continent"]],
                                                               start_year = est_parameters[["start_year"]],
                                                               end_year = est_parameters[["end_year"]],
                                                               
                                                               aggr_dyn = TRUE,
                                                               exposure_rfs = rfs_exp,
                                                               all_exposures_rfs = est_parameters[["all_exposures_rfs"]],
                                                               
                                                               original_rfs_treatments = c("statute_conv"),
                                                               rfs_lead = est_parameters[["leads"]], 
                                                               rfs_lag = est_parameters[["lags"]],
                                                               rfs_fya = est_parameters[["fya"]],  #ya,
                                                               rfs_pya = est_parameters[["pya"]], #ya,
                                                               
                                                               most_correlated_only = est_parameters[["most_correlated_only"]],
                                                               annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                               
                                                               s_trend = est_parameters[["s_trend"]],
                                                               s_trend_loga = est_parameters[["s_trend_loga"]],
                                                               fe = est_parameters[["fe"]], 
                                                               cluster_var1 = est_parameters[["cluster_var1"]], 
                                                               
                                                               rfs_rando = ""
    )
    
    names(rfs_1lag3lead_notrend_clt10[[CNT]])[elm] <- rfs_exp
    elm <- elm + 1
  }
  
}

### PLOT COEFFICIENTS ### 
# prepare regression outputs in a tidy data frame readable by dwplot

rm(df)
dyn_df_list <- list()

for(CNT in c("all", "America", "Africa", "Asia")){ # 1 and 2 correspond respectively to fya and pya, or aggrleads and aggrlags
  dyn_df <- lapply(rfs_1lag3lead_notrend_clt10[[CNT]], FUN = function(x){as.data.frame(x)[1,] }) %>% bind_rows()
  dyn_df$term <- gsub(pattern = "_X_.*$", x = row.names(dyn_df), replacement = "") # replace everything after and including _X_ with nothing
  dyn_df$model <- CNT # dyn_des[dyn]
  dyn_df_list[[CNT]] <- dyn_df
}
df <- bind_rows(dyn_df_list)
df[df$model == "all", "model"] <- "Global"
#df$term <- rep(eaear_mapmat[,"Crops"], length(dyn_df_list))
df$significant01 <- ""
df[df[,"Pr(>|t|)"] < 0.1, "significant01"] <- "p-value < .1"

names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"
row.names(df) <- NULL

#title <- "Impacts of RFS ramping up on tropical deforestation, via exposures to different land uses, 2001-2019"
title <- paste0("Cummulative effects of realized and expected biofuel mandates on tropical deforestation for major commercial crops, 2011-2019")

# IF ONLY SIGNIFICANT ONES ARE FEATURED
crop_groups <- list(c("Group 1", "Cereals", "Rice"),
                    c("Group 2", "Cotton", "Oil crops"),
                    c("Group 3", "Coconut", "Oil palm"),
                    c("Group 4", "Biomass crops", "Sugar crops"), # "Groundnut","Maize","Rapeseed","Rice","Sorghum","Soy","Sugar (cane or beet)","Sunflower",
                    c("Group 5", "Tobacco", "Citrus"), # "Citrus", "Cotton",
                    c("Group 6", "Banana", "Tea") # "Cocoa", "Coffee", "Rubber", "Cocoa or Coffee
)
{dwplot(df,
        dot_args = list(size = 2.5,aes(shape = model, alpha = significant01)),#, alpha = significant01
        whisker_args = list(size = 1, aes(alpha = significant01)),#alpha = significant01
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(#eaear_Barley = "Barley",
      #eaear_Groundnut = "Groundnut", 
      eaear_Maizegrain = "Maize",
      
      eaear_Fodder = "Pasture",
      
      eaear_Cereals = "Cereals",
      eaear_Rice = "Rice",
      
      eaear_Soy_compo = "Soy", 
      
      eaear_Cotton = "Cotton",
      eaear_Oilfeed_crops = "Oil crops",
      
      #eaear_Rapeseed = "Rapeseed", 
      # eaear_Sorghum2 = "Sorghum", 
      eaear_Biomass = "Biomass crops",
      eaear_Sugarcane = "Sugar crops", 
      #eaear_Sunflower = "Sunflower",
      #eaear_Wheat = "Wheat",
      
      eaear_Coconut = "Coconut", 
      eaear_Oilpalm = "Oil palm",
      
      eaear_Tobacco = "Tobacco", 
      eaear_Citrus = "Citrus",
      
      eaear_Banana = "Banana", 
      eaear_Cocoa_Coffee = "Cocoa or Coffee",
      #eaear_Cocoa  = "Cocoa", 
      #eaear_Coffee = "Coffee", 
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
          legend.title = element_blank()) } %>% 
  # don't run this if we used only limited crops
  add_brackets(crop_groups, face = "italic")


#### ROBUSTNESS CHECKS #### 
# parameters (edited in single robustness checks)
est_parameters <- list(outcome_variable  = "loss_commodity",
                       continent = NULL, # this is filled within the loop 
                       start_year = 2009, 
                       end_year = 2019, 
                       most_correlated_only = FALSE,
                       annual_rfs_controls = FALSE,
                       all_exposures_rfs = eaear_mapmat[,"Crops"],
                       sjpos = FALSE,
                       lags = 1,
                       leads = 3,
                       fya = 0, 
                       pya = 0,
                       lag_controls = NULL,
                       control_pasture = FALSE,
                       s_trend = FALSE, 
                       s_trend_loga = FALSE,
                       fe = "grid_id + country_year",
                       clustering = "oneway",
                       cluster_var1 = "grid_id_10",
                       cluster_var2 = "grid_id_10",
                       distribution = "quasipoisson")

# prepare data common to all robustness checks 
d_all <- main_data

# Keep only in data the useful variables 
d_all <- dplyr::select(d_all, all_of(unique(c("grid_id", "year", "lat", "lon", "continent_name", "country_name", "country_year",
                                              # "fc_2000", "fc_2009", "remaining_fc", # "accu_defo_since2k",
                                              "grid_id_5", "grid_id_10", "grid_id_20", "grid_id_5_year", "grid_id_10_year", "grid_id_20_year",
                                              # "tmf_agri", "tmf_flood", "tmf_plantation", "tmf_deforestation",
                                              "loss_commodity",
                                              "pasture_share_2000",
                                              eaear_mapmat[,"Crops"] ))))

# Merge only the prices needed, not the whole price dataframe
d_all <- left_join(d_all, prices[,c("year", unique(c(all_rfs_treatments)))], by = c("year"))

# - are suitable to crop j 
if(est_parameters[["sjpos"]]){
  d_all <- dplyr::filter(d_all, !!as.symbol(exposure_rfs) > 0)
}

# - are in study period
# if(start_year != 2011 | end_year != 2019){
d_all <- dplyr::filter(d_all, year >= est_parameters[["start_year"]])
d_all <- dplyr::filter(d_all, year <= est_parameters[["end_year"]])
# }

#### Dynamics #### 
## 1 LEAD
est_parameters[["leads"]] <- 1
rfs_1lag1lead_notrend_clt10 <- list(all = list(), 
                                    America = list(), 
                                    Africa = list(), 
                                    Asia = list())


for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters[["continent"]] <- CNT
  
  # - are in study area
  if(est_parameters[["continent"]] != "all"){
    d <- dplyr::filter(d_all, continent_name == est_parameters[["continent"]])
  } else {
    d <- d_all
  }
  
  # Have tmf deforestation to agriculture at least once
  temp_est <- feglm(fml = as.formula(paste0("loss_commodity", " ~ 1 | ", est_parameters[["fe"]])),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    pre_d_clean_agri <- d[unlist(temp_est$obs_selection),]
  }  else { 
    pre_d_clean_agri <- d
  }
  
  rm(d)
  
  # Regressions of tmf_agri
  elm <- 1
  for(rfs_exp in est_parameters[["all_exposures_rfs"]]){#agri_crops eaear_mapmat[,"Crops"]
    rfs_1lag1lead_notrend_clt10[[CNT]][[elm]] <- make_main_reg(pre_process = TRUE,
                                                               pre_processed_data = pre_d_clean_agri, # NOTICE THIS
                                                               
                                                               outcome_variable = "loss_commodity", # AND THIS
                                                               continent = est_parameters[["continent"]],
                                                               start_year = est_parameters[["start_year"]],
                                                               end_year = est_parameters[["end_year"]],
                                                               
                                                               aggr_dyn = TRUE,
                                                               exposure_rfs = rfs_exp,
                                                               all_exposures_rfs = est_parameters[["all_exposures_rfs"]],
                                                               
                                                               original_rfs_treatments = c("statute_conv"),
                                                               rfs_lead = est_parameters[["leads"]], 
                                                               rfs_lag = est_parameters[["lags"]],
                                                               rfs_fya = est_parameters[["fya"]],  #ya,
                                                               rfs_pya = est_parameters[["pya"]], #ya,
                                                               
                                                               most_correlated_only = est_parameters[["most_correlated_only"]],
                                                               annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                               
                                                               s_trend = est_parameters[["s_trend"]],
                                                               s_trend_loga = est_parameters[["s_trend_loga"]],
                                                               fe = est_parameters[["fe"]], 
                                                               cluster_var1 = est_parameters[["cluster_var1"]], 
                                                               
                                                               rfs_rando = ""
    )
    
    names(rfs_1lag1lead_notrend_clt10[[CNT]])[elm] <- rfs_exp
    elm <- elm + 1
  }
  
}

## 2 LEADS
est_parameters[["leads"]] <- 2
rfs_1lag2lead_notrend_clt10 <- list(all = list(), 
                                    America = list(), 
                                    Africa = list(), 
                                    Asia = list())


for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters[["continent"]] <- CNT
  
  # - are in study area
  if(est_parameters[["continent"]] != "all"){
    d <- dplyr::filter(d_all, continent_name == est_parameters[["continent"]])
  } else {
    d <- d_all
  }
  
  # Have tmf deforestation to agriculture at least once
  temp_est <- feglm(fml = as.formula(paste0("loss_commodity", " ~ 1 | ", est_parameters[["fe"]])),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    pre_d_clean_agri <- d[unlist(temp_est$obs_selection),]
  }  else { 
    pre_d_clean_agri <- d
  }
  
  rm(d)
  
  # Regressions of tmf_agri
  elm <- 1
  for(rfs_exp in est_parameters[["all_exposures_rfs"]]){#agri_crops eaear_mapmat[,"Crops"]
    rfs_1lag2lead_notrend_clt10[[CNT]][[elm]] <- make_main_reg(pre_process = TRUE,
                                                               pre_processed_data = pre_d_clean_agri, # NOTICE THIS
                                                               
                                                               outcome_variable = "loss_commodity", # AND THIS
                                                               continent = est_parameters[["continent"]],
                                                               start_year = est_parameters[["start_year"]],
                                                               end_year = est_parameters[["end_year"]],
                                                               
                                                               aggr_dyn = TRUE,
                                                               exposure_rfs = rfs_exp,
                                                               all_exposures_rfs = est_parameters[["all_exposures_rfs"]],
                                                               
                                                               original_rfs_treatments = c("statute_conv"),
                                                               rfs_lead = est_parameters[["leads"]], 
                                                               rfs_lag = est_parameters[["lags"]],
                                                               rfs_fya = est_parameters[["fya"]],  #ya,
                                                               rfs_pya = est_parameters[["pya"]], #ya,
                                                               
                                                               most_correlated_only = est_parameters[["most_correlated_only"]],
                                                               annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                               
                                                               s_trend = est_parameters[["s_trend"]],
                                                               s_trend_loga = est_parameters[["s_trend_loga"]],
                                                               fe = est_parameters[["fe"]], 
                                                               cluster_var1 = est_parameters[["cluster_var1"]], 
                                                               
                                                               rfs_rando = ""
    )
    
    names(rfs_1lag2lead_notrend_clt10[[CNT]])[elm] <- rfs_exp
    elm <- elm + 1
  }
  
}




## 4 LEADS
est_parameters[["leads"]] <- 4
rfs_1lag4lead_notrend_clt10 <- list(all = list(), 
                                    America = list(), 
                                    Africa = list(), 
                                    Asia = list())


for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters[["continent"]] <- CNT
  
  # - are in study area
  if(est_parameters[["continent"]] != "all"){
    d <- dplyr::filter(d_all, continent_name == est_parameters[["continent"]])
  } else {
    d <- d_all
  }
  
  # Have tmf deforestation to agriculture at least once
  temp_est <- feglm(fml = as.formula(paste0("loss_commodity", " ~ 1 | ", est_parameters[["fe"]])),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    pre_d_clean_agri <- d[unlist(temp_est$obs_selection),]
  }  else { 
    pre_d_clean_agri <- d
  }
  
  rm(d)
  
  # Regressions of tmf_agri
  elm <- 1
  for(rfs_exp in est_parameters[["all_exposures_rfs"]]){#agri_crops eaear_mapmat[,"Crops"]
    rfs_1lag4lead_notrend_clt10[[CNT]][[elm]] <- make_main_reg(pre_process = TRUE,
                                                               pre_processed_data = pre_d_clean_agri, # NOTICE THIS
                                                               
                                                               outcome_variable = "loss_commodity", # AND THIS
                                                               continent = est_parameters[["continent"]],
                                                               start_year = est_parameters[["start_year"]],
                                                               end_year = est_parameters[["end_year"]],
                                                               
                                                               aggr_dyn = TRUE,
                                                               exposure_rfs = rfs_exp,
                                                               all_exposures_rfs = est_parameters[["all_exposures_rfs"]],
                                                               
                                                               original_rfs_treatments = c("statute_conv"),
                                                               rfs_lead = est_parameters[["leads"]], 
                                                               rfs_lag = est_parameters[["lags"]],
                                                               rfs_fya = est_parameters[["fya"]],  #ya,
                                                               rfs_pya = est_parameters[["pya"]], #ya,
                                                               
                                                               most_correlated_only = est_parameters[["most_correlated_only"]],
                                                               annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                               
                                                               s_trend = est_parameters[["s_trend"]],
                                                               s_trend_loga = est_parameters[["s_trend_loga"]],
                                                               fe = est_parameters[["fe"]], 
                                                               cluster_var1 = est_parameters[["cluster_var1"]], 
                                                               
                                                               rfs_rando = ""
    )
    
    names(rfs_1lag4lead_notrend_clt10[[CNT]])[elm] <- rfs_exp
    elm <- elm + 1
  }
  
}





## 5 LEADS
est_parameters[["leads"]] <- 5
rfs_1lag5lead_notrend_clt10 <- list(all = list(), 
                                    America = list(), 
                                    Africa = list(), 
                                    Asia = list())


for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters[["continent"]] <- CNT
  
  # - are in study area
  if(est_parameters[["continent"]] != "all"){
    d <- dplyr::filter(d_all, continent_name == est_parameters[["continent"]])
  } else {
    d <- d_all
  }
  
  # Have tmf deforestation to agriculture at least once
  temp_est <- feglm(fml = as.formula(paste0("loss_commodity", " ~ 1 | ", est_parameters[["fe"]])),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    pre_d_clean_agri <- d[unlist(temp_est$obs_selection),]
  }  else { 
    pre_d_clean_agri <- d
  }
  
  rm(d)
  
  # Regressions of tmf_agri
  elm <- 1
  for(rfs_exp in est_parameters[["all_exposures_rfs"]]){#agri_crops eaear_mapmat[,"Crops"]
    rfs_1lag5lead_notrend_clt10[[CNT]][[elm]] <- make_main_reg(pre_process = TRUE,
                                                               pre_processed_data = pre_d_clean_agri, # NOTICE THIS
                                                               
                                                               outcome_variable = "loss_commodity", # AND THIS
                                                               continent = est_parameters[["continent"]],
                                                               start_year = est_parameters[["start_year"]],
                                                               end_year = est_parameters[["end_year"]],
                                                               
                                                               aggr_dyn = TRUE,
                                                               exposure_rfs = rfs_exp,
                                                               all_exposures_rfs = est_parameters[["all_exposures_rfs"]],
                                                               
                                                               original_rfs_treatments = c("statute_conv"),
                                                               rfs_lead = est_parameters[["leads"]], 
                                                               rfs_lag = est_parameters[["lags"]],
                                                               rfs_fya = est_parameters[["fya"]],  #ya,
                                                               rfs_pya = est_parameters[["pya"]], #ya,
                                                               
                                                               most_correlated_only = est_parameters[["most_correlated_only"]],
                                                               annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                               
                                                               s_trend = est_parameters[["s_trend"]],
                                                               s_trend_loga = est_parameters[["s_trend_loga"]],
                                                               fe = est_parameters[["fe"]], 
                                                               cluster_var1 = est_parameters[["cluster_var1"]], 
                                                               
                                                               rfs_rando = ""
    )
    
    names(rfs_1lag5lead_notrend_clt10[[CNT]])[elm] <- rfs_exp
    elm <- elm + 1
  }
  
}

## because it's the end of the series of robustness checks on leads, set the parameter to default
est_parameters[["leads"]] <- 3

## CONTROL FOR 1 MORE LAG ## 
est_parameters[["lags"]] <- 2
est_parameters[["lag_controls"]] <- 2
rfs_1lag3lead_lag2ctrl_notrend_clt10 <- list(all = list(), 
                                             America = list(), 
                                             Africa = list(), 
                                             Asia = list())


for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters[["continent"]] <- CNT
  
  # - are in study area
  if(est_parameters[["continent"]] != "all"){
    d <- dplyr::filter(d_all, continent_name == est_parameters[["continent"]])
  } else {
    d <- d_all
  }
  
  # keep only exogenous RFS years (the longer lags may take earlier values)
  d <- dplyr::filter(d, year >= 2008+est_parameters[["lags"]])
  
  # Have tmf deforestation to agriculture at least once
  temp_est <- feglm(fml = as.formula(paste0("loss_commodity", " ~ 1 | ", est_parameters[["fe"]])),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    pre_d_clean_agri <- d[unlist(temp_est$obs_selection),]
  }  else { 
    pre_d_clean_agri <- d
  }
  
  rm(d)
  
  # Regressions of tmf_agri
  elm <- 1
  for(rfs_exp in est_parameters[["all_exposures_rfs"]]){#agri_crops eaear_mapmat[,"Crops"]
    rfs_1lag3lead_lag2ctrl_notrend_clt10[[CNT]][[elm]] <- make_main_reg(pre_process = TRUE,
                                                                        pre_processed_data = pre_d_clean_agri, # NOTICE THIS
                                                                        
                                                                        outcome_variable = "loss_commodity", # AND THIS
                                                                        continent = est_parameters[["continent"]],
                                                                        start_year = est_parameters[["start_year"]],
                                                                        end_year = est_parameters[["end_year"]],
                                                                        
                                                                        aggr_dyn = TRUE,
                                                                        exposure_rfs = rfs_exp,
                                                                        all_exposures_rfs = est_parameters[["all_exposures_rfs"]],
                                                                        
                                                                        original_rfs_treatments = c("statute_conv"),
                                                                        rfs_lead = est_parameters[["leads"]], 
                                                                        rfs_lag = est_parameters[["lags"]],
                                                                        rfs_fya = est_parameters[["fya"]],  #ya,
                                                                        rfs_pya = est_parameters[["pya"]], #ya,
                                                                        lag_controls = est_parameters[["lag_controls"]],
                                                                        
                                                                        most_correlated_only = est_parameters[["most_correlated_only"]],
                                                                        annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                                        
                                                                        control_pasture = est_parameters[["control_pasture"]],
                                                                        
                                                                        s_trend = est_parameters[["s_trend"]],
                                                                        s_trend_loga = est_parameters[["s_trend_loga"]],
                                                                        fe = est_parameters[["fe"]], 
                                                                        cluster_var1 = est_parameters[["cluster_var1"]], 
                                                                        
                                                                        rfs_rando = ""
    )
    
    names(rfs_1lag3lead_lag2ctrl_notrend_clt10[[CNT]])[elm] <- rfs_exp
    elm <- elm + 1
  }
  
}

est_parameters[["lag_controls"]] <- NULL
est_parameters[["lags"]] <- 1

## ADD t-2 IN CUMULATIVE EFFECTS OF ALL CROPS ## 
est_parameters[["lags"]] <- 2
rfs_2lag3lead_notrend_clt10 <- list(all = list(), 
                                    America = list(), 
                                    Africa = list(), 
                                    Asia = list())


for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters[["continent"]] <- CNT
  
  # - are in study area
  if(est_parameters[["continent"]] != "all"){
    d <- dplyr::filter(d_all, continent_name == est_parameters[["continent"]])
  } else {
    d <- d_all
  }
  
  # keep only exogenous RFS years (the longer lags may take earlier values)
  d <- dplyr::filter(d, year >= 2008+est_parameters[["lags"]])
  
  # Have tmf deforestation to agriculture at least once
  temp_est <- feglm(fml = as.formula(paste0("loss_commodity", " ~ 1 | ", est_parameters[["fe"]])),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    pre_d_clean_agri <- d[unlist(temp_est$obs_selection),]
  }  else { 
    pre_d_clean_agri <- d
  }
  
  rm(d)
  
  # Regressions of tmf_agri
  elm <- 1
  for(rfs_exp in est_parameters[["all_exposures_rfs"]]){#agri_crops eaear_mapmat[,"Crops"]
    rfs_2lag3lead_notrend_clt10[[CNT]][[elm]] <- make_main_reg(pre_process = TRUE,
                                                               pre_processed_data = pre_d_clean_agri, # NOTICE THIS
                                                               
                                                               outcome_variable = "loss_commodity", # AND THIS
                                                               continent = est_parameters[["continent"]],
                                                               start_year = est_parameters[["start_year"]],
                                                               end_year = est_parameters[["end_year"]],
                                                               
                                                               aggr_dyn = TRUE,
                                                               exposure_rfs = rfs_exp,
                                                               all_exposures_rfs = est_parameters[["all_exposures_rfs"]],
                                                               
                                                               original_rfs_treatments = c("statute_conv"),
                                                               rfs_lead = est_parameters[["leads"]], 
                                                               rfs_lag = est_parameters[["lags"]],
                                                               rfs_fya = est_parameters[["fya"]],  #ya,
                                                               rfs_pya = est_parameters[["pya"]], #ya,
                                                               lag_controls = est_parameters[["lag_controls"]],
                                                               
                                                               most_correlated_only = est_parameters[["most_correlated_only"]],
                                                               annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                               
                                                               control_pasture = est_parameters[["control_pasture"]],
                                                               
                                                               s_trend = est_parameters[["s_trend"]],
                                                               s_trend_loga = est_parameters[["s_trend_loga"]],
                                                               fe = est_parameters[["fe"]], 
                                                               cluster_var1 = est_parameters[["cluster_var1"]], 
                                                               
                                                               rfs_rando = ""
    )
    
    names(rfs_2lag3lead_notrend_clt10[[CNT]])[elm] <- rfs_exp
    elm <- elm + 1
  }
  
}

est_parameters[["lags"]] <- 1




#### Trends ####

## LINEAR TREND  
est_parameters[["s_trend"]] <- TRUE
rfs_1lag3lead_lintrend_clt10 <- list(all = list(), 
                                     America = list(), 
                                     Africa = list(), 
                                     Asia = list())


for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters[["continent"]] <- CNT
  
  # - are in study area
  if(est_parameters[["continent"]] != "all"){
    d <- dplyr::filter(d_all, continent_name == est_parameters[["continent"]])
  } else {
    d <- d_all
  }
  
  # Have tmf deforestation to agriculture at least once
  temp_est <- feglm(fml = as.formula(paste0("loss_commodity", " ~ 1 | ", est_parameters[["fe"]])),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    pre_d_clean_agri <- d[unlist(temp_est$obs_selection),]
  }  else { 
    pre_d_clean_agri <- d
  }
  
  rm(d)
  
  # Regressions of tmf_agri
  elm <- 1
  for(rfs_exp in est_parameters[["all_exposures_rfs"]]){#agri_crops eaear_mapmat[,"Crops"]
    rfs_1lag3lead_lintrend_clt10[[CNT]][[elm]] <- make_main_reg(pre_process = TRUE,
                                                                pre_processed_data = pre_d_clean_agri, # NOTICE THIS
                                                                
                                                                outcome_variable = "loss_commodity", # AND THIS
                                                                continent = est_parameters[["continent"]],
                                                                start_year = est_parameters[["start_year"]],
                                                                end_year = est_parameters[["end_year"]],
                                                                
                                                                aggr_dyn = TRUE,
                                                                exposure_rfs = rfs_exp,
                                                                all_exposures_rfs = est_parameters[["all_exposures_rfs"]],
                                                                
                                                                original_rfs_treatments = c("statute_conv"),
                                                                rfs_lead = est_parameters[["leads"]], 
                                                                rfs_lag = est_parameters[["lags"]],
                                                                rfs_fya = est_parameters[["fya"]],  #ya,
                                                                rfs_pya = est_parameters[["pya"]], #ya,
                                                                
                                                                most_correlated_only = est_parameters[["most_correlated_only"]],
                                                                annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                                
                                                                s_trend = est_parameters[["s_trend"]],
                                                                s_trend_loga = est_parameters[["s_trend_loga"]],
                                                                fe = est_parameters[["fe"]], 
                                                                cluster_var1 = est_parameters[["cluster_var1"]], 
                                                                
                                                                rfs_rando = ""
    )
    
    names(rfs_1lag3lead_lintrend_clt10[[CNT]])[elm] <- rfs_exp
    elm <- elm + 1
  }
  
}


## LINEAR TREND  
est_parameters[["s_trend"]] <- FALSE
est_parameters[["s_trend_loga"]] <- TRUE

rfs_1lag3lead_logtrend_clt10 <- list(all = list(), 
                                     America = list(), 
                                     Africa = list(), 
                                     Asia = list())


for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters[["continent"]] <- CNT
  
  # - are in study area
  if(est_parameters[["continent"]] != "all"){
    d <- dplyr::filter(d_all, continent_name == est_parameters[["continent"]])
  } else {
    d <- d_all
  }
  
  # Have tmf deforestation to agriculture at least once
  temp_est <- feglm(fml = as.formula(paste0("loss_commodity", " ~ 1 | ", est_parameters[["fe"]])),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    pre_d_clean_agri <- d[unlist(temp_est$obs_selection),]
  }  else { 
    pre_d_clean_agri <- d
  }
  
  rm(d)
  
  # Regressions of tmf_agri
  elm <- 1
  for(rfs_exp in est_parameters[["all_exposures_rfs"]]){#agri_crops eaear_mapmat[,"Crops"]
    rfs_1lag3lead_logtrend_clt10[[CNT]][[elm]] <- make_main_reg(pre_process = TRUE,
                                                                pre_processed_data = pre_d_clean_agri, # NOTICE THIS
                                                                
                                                                outcome_variable = "loss_commodity", # AND THIS
                                                                continent = est_parameters[["continent"]],
                                                                start_year = est_parameters[["start_year"]],
                                                                end_year = est_parameters[["end_year"]],
                                                                
                                                                aggr_dyn = TRUE,
                                                                exposure_rfs = rfs_exp,
                                                                all_exposures_rfs = est_parameters[["all_exposures_rfs"]],
                                                                
                                                                original_rfs_treatments = c("statute_conv"),
                                                                rfs_lead = est_parameters[["leads"]], 
                                                                rfs_lag = est_parameters[["lags"]],
                                                                rfs_fya = est_parameters[["fya"]],  #ya,
                                                                rfs_pya = est_parameters[["pya"]], #ya,
                                                                
                                                                most_correlated_only = est_parameters[["most_correlated_only"]],
                                                                annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                                
                                                                s_trend = est_parameters[["s_trend"]],
                                                                s_trend_loga = est_parameters[["s_trend_loga"]],
                                                                fe = est_parameters[["fe"]], 
                                                                cluster_var1 = est_parameters[["cluster_var1"]], 
                                                                
                                                                rfs_rando = ""
    )
    
    names(rfs_1lag3lead_logtrend_clt10[[CNT]])[elm] <- rfs_exp
    elm <- elm + 1
  }
  
}


#### Clustering #### 
## CLUSTER 5 
est_parameters[["cluster_var1"]] <- "grid_id_5"
rfs_1lag3lead_notrend_clt5 <- list(all = list(), 
                                   America = list(), 
                                   Africa = list(), 
                                   Asia = list())


for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters[["continent"]] <- CNT
  
  # - are in study area
  if(est_parameters[["continent"]] != "all"){
    d <- dplyr::filter(d_all, continent_name == est_parameters[["continent"]])
  } else {
    d <- d_all
  }
  
  # Have tmf deforestation to agriculture at least once
  temp_est <- feglm(fml = as.formula(paste0("loss_commodity", " ~ 1 | ", est_parameters[["fe"]])),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    pre_d_clean_agri <- d[unlist(temp_est$obs_selection),]
  }  else { 
    pre_d_clean_agri <- d
  }
  
  rm(d)
  
  # Regressions of tmf_agri
  elm <- 1
  for(rfs_exp in est_parameters[["all_exposures_rfs"]]){#agri_crops eaear_mapmat[,"Crops"]
    rfs_1lag3lead_notrend_clt5[[CNT]][[elm]] <- make_main_reg(pre_process = TRUE,
                                                              pre_processed_data = pre_d_clean_agri, # NOTICE THIS
                                                              
                                                              outcome_variable = "loss_commodity", # AND THIS
                                                              continent = est_parameters[["continent"]],
                                                              start_year = est_parameters[["start_year"]],
                                                              end_year = est_parameters[["end_year"]],
                                                              
                                                              aggr_dyn = TRUE,
                                                              exposure_rfs = rfs_exp,
                                                              all_exposures_rfs = est_parameters[["all_exposures_rfs"]],
                                                              
                                                              original_rfs_treatments = c("statute_conv"),
                                                              rfs_lead = est_parameters[["leads"]], 
                                                              rfs_lag = est_parameters[["lags"]],
                                                              rfs_fya = est_parameters[["fya"]],  #ya,
                                                              rfs_pya = est_parameters[["pya"]], #ya,
                                                              
                                                              most_correlated_only = est_parameters[["most_correlated_only"]],
                                                              annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                              
                                                              s_trend = est_parameters[["s_trend"]],
                                                              s_trend_loga = est_parameters[["s_trend_loga"]],
                                                              fe = est_parameters[["fe"]], 
                                                              cluster_var1 = est_parameters[["cluster_var1"]], 
                                                              
                                                              rfs_rando = ""
    )
    
    names(rfs_1lag3lead_notrend_clt5[[CNT]])[elm] <- rfs_exp
    elm <- elm + 1
  }
  
}


## CLUSTER 20 
est_parameters[["cluster_var1"]] <- "grid_id_20"
rfs_1lag3lead_notrend_clt20 <- list(all = list(), 
                                    America = list(), 
                                    Africa = list(), 
                                    Asia = list())


for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters[["continent"]] <- CNT
  
  # - are in study area
  if(est_parameters[["continent"]] != "all"){
    d <- dplyr::filter(d_all, continent_name == est_parameters[["continent"]])
  } else {
    d <- d_all
  }
  
  # Have tmf deforestation to agriculture at least once
  temp_est <- feglm(fml = as.formula(paste0("loss_commodity", " ~ 1 | ", est_parameters[["fe"]])),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    pre_d_clean_agri <- d[unlist(temp_est$obs_selection),]
  }  else { 
    pre_d_clean_agri <- d
  }
  
  rm(d)
  
  # Regressions of tmf_agri
  elm <- 1
  for(rfs_exp in est_parameters[["all_exposures_rfs"]]){#agri_crops eaear_mapmat[,"Crops"]
    rfs_1lag3lead_notrend_clt20[[CNT]][[elm]] <- make_main_reg(pre_process = TRUE,
                                                               pre_processed_data = pre_d_clean_agri, # NOTICE THIS
                                                               
                                                               outcome_variable = "loss_commodity", # AND THIS
                                                               continent = est_parameters[["continent"]],
                                                               start_year = est_parameters[["start_year"]],
                                                               end_year = est_parameters[["end_year"]],
                                                               
                                                               aggr_dyn = TRUE,
                                                               exposure_rfs = rfs_exp,
                                                               all_exposures_rfs = est_parameters[["all_exposures_rfs"]],
                                                               
                                                               original_rfs_treatments = c("statute_conv"),
                                                               rfs_lead = est_parameters[["leads"]], 
                                                               rfs_lag = est_parameters[["lags"]],
                                                               rfs_fya = est_parameters[["fya"]],  #ya,
                                                               rfs_pya = est_parameters[["pya"]], #ya,
                                                               
                                                               most_correlated_only = est_parameters[["most_correlated_only"]],
                                                               annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                               
                                                               s_trend = est_parameters[["s_trend"]],
                                                               s_trend_loga = est_parameters[["s_trend_loga"]],
                                                               fe = est_parameters[["fe"]], 
                                                               cluster_var1 = est_parameters[["cluster_var1"]], 
                                                               
                                                               rfs_rando = ""
    )
    
    names(rfs_1lag3lead_notrend_clt20[[CNT]])[elm] <- rfs_exp
    elm <- elm + 1
  }
  
}

## TWO WAYS CLUSTERING (EXPOSURE LEVEL & GRID_10-YEAR)
est_parameters[["clustering"]] <- "exposure2ways"
est_parameters[["cluster_var2"]] <- "grid_id_10_year"

rfs_1lag3lead_notrend_clt2way <- list(all = list(), 
                                      America = list(), 
                                      Africa = list(), 
                                      Asia = list())


for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters[["continent"]] <- CNT
  
  # - are in study area
  if(est_parameters[["continent"]] != "all"){
    d <- dplyr::filter(d_all, continent_name == est_parameters[["continent"]])
  } else {
    d <- d_all
  }
  
  # Have tmf deforestation to agriculture at least once
  temp_est <- feglm(fml = as.formula(paste0("loss_commodity", " ~ 1 | ", est_parameters[["fe"]])),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    pre_d_clean_agri <- d[unlist(temp_est$obs_selection),]
  }  else { 
    pre_d_clean_agri <- d
  }
  
  rm(d)
  
  # Regressions of tmf_agri
  elm <- 1
  for(rfs_exp in est_parameters[["all_exposures_rfs"]]){#agri_crops eaear_mapmat[,"Crops"]
    rfs_1lag3lead_notrend_clt2way[[CNT]][[elm]] <- make_main_reg(pre_process = TRUE,
                                                                 pre_processed_data = pre_d_clean_agri, # NOTICE THIS
                                                                 
                                                                 outcome_variable = "loss_commodity", # AND THIS
                                                                 continent = est_parameters[["continent"]],
                                                                 start_year = est_parameters[["start_year"]],
                                                                 end_year = est_parameters[["end_year"]],
                                                                 
                                                                 aggr_dyn = TRUE,
                                                                 exposure_rfs = rfs_exp,
                                                                 all_exposures_rfs = est_parameters[["all_exposures_rfs"]],
                                                                 
                                                                 original_rfs_treatments = c("statute_conv"),
                                                                 rfs_lead = est_parameters[["leads"]], 
                                                                 rfs_lag = est_parameters[["lags"]],
                                                                 rfs_fya = est_parameters[["fya"]],  #ya,
                                                                 rfs_pya = est_parameters[["pya"]], #ya,
                                                                 
                                                                 most_correlated_only = est_parameters[["most_correlated_only"]],
                                                                 annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                                 
                                                                 s_trend = est_parameters[["s_trend"]],
                                                                 s_trend_loga = est_parameters[["s_trend_loga"]],
                                                                 fe = est_parameters[["fe"]], 
                                                                 
                                                                 clustering = est_parameters[["clustering"]], 
                                                                 cluster_var1 = est_parameters[["cluster_var1"]], 
                                                                 cluster_var2 = est_parameters[["cluster_var2"]], 
                                                                 
                                                                 rfs_rando = ""
    )
    
    names(rfs_1lag3lead_notrend_clt2way[[CNT]])[elm] <- rfs_exp
    elm <- elm + 1
  }
  
}



## because it's the end of the series of robustness checks on leads, set the parameter to default
est_parameters[["clustering"]] <- "oneway"
est_parameters[["cluster_var1"]] <- "grid_id_10" 
est_parameters[["cluster_var2"]] <- "grid_id_10"



#### Pasture control #### 
est_parameters[["control_pasture"]] <- TRUE
rfs_1lag3lead_notrend_clt10_pstctrl <- list(all = list(), 
                                            America = list(), 
                                            Africa = list(), 
                                            Asia = list())


for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters[["continent"]] <- CNT
  
  # - are in study area
  if(est_parameters[["continent"]] != "all"){
    d <- dplyr::filter(d_all, continent_name == est_parameters[["continent"]])
  } else {
    d <- d_all
  }
  
  # Have tmf deforestation to agriculture at least once
  temp_est <- feglm(fml = as.formula(paste0("loss_commodity", " ~ 1 | ", est_parameters[["fe"]])),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    pre_d_clean_agri <- d[unlist(temp_est$obs_selection),]
  }  else { 
    pre_d_clean_agri <- d
  }
  
  rm(d)
  
  # Regressions of tmf_agri
  elm <- 1
  for(rfs_exp in est_parameters[["all_exposures_rfs"]]){#agri_crops eaear_mapmat[,"Crops"]
    rfs_1lag3lead_notrend_clt10_pstctrl[[CNT]][[elm]] <- make_main_reg(pre_process = TRUE,
                                                                       pre_processed_data = pre_d_clean_agri, # NOTICE THIS
                                                                       
                                                                       outcome_variable = "loss_commodity", # AND THIS
                                                                       continent = est_parameters[["continent"]],
                                                                       start_year = est_parameters[["start_year"]],
                                                                       end_year = est_parameters[["end_year"]],
                                                                       
                                                                       aggr_dyn = TRUE,
                                                                       exposure_rfs = rfs_exp,
                                                                       all_exposures_rfs = est_parameters[["all_exposures_rfs"]],
                                                                       
                                                                       original_rfs_treatments = c("statute_conv"),
                                                                       rfs_lead = est_parameters[["leads"]], 
                                                                       rfs_lag = est_parameters[["lags"]],
                                                                       rfs_fya = est_parameters[["fya"]],  #ya,
                                                                       rfs_pya = est_parameters[["pya"]], #ya,
                                                                       
                                                                       most_correlated_only = est_parameters[["most_correlated_only"]],
                                                                       annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                                       
                                                                       control_pasture = est_parameters[["control_pasture"]],
                                                                       
                                                                       s_trend = est_parameters[["s_trend"]],
                                                                       s_trend_loga = est_parameters[["s_trend_loga"]],
                                                                       fe = est_parameters[["fe"]], 
                                                                       cluster_var1 = est_parameters[["cluster_var1"]], 
                                                                       
                                                                       rfs_rando = ""
    )
    
    names(rfs_1lag3lead_notrend_clt10_pstctrl[[CNT]])[elm] <- rfs_exp
    elm <- elm + 1
  }
  
}


est_parameters[["control_pasture"]] <- FALSE


saveRDS(rfs_1lag1lead_notrend_clt10, here("temp_data", "reg_results", "rfs", "result_lists", 
                                          "rfs_1lag1lead_notrend_clt10.Rdata"))
saveRDS(rfs_1lag2lead_notrend_clt10, here("temp_data", "reg_results", "rfs", "result_lists", 
                                          "rfs_1lag2lead_notrend_clt10.Rdata"))
saveRDS(rfs_1lag4lead_notrend_clt10, here("temp_data", "reg_results", "rfs", "result_lists", 
                                          "rfs_1lag4lead_notrend_clt10.Rdata"))
saveRDS(rfs_1lag5lead_notrend_clt10, here("temp_data", "reg_results", "rfs", "result_lists", 
                                          "rfs_1lag5lead_notrend_clt10.Rdata"))

saveRDS(rfs_1lag3lead_lag2ctrl_notrend_clt10, here("temp_data", "reg_results", "rfs", "result_lists", 
                                                   "rfs_1lag3lead_lag2ctrl_notrend_clt10.Rdata"))
saveRDS(rfs_2lag3lead_notrend_clt10, here("temp_data", "reg_results", "rfs", "result_lists", 
                                          "rfs_2lag3lead_notrend_clt10.Rdata"))

saveRDS(rfs_1lag3lead_lintrend_clt10, here("temp_data", "reg_results", "rfs", "result_lists", 
                                           "rfs_1lag3lead_lintrend_clt10.Rdata"))
saveRDS(rfs_1lag3lead_logtrend_clt10, here("temp_data", "reg_results", "rfs", "result_lists", 
                                           "rfs_1lag3lead_logtrend_clt10.Rdata"))

saveRDS(rfs_1lag3lead_notrend_clt5, here("temp_data", "reg_results", "rfs", "result_lists", 
                                         "rfs_1lag3lead_notrend_clt5.Rdata"))
saveRDS(rfs_1lag3lead_notrend_clt20, here("temp_data", "reg_results", "rfs", "result_lists", 
                                          "rfs_1lag3lead_notrend_clt20.Rdata"))
saveRDS(rfs_1lag3lead_notrend_clt2way, here("temp_data", "reg_results", "rfs", "result_lists", 
                                            "rfs_1lag3lead_notrend_clt2way.Rdata"))

saveRDS(rfs_1lag3lead_notrend_clt10_pstctrl, here("temp_data", "reg_results", "rfs", "result_lists", 
                                                  "rfs_1lag3lead_notrend_clt10_pstctrl.Rdata"))



#### PLOT MAIN FILTERED BY ROBUSTNESS #### 
# gather all robustness results in a list 
robu_list <- list(#rfs_1lag1lead_notrend_clt10, 
  rfs_1lag2lead_notrend_clt10, 
  rfs_1lag4lead_notrend_clt10, 
  #rfs_1lag5lead_notrend_clt10,
  rfs_1lag3lead_lag2ctrl_notrend_clt10,
  rfs_2lag3lead_notrend_clt10,
  rfs_1lag3lead_lintrend_clt10, 
  rfs_1lag3lead_logtrend_clt10, 
  rfs_1lag3lead_notrend_clt5, 
  rfs_1lag3lead_notrend_clt20, 
  rfs_1lag3lead_notrend_clt2way, 
  rfs_1lag3lead_notrend_clt10_pstctrl)

rm(df)
dyn_df_list <- list()

for(CNT in c("all", "America", "Africa", "Asia")){ # 1 and 2 correspond respectively to fya and pya, or aggrleads and aggrlags
  # This extracts the main results 
  dyn_df <- lapply(rfs_1lag3lead_notrend_clt10[[CNT]], FUN = function(x){as.data.frame(x)[1,] }) %>% bind_rows()
  dyn_df$term <- gsub(pattern = "_X_.*$", x = row.names(dyn_df), replacement = "") # replace everything after and including _X_ with nothing
  dyn_df$model <- CNT # dyn_des[dyn]
  
  # Now, extract estimates and p-values from every robustness exercise
  rob_no <- 1
  for(rc_res in robu_list){
    rob_df <- lapply(rc_res[[CNT]], FUN = function(x){as.data.frame(x)[1,] }) %>% bind_rows()
    
    # new column for the validity of crop estimates to this particular robustness check
    dyn_df[, paste0("rob_valid_", rob_no)] <- NA
    # fill it for every crop estimate
    for(crop_est in row.names(dyn_df)){
      dyn_df[crop_est, paste0("rob_valid_", rob_no)] <- (dyn_df[crop_est,"Estimate"]*rob_df[crop_est, "Estimate"] > 0 & # estimates should be of same signe
                                                           rob_df[crop_est, "Pr(>|t|)"] < 0.1) # AND p-value should be lower than 10%
    }
    
    rob_no <- rob_no + 1
  }
  rm(rob_df, rob_no, crop_est, rc_res)
  dyn_df_list[[CNT]] <- dyn_df
  
}

df <- bind_rows(dyn_df_list)
df[df$model == "all", "model"] <- "All"

df <- dplyr::mutate(df, rob_valid = rowMeans(across(.cols = contains("rob_valid_"), 
                                                    na.rm = FALSE)))

df$significant01 <- ""
df[df[,"Pr(>|t|)"] < 0.1, "significant01"] <- "p-value < .1"

df$highlight <- ""
df[df$significant01 == "p-value < .1" & df$rob_valid==1, "highlight"] <- "p-value < .1 & robust" 


names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"
row.names(df) <- NULL

#title <- "Impacts of RFS ramping up on tropical deforestation, via exposures to different land uses, 2001-2019"
title <- paste0("Cummulative effects of realized and expected biofuel mandates on tropical deforestation for major commercial crops, 2011-2019")

# IF ONLY SIGNIFICANT ONES ARE FEATURED
crop_groups <- list(c("Group 1", "Cereals", "Rice"),
                    c("Group 2", "Cotton", "Oil crops"),
                    c("Group 3", "Coconut", "Oil palm"),
                    c("Group 4", "Biomass crops", "Sugar crops"), # "Groundnut","Maize","Rapeseed","Rice","Sorghum","Soy","Sugar (cane or beet)","Sunflower",
                    c("Group 5", "Tobacco", "Citrus"), # "Citrus", "Cotton",
                    c("Group 6", "Banana", "Tea") # "Cocoa", "Coffee", "Rubber", "Cocoa or Coffee
)
{dwplot(df,
        dot_args = list(size = 2.5,aes(shape = model, alpha = highlight)),#, alpha = significant01
        whisker_args = list(size = 1, aes(alpha = highlight)),#alpha = significant01
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(#eaear_Barley = "Barley",
      #eaear_Groundnut = "Groundnut", 
      eaear_Maizegrain = "Maize",
      
      eaear_Fodder = "Pasture",
      
      eaear_Cereals = "Cereals",
      eaear_Rice = "Rice",
      
      eaear_Soy_compo = "Soy", 
      
      eaear_Cotton = "Cotton",
      eaear_Oilfeed_crops = "Oil crops",
      
      #eaear_Rapeseed = "Rapeseed", 
      # eaear_Sorghum2 = "Sorghum", 
      eaear_Biomass = "Biomass crops",
      eaear_Sugarcane = "Sugar crops", 
      #eaear_Sunflower = "Sunflower",
      #eaear_Wheat = "Wheat",
      
      eaear_Coconut = "Coconut", 
      eaear_Oilpalm = "Oil palm",
      
      eaear_Tobacco = "Tobacco", 
      eaear_Citrus = "Citrus",
      
      eaear_Banana = "Banana", 
      eaear_Cocoa_Coffee = "Cocoa or Coffee",
      #eaear_Cocoa  = "Cocoa", 
      #eaear_Coffee = "Coffee", 
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
          legend.title = element_blank()) } %>% 
  # don't run this if we used only limited crops
  add_brackets(crop_groups, face = "italic")


### LIMITED PLOT ### 
df_lmt <- df[df$term %in% limited_crops1,]
{dwplot(df_lmt,
        dot_args = list(size = 2.5,aes(shape = model, alpha = highlight)),#, alpha = significant01
        whisker_args = list(size = 1, aes(alpha = highlight)),#alpha = significant01
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(#eaear_Barley = "Barley",
      #eaear_Groundnut = "Groundnut", 
      eaear_Maizegrain = "Maize",
      
      eaear_Fodder = "Pasture",
      
      eaear_Cereals = "Cereals",
      eaear_Rice = "Rice",
      
      eaear_Soy_compo = "Soy", 
      
      eaear_Cotton = "Cotton",
      eaear_Oilfeed_crops = "Oil crops",
      
      #eaear_Rapeseed = "Rapeseed", 
      # eaear_Sorghum2 = "Sorghum", 
      eaear_Biomass = "Biomass crops",
      eaear_Sugarcane = "Sugar crops", 
      #eaear_Sunflower = "Sunflower",
      #eaear_Wheat = "Wheat",
      
      eaear_Coconut = "Coconut", 
      eaear_Oilpalm = "Oil palm",
      
      eaear_Tobacco = "Tobacco", 
      eaear_Citrus = "Citrus",
      
      eaear_Banana = "Banana", 
      eaear_Cocoa_Coffee = "Cocoa or Coffee",
      #eaear_Cocoa  = "Cocoa", 
      #eaear_Coffee = "Coffee", 
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




### PLOT AGGREGATED LEAD AND LAG EFFECTS ####
rm(dyn_df)
dyn_df_list <- list()
for(CNT in c("America", "Africa", "Asia")){
  dyn_df <- lapply(rfs_1lag3lead_notrend_clt10[[CNT]], FUN = function(x){as.data.frame(x)[c(1:3),]}) %>% bind_rows()
  
  dyn_df$term <- gsub(pattern = "_X_.*$", x = row.names(dyn_df), replacement = "") # replace everything after and including _X_ with nothing
  
  dyn_df$model <- "Via realized mediation"
  dyn_df$model[grepl("aggrleads", row.names(dyn_df))] <- "Via anticipations"
  dyn_df$model[grepl("aggrall", row.names(dyn_df))] <- "Cumulative effect"
  #dyn_df$model <- factor(dyn_df$model, levels = c("Via realized mediation", "Via anticipations", "Cumulative effects" ))
  
  dyn_df$continent <- CNT
  
  dyn_df_list[[CNT]] <- dyn_df
}

dyn_df <- bind_rows(dyn_df_list)

dyn_df$continent <- factor(dyn_df$continent, levels = c("America", "Africa", "Asia"))

dyn_df$significant01 <- if_else(dyn_df[,"Pr(>|t|)"] < 0.1, "p-value < .1", "")
# dyn_df <- dyn_df[dyn_df$significant01 != "", ]

names(dyn_df)[names(dyn_df)=="Estimate"] <- "estimate"
names(dyn_df)[names(dyn_df)=="Std. Error"] <- "std.error"
row.names(dyn_df) <- NULL

dyn_df <- left_join(dyn_df, df[,c("term", "model", "highlight")], 
                    by = c("term" = "term", "continent" = "model"))

# Limite to main crops of interest for conciseness
dyn_df_lmt <- dyn_df[dyn_df$term %in% limited_crops1,]

# adhoc remove africa, since there is nothing to highlight in current specification
dyn_df_lmt <- dplyr::filter(dyn_df_lmt, continent != "Africa")

{dwplot(dyn_df_lmt,
        dot_args = list(size = 2.5,aes(shape = model, col = model, alpha = highlight)),#,  shape = direction,
        whisker_args = list(size = 1, aes(col = model, alpha = highlight)),#alpha = significant01 linetype = model, 
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>%  
    relabel_predictors(c(#eaear_Barley = "Barley",
      #eaear_Groundnut = "Groundnut", 
      eaear_Maizegrain = "Maize",
      
      eaear_Fodder = "Pasture",
      
      eaear_Cereals = "Cereals",
      eaear_Rice = "Rice",
      
      eaear_Soy_compo = "Soy", 
      
      eaear_Cotton = "Cotton",
      eaear_Oilfeed_crops = "Oil crops",
      
      #eaear_Rapeseed = "Rapeseed", 
      # eaear_Sorghum2 = "Sorghum", 
      eaear_Biomass = "Biomass crops",
      eaear_Sugarcane = "Sugar crops", 
      #eaear_Sunflower = "Sunflower",
      #eaear_Wheat = "Wheat",
      
      eaear_Coconut = "Coconut", 
      eaear_Oilpalm = "Oil palm",
      
      eaear_Tobacco = "Tobacco", 
      eaear_Citrus = "Citrus",
      
      eaear_Banana = "Banana", 
      eaear_Cocoa_Coffee = "Cocoa or Coffee",
      #eaear_Cocoa  = "Cocoa", 
      #eaear_Coffee = "Coffee", 
      eaear_Rubber = "Rubber", 
      eaear_Tea = "Tea"
    )) +
    #scale_color_manual(values=wes_palette(n=4, name="Darjeeling1", type = "discrete")) +
    scale_color_brewer(type = "div", palette="Dark2") +
    guides(
      shape = guide_legend("", reverse = T), 
      colour = guide_legend("", reverse = T)
    ) +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    facet_grid(cols = vars(continent))+
    #ggtitle(title) +
    #labs(title = "", subtitle = "", caption = "", tag = CNT) +
    theme(plot.title = element_text(face="bold", size=c(10)),
          legend.position = c(0.85, 0.05),
          legend.justification = c(0, 0), 
          legend.background = element_rect(colour="grey80"),
          legend.title = element_blank()) }







#### PLOT AGGR LEAD AND LAGS WHEN t-2 IN CUMULATIVE EFFECTS ####
rm(dyn_df)
dyn_df_list <- list()
for(CNT in c("America", "Africa", "Asia")){
  dyn_df <- lapply(rfs_2lag3lead_notrend_clt10[[CNT]], FUN = function(x){as.data.frame(x)[c(1:3),]}) %>% bind_rows()
  
  dyn_df$term <- gsub(pattern = "_X_.*$", x = row.names(dyn_df), replacement = "") # replace everything after and including _X_ with nothing
  
  dyn_df$model <- "Via realized mediation"
  dyn_df$model[grepl("aggrleads", row.names(dyn_df))] <- "Via anticipations"
  dyn_df$model[grepl("aggrall", row.names(dyn_df))] <- "Cumulative effect"
  #dyn_df$model <- factor(dyn_df$model, levels = c("Via realized mediation", "Via anticipations", "Cumulative effects" ))
  
  dyn_df$continent <- CNT
  
  dyn_df_list[[CNT]] <- dyn_df
}

dyn_df <- bind_rows(dyn_df_list)

dyn_df$continent <- factor(dyn_df$continent, levels = c("America", "Africa", "Asia"))

dyn_df$significant01 <- if_else(dyn_df[,"Pr(>|t|)"] < 0.1, "p-value < .1", "")
# dyn_df <- dyn_df[dyn_df$significant01 != "", ]

names(dyn_df)[names(dyn_df)=="Estimate"] <- "estimate"
names(dyn_df)[names(dyn_df)=="Std. Error"] <- "std.error"
row.names(dyn_df) <- NULL

# dyn_df <- left_join(dyn_df, df[,c("term", "model", "highlight")], 
#                     by = c("term" = "term", "continent" = "model"))

# Limite to main crops of interest for conciseness
dyn_df_lmt <- dyn_df[dyn_df$term %in% limited_crops1,]

# adhoc remove africa, since there is nothing to highlight in current specification
dyn_df_lmt <- dplyr::filter(dyn_df_lmt, continent != "Africa")

{dwplot(dyn_df_lmt,
        dot_args = list(size = 2.5,aes(shape = model, col = model, alpha = significant01)),#,  shape = direction,
        whisker_args = list(size = 1, aes(col = model, alpha = significant01)),#alpha = significant01 linetype = model, 
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>%  
    relabel_predictors(c(#eaear_Barley = "Barley",
      #eaear_Groundnut = "Groundnut", 
      eaear_Maizegrain = "Maize",
      
      eaear_Fodder = "Pasture",
      
      eaear_Cereals = "Cereals",
      eaear_Rice = "Rice",
      
      eaear_Soy_compo = "Soy", 
      
      eaear_Cotton = "Cotton",
      eaear_Oilfeed_crops = "Oil crops",
      
      #eaear_Rapeseed = "Rapeseed", 
      # eaear_Sorghum2 = "Sorghum", 
      eaear_Biomass = "Biomass crops",
      eaear_Sugarcane = "Sugar crops", 
      #eaear_Sunflower = "Sunflower",
      #eaear_Wheat = "Wheat",
      
      eaear_Coconut = "Coconut", 
      eaear_Oilpalm = "Oil palm",
      
      eaear_Tobacco = "Tobacco", 
      eaear_Citrus = "Citrus",
      
      eaear_Banana = "Banana", 
      eaear_Cocoa_Coffee = "Cocoa or Coffee",
      #eaear_Cocoa  = "Cocoa", 
      #eaear_Coffee = "Coffee", 
      eaear_Rubber = "Rubber", 
      eaear_Tea = "Tea"
    )) +
    #scale_color_manual(values=wes_palette(n=4, name="Darjeeling1", type = "discrete")) +
    scale_color_brewer(type = "div", palette="Dark2") +
    guides(
      shape = guide_legend("", reverse = T), 
      colour = guide_legend("", reverse = T)
    ) +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    facet_grid(cols = vars(continent))+
    #ggtitle(title) +
    #labs(title = "", subtitle = "", caption = "", tag = CNT) +
    theme(plot.title = element_text(face="bold", size=c(10)),
          legend.position = c(0.85, 0.05),
          legend.justification = c(0, 0), 
          legend.background = element_rect(colour="grey80"),
          legend.title = element_blank()) }





#### JOINT ESTIMATION #### 

rfs_1lag3lead_notrend_clt10_joint <- list(all = list(), 
                                          America = list(), 
                                          Africa = list(), 
                                          Asia = list())

# prepare data set in advance, to repeat in all regressions
for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  elm <- 1
  rfs_1lag3lead_notrend_clt10_joint[[CNT]][[elm]] <- make_main_reg(pre_process = FALSE, # NOTICE THIS
                                                                   
                                                                   outcome_variable = "loss_commodity", 
                                                                   continent = CNT, # AND THIS
                                                                   start_year = est_parameters[["start_year"]],
                                                                   end_year = est_parameters[["end_year"]],
                                                                   
                                                                   aggr_dyn = TRUE,
                                                                   exposure_rfs = est_parameters[["all_exposures_rfs"]], # AND THIS
                                                                   all_exposures_rfs = est_parameters[["all_exposures_rfs"]],
                                                                   
                                                                   original_rfs_treatments = c("statute_conv"),
                                                                   rfs_lead = est_parameters[["leads"]], 
                                                                   rfs_lag = est_parameters[["lags"]],
                                                                   rfs_fya = est_parameters[["fya"]],  #ya,
                                                                   rfs_pya = est_parameters[["pya"]], #ya,
                                                                   
                                                                   control_all_absolute_rfs = FALSE, # AND THIS
                                                                   most_correlated_only = est_parameters[["most_correlated_only"]],
                                                                   annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                                   
                                                                   s_trend = est_parameters[["s_trend"]],
                                                                   s_trend_loga = est_parameters[["s_trend_loga"]],
                                                                   fe = est_parameters[["fe"]], 
                                                                   cluster_var1 = est_parameters[["cluster_var1"]], 
                                                                   
                                                                   rfs_rando = ""
  )
}    

### PLOT COEFFICIENTS ### 
# prepare regression outputs in a tidy data frame readable by dwplot

rm(df)
dyn_df_list <- list()

for(CNT in c("all", "America", "Africa", "Asia")){ # 1 and 2 correspond respectively to fya and pya, or aggrleads and aggrlags
  dyn_df <- lapply(rfs_1lag3lead_notrend_clt10_joint[[CNT]], FUN = function(x){x <- as.data.frame(x) 
  return(x[grepl("_aggrall", row.names(x)),])}) %>% bind_rows()
  dyn_df$term <- gsub(pattern = "_X_.*$", x = row.names(dyn_df), replacement = "") # replace everything after and including _X_ with nothing
  dyn_df$model <- CNT # dyn_des[dyn]
  dyn_df_list[[CNT]] <- dyn_df
}
df <- bind_rows(dyn_df_list)
df[df$model == "all", "model"] <- "Global"
#df$term <- rep(eaear_mapmat[,"Crops"], length(dyn_df_list))
df$significant01 <- ""
df[df[,"Pr(>|t|)"] < 0.1, "significant01"] <- "p-value < .1"

names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"
row.names(df) <- NULL

df_lmt <- df[df$term %in% limited_crops1,]
{dwplot(df,
        dot_args = list(size = 2.5,aes(shape = model, alpha = significant01)),#, alpha = significant01
        whisker_args = list(size = 1, aes(alpha = significant01)),#alpha = significant01
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(#eaear_Barley = "Barley",
      #eaear_Groundnut = "Groundnut", 
      eaear_Maizegrain = "Maize",
      
      eaear_Fodder = "Pasture",
      
      eaear_Cereals = "Cereals",
      eaear_Rice = "Rice",
      
      eaear_Soy_compo = "Soy", 
      
      eaear_Cotton = "Cotton",
      eaear_Oilfeed_crops = "Oil crops",
      
      #eaear_Rapeseed = "Rapeseed", 
      # eaear_Sorghum2 = "Sorghum", 
      eaear_Biomass = "Biomass crops",
      eaear_Sugarcane = "Sugar crops", 
      #eaear_Sunflower = "Sunflower",
      #eaear_Wheat = "Wheat",
      
      eaear_Coconut = "Coconut", 
      eaear_Oilpalm = "Oil palm",
      
      eaear_Tobacco = "Tobacco", 
      eaear_Citrus = "Citrus",
      
      eaear_Banana = "Banana", 
      eaear_Cocoa_Coffee = "Cocoa or Coffee",
      #eaear_Cocoa  = "Cocoa", 
      #eaear_Coffee = "Coffee", 
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




#### GAUSSIAN #### 
est_parameters[["distribution"]] <- "gaussian"
rfs_1lag3lead_notrend_clt10_gaussian <- list(all = list(), 
                                             America = list(), 
                                             Africa = list(), 
                                             Asia = list())


for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters[["continent"]] <- CNT
  
  elm <- 1
  for(rfs_exp in est_parameters[["all_exposures_rfs"]]){#agri_crops eaear_mapmat[,"Crops"]
    rfs_1lag3lead_notrend_clt10_gaussian[[CNT]][[elm]] <- make_main_reg(pre_process = FALSE,
                                                                        
                                                                        distribution = est_parameters[["distribution"]],
                                                                        
                                                                        outcome_variable = "loss_commodity", 
                                                                        continent = est_parameters[["continent"]],
                                                                        start_year = est_parameters[["start_year"]],
                                                                        end_year = est_parameters[["end_year"]],
                                                                        
                                                                        aggr_dyn = TRUE,
                                                                        exposure_rfs = rfs_exp,
                                                                        all_exposures_rfs = est_parameters[["all_exposures_rfs"]],
                                                                        
                                                                        original_rfs_treatments = c("statute_conv"),
                                                                        rfs_lead = est_parameters[["leads"]], 
                                                                        rfs_lag = est_parameters[["lags"]],
                                                                        rfs_fya = est_parameters[["fya"]],  #ya,
                                                                        rfs_pya = est_parameters[["pya"]], #ya,
                                                                        
                                                                        most_correlated_only = est_parameters[["most_correlated_only"]],
                                                                        annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                                        
                                                                        s_trend = est_parameters[["s_trend"]],
                                                                        s_trend_loga = est_parameters[["s_trend_loga"]],
                                                                        fe = est_parameters[["fe"]], 
                                                                        cluster_var1 = est_parameters[["cluster_var1"]], 
                                                                        
                                                                        rfs_rando = ""
    )
    
    names(rfs_1lag3lead_notrend_clt10_gaussian[[CNT]])[elm] <- rfs_exp
    elm <- elm + 1
  }
  
}
est_parameters[["distribution"]] <- "quasipoisson"

rm(df)
dyn_df_list <- list()

for(CNT in c("all", "America", "Africa", "Asia")){ # 1 and 2 correspond respectively to fya and pya, or aggrleads and aggrlags
  dyn_df <- lapply(rfs_1lag3lead_notrend_clt10_gaussian[[CNT]], FUN = function(x){as.data.frame(x)[1,] }) %>% bind_rows()
  dyn_df$term <- gsub(pattern = "_X_.*$", x = row.names(dyn_df), replacement = "") # replace everything after and including _X_ with nothing
  dyn_df$model <- CNT # dyn_des[dyn]
  dyn_df_list[[CNT]] <- dyn_df
}
df <- bind_rows(dyn_df_list)
df[df$model == "all", "model"] <- "Global"
#df$term <- rep(eaear_mapmat[,"Crops"], length(dyn_df_list))
df$significant01 <- ""
df[df[,"Pr(>|t|)"] < 0.1, "significant01"] <- "p-value < .1"

names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"
row.names(df) <- NULL


df_lmt <- df[df$term %in% limited_crops1,]
{dwplot(df_lmt,
        dot_args = list(size = 2.5,aes(shape = model, alpha = significant01)),#, alpha = significant01
        whisker_args = list(size = 1, aes(alpha = significant01)),#alpha = significant01
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(#eaear_Barley = "Barley",
      #eaear_Groundnut = "Groundnut", 
      eaear_Maizegrain = "Maize",
      
      eaear_Fodder = "Pasture",
      
      eaear_Cereals = "Cereals",
      eaear_Rice = "Rice",
      
      eaear_Soy_compo = "Soy", 
      
      eaear_Cotton = "Cotton",
      eaear_Oilfeed_crops = "Oil crops",
      
      #eaear_Rapeseed = "Rapeseed", 
      # eaear_Sorghum2 = "Sorghum", 
      eaear_Biomass = "Biomass crops",
      eaear_Sugarcane = "Sugar crops", 
      #eaear_Sunflower = "Sunflower",
      #eaear_Wheat = "Wheat",
      
      eaear_Coconut = "Coconut", 
      eaear_Oilpalm = "Oil palm",
      
      eaear_Tobacco = "Tobacco", 
      eaear_Citrus = "Citrus",
      
      eaear_Banana = "Banana", 
      eaear_Cocoa_Coffee = "Cocoa or Coffee",
      #eaear_Cocoa  = "Cocoa", 
      #eaear_Coffee = "Coffee", 
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





#### YEAR FE #### 
est_parameters[["fe"]] <- "grid_id + year"
rfs_1lag3lead_notrend_clt10_yearFE <- list(all = list(), 
                                           America = list(), 
                                           Africa = list(), 
                                           Asia = list())


for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters[["continent"]] <- CNT
  
  elm <- 1
  for(rfs_exp in est_parameters[["all_exposures_rfs"]]){#agri_crops eaear_mapmat[,"Crops"]
    rfs_1lag3lead_notrend_clt10_yearFE[[CNT]][[elm]] <- make_main_reg(pre_process = FALSE,
                                                                      
                                                                      outcome_variable = "loss_commodity", 
                                                                      continent = est_parameters[["continent"]],
                                                                      start_year = est_parameters[["start_year"]],
                                                                      end_year = est_parameters[["end_year"]],
                                                                      
                                                                      aggr_dyn = TRUE,
                                                                      exposure_rfs = rfs_exp,
                                                                      all_exposures_rfs = est_parameters[["all_exposures_rfs"]],
                                                                      
                                                                      original_rfs_treatments = c("statute_conv"),
                                                                      rfs_lead = est_parameters[["leads"]], 
                                                                      rfs_lag = est_parameters[["lags"]],
                                                                      rfs_fya = est_parameters[["fya"]],  #ya,
                                                                      rfs_pya = est_parameters[["pya"]], #ya,
                                                                      
                                                                      most_correlated_only = est_parameters[["most_correlated_only"]],
                                                                      annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                                      
                                                                      s_trend = est_parameters[["s_trend"]],
                                                                      s_trend_loga = est_parameters[["s_trend_loga"]],
                                                                      fe = est_parameters[["fe"]], 
                                                                      cluster_var1 = est_parameters[["cluster_var1"]], 
                                                                      
                                                                      rfs_rando = ""
    )
    
    names(rfs_1lag3lead_notrend_clt10_yearFE[[CNT]])[elm] <- rfs_exp
    elm <- elm + 1
  }
  
}
est_parameters[["fe"]] <- "grid_id + country_year"

rm(df)
dyn_df_list <- list()

for(CNT in c("all", "America", "Africa", "Asia")){ # 1 and 2 correspond respectively to fya and pya, or aggrleads and aggrlags
  dyn_df <- lapply(rfs_1lag3lead_notrend_clt10_yearFE[[CNT]], FUN = function(x){as.data.frame(x)[1,] }) %>% bind_rows()
  dyn_df$term <- gsub(pattern = "_X_.*$", x = row.names(dyn_df), replacement = "") # replace everything after and including _X_ with nothing
  dyn_df$model <- CNT # dyn_des[dyn]
  dyn_df_list[[CNT]] <- dyn_df
}
df <- bind_rows(dyn_df_list)
df[df$model == "all", "model"] <- "Global"
#df$term <- rep(eaear_mapmat[,"Crops"], length(dyn_df_list))
df$significant01 <- ""
df[df[,"Pr(>|t|)"] < 0.1, "significant01"] <- "p-value < .1"

names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"
row.names(df) <- NULL


df_lmt <- df[df$term %in% limited_crops1,]
{dwplot(df_lmt,
        dot_args = list(size = 2.5,aes(shape = model, alpha = significant01)),#, alpha = significant01
        whisker_args = list(size = 1, aes(alpha = significant01)),#alpha = significant01
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(#eaear_Barley = "Barley",
      #eaear_Groundnut = "Groundnut", 
      eaear_Maizegrain = "Maize",
      
      eaear_Fodder = "Pasture",
      
      eaear_Cereals = "Cereals",
      eaear_Rice = "Rice",
      
      eaear_Soy_compo = "Soy", 
      
      eaear_Cotton = "Cotton",
      eaear_Oilfeed_crops = "Oil crops",
      
      #eaear_Rapeseed = "Rapeseed", 
      # eaear_Sorghum2 = "Sorghum", 
      eaear_Biomass = "Biomass crops",
      eaear_Sugarcane = "Sugar crops", 
      #eaear_Sunflower = "Sunflower",
      #eaear_Wheat = "Wheat",
      
      eaear_Coconut = "Coconut", 
      eaear_Oilpalm = "Oil palm",
      
      eaear_Tobacco = "Tobacco", 
      eaear_Citrus = "Citrus",
      
      eaear_Banana = "Banana", 
      eaear_Cocoa_Coffee = "Cocoa or Coffee",
      #eaear_Cocoa  = "Cocoa", 
      #eaear_Coffee = "Coffee", 
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

#### DYNAMICS DESCRIPTIVES #### 
## RFS time series (un-interacted instrument) ##
rfs <- dplyr::arrange(rfs, year)
ts_rfs <- rfs[rfs$year>=2008,c("statute_conv")]
plot(ts_rfs)
lines(lowess(ts_rfs), col = 3)

# Dickey Fuller tests
adf_rfs <- ur.df(ts_rfs, 
                 type = "trend")

m3 <- cbind(
  t(adf_rfs@teststat),
  adf_rfs@cval)

# INTERPRETTION.
# First, note that in the model tested, lag = 0 means that there IS a unit root (see Wikipedia intuition section). 
# While trend and intercept terms = 0 reflect there is NO trend or intercept. 

# phi2 is the statistic for the null that all three of trend, intercept, AND lag are (jointly, it is) null (= 0) 
# Rejecting the null, i.e. phi2 > critical value means AT LEAST one of the three is non zero. (Does not say if it's the trend/intercept, or the lag). 
# Not rejecting it means that possibly there is a unit root (lag=0) and no trend and no intercept. --> Check phi3
# phi3: rejecting the null, i.e. phi3 > cval means AT LEAST one of drift OR lag is non zero. 
# Not rejecting means that possibly there is a unit root (lag=0) AND no intercept. --> Check tau3
# tau3: rejecting the null that lag = 0 means that there ISN'T a unit root. 
# Not rejecting it means there may be a unit root. 

# Note from Wikipedia that "As this test is asymmetrical, we are only concerned with negative values of our test statistic. 
# If the calculated test statistic is less (more negative) than the critical value, then the null hypothesis is rejected and no unit root is present."

# Here, all tests fail to reject the null. This is most likely due to the length of the time series that underpowers the tests. 
# If well powered, this could mean that there is no trend nor intercept, but there is a unit root. 

# That there might be a unit root is confirmed by the least specified test: 
adf_rfs <- ur.df(ts_rfs, 
                 type = "none")

m1 <- cbind(
  t(adf_rfs@teststat),
  adf_rfs@cval)

# However, the presence of a unit root and the absence of an intercept can be rejected at 5% confidence in this specification
# meaning that once one does not control for a trend, a drift appears in the data. But specifying a trend, the trend is not significant. 
# So we would say the drift-only is the best specification. We stop here. And once we stop at the drift, the presence of a unit root can be rejected. 
# Meaning that the data is stationary, and the only pattern that should be specified is an intercept. 
adf_rfs <- ur.df(ts_rfs, 
                 type = "drift")
m2 <- cbind(
  t(adf_rfs@teststat),
  adf_rfs@cval)


### PLOT ### 
adf_mat <- rbind(m3, m2, m1)
colnames(adf_mat) <- NULL

options(knitr.table.format = "latex")
kable(adf_mat, booktabs = T, align = "r",
      caption = "Augmented Dickey-Fuller tests of the RFS statutory mandates on conventional biofuels, 2008-2022") %>%
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Statistic" = 1,
                     "0.01" = 1,
                     "0.05" = 1,
                     "0.1" = 1),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 2,
                     "Critical values" = 3),
                   align = "c",
                   strikeout = F) %>%
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(adf_mat))),
              width = "4em",
              latex_valign = "c")

# On KPSS: 
# Here, the null hypothesis IS stationarity. the drift only is type "mu", and the drift + trend is "tau". 
# Use 1 lag as the default in ADF above
# "A major disadvantage for the KPSS test is that it has a high rate of Type I errors (it tends to reject the null hypothesis too often)." Combine with ADF 
ur.kpss(ts_rfs, type = "mu", use.lag = 1) %>% summary()

ur.kpss(diff(ts_rfs, differences = 1), 
        type = "mu", use.lag = 1) %>% summary() 

ur.kpss(diff(ts_rfs, differences = 2), 
        type = "mu", use.lag = 1) %>% summary() # --> second differencing necessary

# Only once the TS is differenced twice, we cannot reject that it is stationary 
# similar when specifying a trend in the deterministic component. 

## Outcome time series (summed within years) ##
cropexposed <- main_data[main_data$year >= 2008,]
cropexposed <- dplyr::arrange(cropexposed, year)
ts_defo <- ddply(cropexposed, "year", summarise, 
                 pantrop_defo = sum(loss_commodity, na.rm = TRUE))

plot(ts_defo)
lines(lowess(ts_defo), col = 3)

# Dickey Fuller tests
adf_defo <- ur.df(ts_defo$pantrop_defo, 
                  type = "trend")
cbind(
  t(adf_defo@teststat),
  adf_defo@cval)

adf_defo <- ur.df(ts_defo$pantrop_defo, 
                  type = "drift")
cbind(
  t(adf_defo@teststat),
  adf_defo@cval)

adf_defo <- ur.df(ts_defo$pantrop_defo, 
                  type = "none")
cbind(
  t(adf_defo@teststat),
  adf_defo@cval)

# what the tests tell is that we cannot reject that there is a unit root nor that there is no trend nor intercept. 
# So there might be a unit root and no trend nor intercept, but failure to reject is most probably due to underpower. 



### Fisher-type panel unit root test ### 
CNT <- "America"
for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters[["continent"]] <- CNT
  
  d <- main_data
  
  # Keep only in data the useful variables 
  d <- dplyr::select(d, all_of(unique(c("grid_id", "year", "lat", "lon", "continent_name", "country_name", "country_year",
                                        # "fc_2000", "fc_2009", "remaining_fc", # "accu_defo_since2k",
                                        "grid_id_5", "grid_id_10", "grid_id_20", "grid_id_5_year", "grid_id_10_year", "grid_id_20_year",
                                        # "tmf_agri", "tmf_flood", "tmf_plantation", "tmf_deforestation",
                                        "loss_commodity",
                                        "pasture_share_2000",
                                        eaear_mapmat[,"Crops"] ))))
  
  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", unique(c(all_rfs_treatments)))], by = c("year"))
  
  
  # - are in study period
  # if(start_year != 2011 | end_year != 2019){
  d <- dplyr::filter(d, year >= est_parameters[["start_year"]])
  d <- dplyr::filter(d, year <= est_parameters[["end_year"]])
  # }
  
  # - are in study area
  if(est_parameters[["continent"]] != "all"){
    d <- dplyr::filter(d, continent_name == est_parameters[["continent"]])
  }
  
  used_vars <- unique(c("grid_id", "year", "lat", "lon","continent_name",  "country_year",  #"country_name",  "remaining_fc", "accu_defo_since2k", # "sj_year",
                        "grid_id_5", "grid_id_10", "grid_id_20", "grid_id_5_year", "grid_id_10_year", "grid_id_20_year",
                        "loss_commodity",# "tmf_agri", "tmf_flood", "tmf_plantation",
                        eaear_mapmat[,"Crops"], "statute_conv"))    
  
  
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
  
  # Have tmf deforestation to agriculture at least once
  temp_est <- feglm(fml = as.formula(paste0("loss_commodity", " ~ 1 | ", est_parameters[["fe"]])),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    d_clean <- d[unlist(temp_est$obs_selection),]
  }  else { 
    d_clean <- d
  }
  
  # make one "interacted instrument", i.e. one of our regressors. 
  # note that the results of the tests are independent of the crop chosen for exposure.
  d_clean <- dplyr::mutate(d_clean, eaear_Soy_compo_X_statute_conv = eaear_Soy_compo*statute_conv)
  
  ### Fisher-type panel unit root test  
  # We need to convert the data to a pdata.frame format, to use this interface in purtest 
  pd_clean <- pdata.frame(d_clean[,c("grid_id", "year", "loss_commodity", "eaear_Soy_compo_X_statute_conv")], index = c("grid_id", "year"))
  # wide_d_clean <- reshape(pd_clean, direction = "wide", v.names = c(outcome_variable, regressors[1]),
  #                         timevar = "grid_id", idvar = "year", times = "grid_id")
  anyNA(pd_clean)
  purt <- list()
  correction <- "1"
  for(correction in c("0", "1", "trend")){
    purt[[correction]] <- plm::purtest(data = pd_clean,
                                       object = as.formula(paste0("eaear_Soy_compo_X_statute_conv ~ ", correction)),
                                       #object = pd_clean$eaear_Soy_compo_X_statute_conv, 
                                       pmax = 1, # length(unique(d_clean$year)) - 1,
                                       exo = "trend",
                                       test = "invnormal", 
                                       na.rm = TRUE)
  }
  
}









#### TMF PLANTATIONS #### 
est_parameters[["outcome_variable"]] <- "tmf_plantation"
rfs_TMF_1lag3lead_notrend_clt10 <- list(all = list(), 
                                        America = list(), 
                                        Africa = list(), 
                                        Asia = list())

d_all  <- readRDS(here("temp_data", "merged_datasets", "tmf_aoi", "tmf_aeay_pantrop_long_final_1990_2020.Rdata"))
# release some memory upfront
d_all <- dplyr::filter(d_all, year >= 2008, year <= 2019)

d_all <- dplyr::mutate(d_all, tmf_deforestation = tmf_agri + tmf_plantation)

# Keep only in data the useful variables 
d_all <- dplyr::select(d_all, all_of(unique(c("grid_id", "year", "lat", "lon", "continent_name", "country_name", "country_year",
                                              # "fc_2000", "fc_2009", "remaining_fc", # "accu_defo_since2k",
                                              "grid_id_5", "grid_id_10", "grid_id_20", "grid_id_5_year", "grid_id_10_year", "grid_id_20_year",
                                              "tmf_agri", "tmf_flood", "tmf_plantation",
                                              "pasture_share_2000",
                                              est_parameters[["all_exposures_rfs"]] )))) #sj, 

# Merge only the prices needed, not the whole price dataframe
d_all <- left_join(d_all, prices[,c("year", unique(c(all_rfs_treatments)))], by = c("year"))

# - are suitable to crop j 
if(est_parameters[["sjpos"]]){
  d_all <- dplyr::filter(d_all, !!as.symbol(exposure_rfs) > 0)
}

# - are in study period
d_all <- dplyr::filter(d_all, year >= est_parameters[["start_year"]])
d_all <- dplyr::filter(d_all, year <= est_parameters[["end_year"]])


for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters[["continent"]] <- CNT
  
  # - are in study area
  if(est_parameters[["continent"]] != "all"){
    d <- dplyr::filter(d_all, continent_name == est_parameters[["continent"]])
  } else {
    d <- d_all
  }
  
  # Have tmf deforestation to agriculture at least once
  temp_est <- feglm(fml = as.formula(paste0(est_parameters[["outcome_variable"]], " ~ 1 | ", est_parameters[["fe"]])),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    pre_d_clean_agri <- d[unlist(temp_est$obs_selection),]
  }  else { 
    pre_d_clean_agri <- d
  }
  
  rm(d)
  
  elm <- 1
  for(rfs_exp in plantation_crops){#agri_crops eaear_mapmat[,"Crops"]
    rfs_TMF_1lag3lead_notrend_clt10[[CNT]][[elm]] <- make_main_reg(pre_process = TRUE,
                                                                   pre_processed_data = pre_d_clean_agri, # NOTICE THIS
                                                                   outcome_variable = est_parameters[["outcome_variable"]], 
                                                                   
                                                                   continent = est_parameters[["continent"]],
                                                                   start_year = est_parameters[["start_year"]],
                                                                   end_year = est_parameters[["end_year"]],
                                                                   
                                                                   aggr_dyn = TRUE,
                                                                   exposure_rfs = rfs_exp,
                                                                   all_exposures_rfs = est_parameters[["all_exposures_rfs"]],
                                                                   
                                                                   original_rfs_treatments = c("statute_conv"),
                                                                   rfs_lead = est_parameters[["leads"]], 
                                                                   rfs_lag = est_parameters[["lags"]],
                                                                   rfs_fya = est_parameters[["fya"]],  #ya,
                                                                   rfs_pya = est_parameters[["pya"]], #ya,
                                                                   
                                                                   most_correlated_only = est_parameters[["most_correlated_only"]],
                                                                   annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                                   
                                                                   s_trend = est_parameters[["s_trend"]],
                                                                   s_trend_loga = est_parameters[["s_trend_loga"]],
                                                                   fe = est_parameters[["fe"]], 
                                                                   cluster_var1 = est_parameters[["cluster_var1"]], 
                                                                   
                                                                   rfs_rando = ""
    )
    
    names(rfs_TMF_1lag3lead_notrend_clt10[[CNT]])[elm] <- rfs_exp
    elm <- elm + 1
  }
  
}
est_parameters[["outcome_variable"]] <- "loss_commodity"




#### TMF WORLD #### 
rfs_tmf_2lag3lead_notrend_joint_clt10 <- list(all = list(), 
                                              America = list(), 
                                              Africa = list(), 
                                              Asia = list())
for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters <- list(outcome_variable  = "tmf_agri",
                         continent = CNT,
                         start_year = 2010, 
                         end_year = 2019, 
                         most_correlated_only = FALSE,
                         annual_rfs_controls = FALSE,
                         sjpos = FALSE,
                         lags = 2,
                         leads = 3,
                         fya = 0, 
                         pya = 0,
                         s_trend = FALSE, 
                         s_trend_loga = FALSE,
                         fe = "grid_id + country_year",
                         clustering = "oneway",
                         cluster_var1 = "grid_id_10",
                         cluster_var2 = "grid_id_10",
                         distribution = "quasipoisson")
  
  rfs_tmf_2lag3lead_notrend_joint_clt10[[CNT]][["tmf_agri"]] <- make_main_reg(pre_process = FALSE,# NOTICE THIS
                                                                              
                                                                              outcome_variable = est_parameters[["outcome_variable"]],
                                                                              continent = est_parameters[["continent"]],
                                                                              start_year = est_parameters[["start_year"]],
                                                                              end_year = est_parameters[["end_year"]],
                                                                              
                                                                              aggr_dyn = TRUE,
                                                                              exposure_rfs = agri_crops, # AND THIS
                                                                              control_all_absolute_rfs = FALSE, # AND THIS
                                                                              
                                                                              original_rfs_treatments = c("statute_conv"),
                                                                              rfs_lead = est_parameters[["leads"]], 
                                                                              rfs_lag = est_parameters[["lags"]],
                                                                              rfs_fya = est_parameters[["fya"]],  #ya,
                                                                              rfs_pya = est_parameters[["pya"]], #ya,
                                                                              
                                                                              most_correlated_only = est_parameters[["most_correlated_only"]],
                                                                              annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                                              
                                                                              s_trend = est_parameters[["s_trend"]],
                                                                              s_trend_loga = est_parameters[["s_trend_loga"]],
                                                                              fe = est_parameters[["fe"]], 
                                                                              clustering = est_parameters[["clustering"]],
                                                                              cluster_var1 = est_parameters[["cluster_var1"]], 
                                                                              cluster_var2 = est_parameters[["cluster_var2"]], 
                                                                              rfs_rando = ""
  )
  
  rfs_tmf_2lag3lead_notrend_joint_clt10[[CNT]][["tmf_plantation"]] <- make_main_reg(pre_process = FALSE,# NOTICE THIS
                                                                                    
                                                                                    outcome_variable = est_parameters[["outcome_variable"]],
                                                                                    continent = est_parameters[["continent"]],
                                                                                    start_year = est_parameters[["start_year"]],
                                                                                    end_year = est_parameters[["end_year"]],
                                                                                    
                                                                                    aggr_dyn = TRUE,
                                                                                    exposure_rfs = plantation_crops, # AND THIS
                                                                                    control_all_absolute_rfs = FALSE, # AND THIS
                                                                                    
                                                                                    original_rfs_treatments = c("statute_conv"),
                                                                                    rfs_lead = est_parameters[["leads"]], 
                                                                                    rfs_lag = est_parameters[["lags"]],
                                                                                    rfs_fya = est_parameters[["fya"]],  #ya,
                                                                                    rfs_pya = est_parameters[["pya"]], #ya,
                                                                                    
                                                                                    most_correlated_only = est_parameters[["most_correlated_only"]],
                                                                                    annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                                                    
                                                                                    s_trend = est_parameters[["s_trend"]],
                                                                                    s_trend_loga = est_parameters[["s_trend_loga"]],
                                                                                    fe = est_parameters[["fe"]], 
                                                                                    clustering = est_parameters[["clustering"]],
                                                                                    cluster_var1 = est_parameters[["cluster_var1"]], 
                                                                                    cluster_var2 = est_parameters[["cluster_var2"]], 
                                                                                    rfs_rando = ""
  )
}

rm(df)
dyn_df_list <- list()

for(CNT in c("all", "America", "Africa", "Asia")){ # 1 and 2 correspond respectively to fya and pya, or aggrleads and aggrlags
  dyn_df <- lapply(rfs_tmf_2lag3lead_notrend_joint_clt10[[CNT]], FUN = function(x){x <- as.data.frame(x) 
  return(x[grepl("_aggrall", row.names(x)),])}) %>% bind_rows()
  dyn_df$term <- gsub(pattern = "_X_.*$", x = row.names(dyn_df), replacement = "") # replace everything after and including _X_ with nothing
  dyn_df$model <- CNT # dyn_des[dyn]
  dyn_df_list[[CNT]] <- dyn_df
}
df <- bind_rows(dyn_df_list)
df[df$model == "all", "model"] <- "Global"
#df$term <- rep(eaear_mapmat[,"Crops"], length(dyn_df_list))
df$significant01 <- ""
df[df[,"Pr(>|t|)"] < 0.1, "significant01"] <- "p-value < .1"
# df[df[,"Pr(>|t|)"] < 0.05, "significant01"] <- "p-value < .05"
# df[df[,"Pr(>|t|)"] < 0.01, "significant01"] <- "p-value < .01"
# df$significant01 <- factor(df$significant01, levels = c("", "p-value < .1", "p-value < .05", "p-value < .01"))
#df <- df[df$significant01 != "", ]
#df$term <- sub(pattern = ".+?(_)", x = row.names(df), replacement = "") # replace everyting before the first underscore with nothing

names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"
row.names(df) <- NULL

#title <- "Impacts of RFS ramping up on tropical deforestation, via exposures to different land uses, 2001-2019"
title <- paste0("Cummulative effects of realized and expected biofuel mandates on tropical deforestation for major commercial crops, 2011-2019")
# If we want to add brackets on y axis to group k commodities. But not necessarily relevant, as some crops as in several categories. 
# A list of brackets; each element of the list should be a character vector consisting of 
# (1) a label for the bracket, (2) the name of the topmost variable to be enclosed by the bracket, 
# and (3) the name of the bottom most variable to be enclosed by the bracket.
# IF ALL CROPS ARE FEATURED: 
# crop_groups <- list(c("Group 1", "Cereals", "Sugar crops"), # "Groundnut","Maize","Rapeseed","Rice","Sorghum","Soy","Sugar (cane or beet)","Sunflower",
#                     c("Group 2", "Pasture", "Tobacco"), # "Citrus", "Cotton", 
#                     c("Group 3", "Coconut", "Oil palm"), 
#                     c("Group 4", "Banana", "Tea") # "Cocoa", "Coffee", "Rubber", 
# )
# IF ONLY SIGNIFICANT ONES ARE FEATURED
crop_groups <- list(c("Group 1", "Biomass crops", "Sugar crops"), # "Groundnut","Maize","Rapeseed","Rice","Sorghum","Soy","Sugar (cane or beet)","Sunflower",
                    c("Group 2", "Coconut", "Oil palm"),
                    c("Group 3", "Citrus", "Tobacco"), # "Citrus", "Cotton",
                    c("Group 4", "Banana", "Tea") # "Cocoa", "Coffee", "Rubber", "Cocoa or Coffee
)
{dwplot(df,
        dot_args = list(size = 2.5,aes(shape = model, alpha = significant01)),#, alpha = significant01
        whisker_args = list(size = 1, aes(alpha = significant01)),#alpha = significant01
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(#eaear_Barley = "Barley",
      #eaear_Groundnut = "Groundnut", 
      eaear_Maizegrain = "Maize",
      
      eaear_Fodder = "Pasture",
      
      eaear_Biomass = "Biomass crops",
      eaear_Cotton = "Cotton",
      eaear_Cereals = "Cereals",
      eaear_Oilfeed_crops = "Oil crops",
      #eaear_Rapeseed = "Rapeseed", 
      eaear_Rice = "Rice",
      eaear_Soy_compo = "Soy", 
      # eaear_Sorghum2 = "Sorghum", 
      eaear_Sugarcane = "Sugar crops", 
      #eaear_Sunflower = "Sunflower",
      #eaear_Wheat = "Wheat",
      
      eaear_Coconut = "Coconut", 
      eaear_Oilpalm = "Oil palm",
      
      eaear_Citrus = "Citrus",
      eaear_Tobacco = "Tobacco", 
      
      eaear_Banana = "Banana", 
      eaear_Cocoa_Coffee = "Cocoa or Coffee",
      #eaear_Cocoa  = "Cocoa", 
      #eaear_Coffee = "Coffee", 
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
          legend.title = element_blank()) } %>%  add_brackets(crop_groups, face = "italic")

#### joint estimation #### 
rfs_1lag3lead_notrend_joint <- list(all = list(), 
                                    America = list(), 
                                    Africa = list(), 
                                    Asia = list())
elm <- 1
for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters <- list(outcome_variable  = "loss_commodity",
                         continent = CNT,
                         start_year = 2010, 
                         end_year = 2019, 
                         most_correlated_only = FALSE,
                         annual_rfs_controls = FALSE,
                         sjpos = FALSE,
                         lags = 1,
                         leads = 3,
                         fya = 0, 
                         pya = 0,
                         s_trend = FALSE, 
                         s_trend_loga = FALSE,
                         fe = "grid_id + country_year",
                         clustering = "oneway",
                         cluster_var1 = "grid_id_10",
                         cluster_var2 = "grid_id_10",
                         distribution = "quasipoisson")
  
  rfs_1lag3lead_notrend_joint[[CNT]][[elm]] <- make_main_reg(pre_process = FALSE,# NOTICE THIS
                                                             
                                                             outcome_variable = est_parameters[["outcome_variable"]],
                                                             continent = est_parameters[["continent"]],
                                                             start_year = est_parameters[["start_year"]],
                                                             end_year = est_parameters[["end_year"]],
                                                             
                                                             aggr_dyn = TRUE,
                                                             exposure_rfs = eaear_mapmat[,"Crops"], # AND THIS
                                                             control_all_absolute_rfs = FALSE, # AND THIS
                                                             
                                                             original_rfs_treatments = c("statute_conv"),
                                                             rfs_lead = est_parameters[["leads"]], 
                                                             rfs_lag = est_parameters[["lags"]],
                                                             rfs_fya = est_parameters[["fya"]],  #ya,
                                                             rfs_pya = est_parameters[["pya"]], #ya,
                                                             
                                                             most_correlated_only = est_parameters[["most_correlated_only"]],
                                                             annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                             
                                                             s_trend = est_parameters[["s_trend"]],
                                                             s_trend_loga = est_parameters[["s_trend_loga"]],
                                                             fe = est_parameters[["fe"]], 
                                                             clustering = est_parameters[["clustering"]],
                                                             cluster_var1 = est_parameters[["cluster_var1"]], 
                                                             cluster_var2 = est_parameters[["cluster_var2"]], 
                                                             rfs_rando = ""
  )
  
  
}



rm(df)
dyn_df_list <- list()

for(CNT in c("all", "America", "Africa", "Asia")){ # 1 and 2 correspond respectively to fya and pya, or aggrleads and aggrlags
  dyn_df <- lapply(rfs_1lag3lead_notrend_joint[[CNT]], FUN = function(x){x <- as.data.frame(x) 
  return(x[grepl("_aggrall", row.names(x)),])}) %>% bind_rows()
  dyn_df$term <- gsub(pattern = "_X_.*$", x = row.names(dyn_df), replacement = "") # replace everything after and including _X_ with nothing
  dyn_df$model <- CNT # dyn_des[dyn]
  dyn_df_list[[CNT]] <- dyn_df
}
df <- bind_rows(dyn_df_list)
df[df$model == "all", "model"] <- "Global"
#df$term <- rep(eaear_mapmat[,"Crops"], length(dyn_df_list))
df$significant01 <- ""
df[df[,"Pr(>|t|)"] < 0.1, "significant01"] <- "p-value < .1"
# df[df[,"Pr(>|t|)"] < 0.05, "significant01"] <- "p-value < .05"
# df[df[,"Pr(>|t|)"] < 0.01, "significant01"] <- "p-value < .01"
# df$significant01 <- factor(df$significant01, levels = c("", "p-value < .1", "p-value < .05", "p-value < .01"))
#df <- df[df$significant01 != "", ]
#df$term <- sub(pattern = ".+?(_)", x = row.names(df), replacement = "") # replace everyting before the first underscore with nothing

names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"
row.names(df) <- NULL

#title <- "Impacts of RFS ramping up on tropical deforestation, via exposures to different land uses, 2001-2019"
title <- paste0("Cummulative effects of realized and expected biofuel mandates on tropical deforestation for major commercial crops, 2011-2019")

# IF ONLY SIGNIFICANT ONES ARE FEATURED
crop_groups <- list(c("Group 1", "Cereals", "Rice"),
                    c("Group 2", "Cotton", "Oil crops"),
                    c("Group 3", "Coconut", "Oil palm"),
                    c("Group 4", "Biomass crops", "Sugar crops"), # "Groundnut","Maize","Rapeseed","Rice","Sorghum","Soy","Sugar (cane or beet)","Sunflower",
                    c("Group 5", "Tobacco", "Citrus"), # "Citrus", "Cotton",
                    c("Group 6", "Banana", "Tea") # "Cocoa", "Coffee", "Rubber", "Cocoa or Coffee
)
{dwplot(df,
        dot_args = list(size = 2.5,aes(shape = model, alpha = significant01)),#, alpha = significant01
        whisker_args = list(size = 1, aes(alpha = significant01)),#alpha = significant01
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(#eaear_Barley = "Barley",
      #eaear_Groundnut = "Groundnut", 
      eaear_Maizegrain = "Maize",
      
      eaear_Fodder = "Pasture",
      
      eaear_Cereals = "Cereals",
      eaear_Rice = "Rice",
      
      eaear_Soy_compo = "Soy", 
      
      eaear_Cotton = "Cotton",
      eaear_Oilfeed_crops = "Oil crops",
      
      #eaear_Rapeseed = "Rapeseed", 
      # eaear_Sorghum2 = "Sorghum", 
      eaear_Biomass = "Biomass crops",
      eaear_Sugarcane = "Sugar crops", 
      #eaear_Sunflower = "Sunflower",
      #eaear_Wheat = "Wheat",
      
      eaear_Coconut = "Coconut", 
      eaear_Oilpalm = "Oil palm",
      
      eaear_Tobacco = "Tobacco", 
      eaear_Citrus = "Citrus",
      
      eaear_Banana = "Banana", 
      eaear_Cocoa_Coffee = "Cocoa or Coffee",
      #eaear_Cocoa  = "Cocoa", 
      #eaear_Coffee = "Coffee", 
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
          legend.title = element_blank()) } %>%  add_brackets(crop_groups, face = "italic")


#### limited crops #### 
rfs_1lag3lead_logtrend_lmtrub_clt10 <- list(all = list(), 
                                            America = list(), 
                                            Africa = list(), 
                                            Asia = list())

# prepare data set in advance, to repeat in all regressions
CNT <- "America"

for(CNT in c("all", "America", "Africa", "Asia")){#, , "all",  "Africa", "Asia" "all", "America", "Africa", 
  
  est_parameters <- list(outcome_variable  = "loss_commodity",
                         continent = CNT,
                         start_year = 2009, 
                         end_year = 2019, 
                         most_correlated_only = FALSE,
                         annual_rfs_controls = FALSE,
                         sjpos = FALSE,
                         lags = 1,
                         leads = 3,
                         fya = 0, 
                         pya = 0,
                         s_trend = FALSE, 
                         s_trend_loga = TRUE,
                         fe = "grid_id + country_year",
                         cluster_var1 = "grid_id_10",
                         distribution = "quasipoisson")
  d <- main_data
  
  # Keep only in data the useful variables 
  d <- dplyr::select(d, all_of(unique(c("grid_id", "year", "lat", "lon", "continent_name", "country_name", "country_year",
                                        # "fc_2000", "fc_2009", "remaining_fc", # "accu_defo_since2k",
                                        "grid_id_5", "grid_id_10", "grid_id_20", "grid_id_5_year", "grid_id_10_year", "grid_id_20_year",
                                        # "tmf_agri", "tmf_flood", "tmf_plantation", "tmf_deforestation",
                                        "loss_commodity",
                                        "pasture_share_2000",
                                        eaear_mapmat[,"Crops"] ))))
  
  # Merge only the prices needed, not the whole price dataframe
  d <- left_join(d, prices[,c("year", unique(c(all_rfs_treatments)))], by = c("year"))
  
  # - are suitable to crop j 
  if(est_parameters[["sjpos"]]){
    d <- dplyr::filter(d, !!as.symbol(exposure_rfs) > 0)
  }
  
  # - are in study period
  # if(start_year != 2011 | end_year != 2019){
  d <- dplyr::filter(d, year >= est_parameters[["start_year"]])
  d <- dplyr::filter(d, year <= est_parameters[["end_year"]])
  # }
  
  # - are in study area
  if(est_parameters[["continent"]] != "all"){
    d <- dplyr::filter(d, continent_name == est_parameters[["continent"]])
  }
  
  # Have tmf deforestation to agriculture at least once
  temp_est <- feglm(fml = as.formula(paste0("loss_commodity", " ~ 1 | ", est_parameters[["fe"]])),
                    data = d,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    pre_d_clean_agri <- d[unlist(temp_est$obs_selection),]
  }  else { 
    pre_d_clean_agri <- d
  }
  
  # # Have some tmf deforestation to plantation at least once
  # temp_est <- feglm(fml = as.formula(paste0("loss_commodity", " ~ 1 | ", est_parameters[["fe"]])),
  #                   data = d,
  #                   family = "poisson")
  # # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  # if(length(temp_est$obs_selection)>0){
  #   pre_d_clean_plantation <- d[unlist(temp_est$obs_selection),]
  # }  else { 
  #   pre_d_clean_plantation <- d
  # }
  
  rm(d)
  
  # Regressions of tmf_agri
  elm <- 1
  for(rfs_exp in limited_crops1){#agri_crops eaear_mapmat[,"Crops"]
    rfs_1lag3lead_logtrend_lmtrub_clt10[[CNT]][[elm]] <- make_main_reg(pre_process = TRUE,
                                                                       pre_processed_data = pre_d_clean_agri, # NOTICE THIS
                                                                       
                                                                       outcome_variable = "loss_commodity", # AND THIS
                                                                       continent = est_parameters[["continent"]],
                                                                       start_year = est_parameters[["start_year"]],
                                                                       end_year = est_parameters[["end_year"]],
                                                                       
                                                                       aggr_dyn = TRUE,
                                                                       exposure_rfs = rfs_exp,
                                                                       
                                                                       original_rfs_treatments = c("statute_conv"),
                                                                       rfs_lead = est_parameters[["leads"]], 
                                                                       rfs_lag = est_parameters[["lags"]],
                                                                       rfs_fya = est_parameters[["fya"]],  #ya,
                                                                       rfs_pya = est_parameters[["pya"]], #ya,
                                                                       
                                                                       most_correlated_only = est_parameters[["most_correlated_only"]],
                                                                       annual_rfs_controls = est_parameters[["annual_rfs_controls"]],
                                                                       
                                                                       s_trend = est_parameters[["s_trend"]],
                                                                       s_trend_loga = est_parameters[["s_trend_loga"]],
                                                                       fe = est_parameters[["fe"]], 
                                                                       cluster_var1 = est_parameters[["cluster_var1"]], 
                                                                       
                                                                       rfs_rando = ""
    )
    
    names(rfs_1lag3lead_logtrend_lmtrub_clt10[[CNT]])[elm] <- rfs_exp
    elm <- elm + 1
  }
  
}

rm(df)
dyn_df_list <- list()

for(CNT in c("all", "America", "Africa", "Asia")){ # 1 and 2 correspond respectively to fya and pya, or aggrleads and aggrlags
  dyn_df <- lapply(rfs_1lag3lead_logtrend_lmtrub_clt10[[CNT]], FUN = function(x){as.data.frame(x)[1,] }) %>% bind_rows()
  dyn_df$term <- gsub(pattern = "_X_.*$", x = row.names(dyn_df), replacement = "") # replace everything after and including _X_ with nothing
  dyn_df$model <- CNT # dyn_des[dyn]
  dyn_df_list[[CNT]] <- dyn_df
}
df <- bind_rows(dyn_df_list)
df[df$model == "all", "model"] <- "Global"
#df$term <- rep(eaear_mapmat[,"Crops"], length(dyn_df_list))
df$significant01 <- ""
df[df[,"Pr(>|t|)"] < 0.1, "significant01"] <- "p-value < .1"

names(df)[names(df)=="Estimate"] <- "estimate"
names(df)[names(df)=="Std. Error"] <- "std.error"
row.names(df) <- NULL


{dwplot(df,
        dot_args = list(size = 2.5,aes(shape = model, alpha = significant01)),#, alpha = significant01
        whisker_args = list(size = 1, aes(alpha = significant01)),#alpha = significant01
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(#eaear_Barley = "Barley",
      #eaear_Groundnut = "Groundnut", 
      eaear_Maizegrain = "Maize",
      
      eaear_Fodder = "Pasture",
      
      eaear_Cereals = "Cereals",
      eaear_Rice = "Rice",
      
      eaear_Soy_compo = "Soy", 
      
      eaear_Cotton = "Cotton",
      eaear_Oilfeed_crops = "Oil crops",
      
      #eaear_Rapeseed = "Rapeseed", 
      # eaear_Sorghum2 = "Sorghum", 
      eaear_Biomass = "Biomass crops",
      eaear_Sugarcane = "Sugar crops", 
      #eaear_Sunflower = "Sunflower",
      #eaear_Wheat = "Wheat",
      
      eaear_Coconut = "Coconut", 
      eaear_Oilpalm = "Oil palm",
      
      eaear_Tobacco = "Tobacco", 
      eaear_Citrus = "Citrus",
      
      eaear_Banana = "Banana", 
      eaear_Cocoa_Coffee = "Cocoa or Coffee",
      #eaear_Cocoa  = "Cocoa", 
      #eaear_Coffee = "Coffee", 
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



