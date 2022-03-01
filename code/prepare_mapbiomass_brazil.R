##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

neededPackages <- c("data.table", "plyr", "tidyr", "dplyr",  "Hmisc", "sjmisc", "stringr",
                    "here", "readstata13", "foreign", "readxl", "writexl",
                    "raster", "rgdal", "sp", "spdep", "sf","gfcanalysis",  "nngeo", # "osrm", "osrmr",
                    "lubridate","exactextractr",
                    "doParallel", "foreach", "snow", 
                    "knitr", "kableExtra",
                    "DataCombine", 
                    "fixest", 
                    "boot", "fwildclusterboot", "sandwich",
                    "ggplot2", "leaflet", "tmap", "dotwhisker")

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
# troublePackages <- c("")
# # Attempt to load packages from user's default libraries. 
# lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

# 3. If the troubling packages could not be loaded ("there is no package called ‘’") 
#   you should try to install them, preferably in their versions stated in the renv.lock file. 
#   see in particular https://rstudio.github.io/renv/articles/renv.html 



### NEW FOLDERS USED IN THIS SCRIPT 
dir.create(here("temp_data", "processed_mapbiomass", "brazil_aoi"), recursive = TRUE)
dir.create(here("temp_data", "processed_lossdrivers", "brazil_aoi"), recursive = TRUE)

### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")


# This is to get the bounding box of Brazil as a region for Mapbiomass data aggregation in GEE 
countries <- st_read(here("input_data", "Global_LSIB_Polygons_Detailed"))
brazil <- countries[countries$COUNTRY_NA == "Brazil", "geometry"]
plot(brazil)
brazil_bb <- st_bbox(brazil) %>% st_as_sfc
plot(brazil_bb, add = T)
brazil_bb

# However, I manually changed the -33.75099 latitude to -30 to match the tropical aoi. 
# Moreover, with the aggregation, the precise extent of the 5km output from GEE is now slightly different. 
#  For these reasons, it is more appropriate to use the extent of the output from GEE in this script.  



# All aggregations are performed to the curtis drivers data, because this is the lowest resolution we will need
drivers <- raster(here("input_data", "curtis", "Goode_FinalClassification_19_05pcnt_prj", "Goode_FinalClassification_19_05pcnt_prj.tif"))
# first crop it to brazil aoi 
drivers <- crop(drivers, brazil_aoi)
# plot(drivers)


#### AGGREGATE AND ALIGNE PASTURES ####
## Brick layers 

# import annual layers of forest loss (in hectares) as computed in GEE (and downloaded from Google Drive to input_data/)
rasterlist <- list.files(path = here("input_data", "MAPBIOMASS"), 
                         pattern = "^MapBiomass60_5km_pasture", 
                         full.names = TRUE) %>% as.list()
pastures <- brick(rasterlist)

brazil_aoi <- extent(pastures)

# values are in hectares of pastures in 5km x 5km grid cells - and THERE ARE NAs 
pastures2001 <- pastures$MapBiomass60_5km_pasture2001
# plot(pastures2001)
valpast01 <- values(pastures2001)
anyNA(valpast01)

writeRaster(pastures[[1:19]], here("temp_data", "processed_mapbiomass", "brazil_aoi", "pastures0119.tif"), 
            overwrite = TRUE)

pastures <- brick(here("temp_data", "processed_mapbiomass", "brazil_aoi", "sugarcane0119.tif"))

# define output file name
aggr_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "aggr_sugarcane0119.tif")

# aggregate it from the ~5km cells to ~10km
raster::aggregate(pastures, fact = c(res(drivers)[1]/res(pastures)[1], res(drivers)[2]/res(pastures)[2]),
                  expand = FALSE,
                  fun = sum,
                  na.rm = FALSE, # NA values are on the eastern band only. Thus they can contaminate aggregation safely. 
                  filename = aggr_output_name,
                  overwrite = TRUE)

# align to DRIVERS exactly 
aggregated <- brick(aggr_output_name)

pastures_resampled_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "resampled_sugarcane0119.tif")

resample(x = aggregated, 
         y = drivers, 
         method = "bilinear", # bilinear or ngb changes nothing 
         filename = pastures_resampled_output_name, 
         overwrite = TRUE)


#### AGGREGATE AND ALIGNE SUGARCANE ####
## Brick layers 

# import annual layers of forest loss (in hectares) as computed in GEE (and downloaded from Google Drive to input_data/)
rasterlist <- list.files(path = here("input_data", "MAPBIOMASS"), 
                         pattern = "^MapBiomass60_5km_sugarcane", 
                         full.names = TRUE) %>% as.list()
sugarcane <- brick(rasterlist)

brazil_aoi <- extent(sugarcane)

# values are in hectares of sugarcane in 5km x 5km grid cells - and THERE ARE NAs 
sugarcane2001 <- sugarcane$MapBiomass60_5km_sugarcane2001
# plot(sugarcane2001)
valpast01 <- values(sugarcane2001)
anyNA(valpast01)

writeRaster(sugarcane[[1:19]], here("temp_data", "processed_mapbiomass", "brazil_aoi", "sugarcane0119.tif"), 
            overwrite = TRUE)

sugarcane <- brick(here("temp_data", "processed_mapbiomass", "brazil_aoi", "sugarcane0119.tif"))

# define output file name
aggr_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "aggr_sugarcane0119.tif")

# aggregate it from the ~5km cells to ~10km
raster::aggregate(sugarcane, fact = c(res(drivers)[1]/res(sugarcane)[1], res(drivers)[2]/res(sugarcane)[2]),
                  expand = FALSE,
                  fun = sum,
                  na.rm = FALSE, # NA values are on the eastern band only. Thus they can contaminate aggregation safely. 
                  filename = aggr_output_name,
                  overwrite = TRUE)

# align to DRIVERS exactly 
aggregated <- brick(aggr_output_name)

sugarcane_resampled_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "resampled_sugarcane0119.tif")

resample(x = aggregated, 
         y = drivers, 
         method = "bilinear", # bilinear or ngb changes nothing 
         filename = sugarcane_resampled_output_name, 
         overwrite = TRUE)



#### AGGREGATE AND ALIGNE GAEZ ####
gaez_dir <- here("temp_data", "GAEZ", "v4", "AEAY_out_density",  "Rain-fed")
gaez <- brick(here(gaez_dir, "high_input_all.tif"))
gaez <- crop(gaez, brazil_aoi)

# gaez$high_input_all.1 %>% values() %>% summary()

# resample directly (without aggregating first) from the ~9km cells to ~10km (aggregate does not work bc resolutions are to close)

# define output file name
gaez_resampled_output_name <- here(gaez_dir, "brazil_resampleddrivers_high_input.tif")

# aligne
resample(x = gaez, 
         y = drivers, 
         method = "ngb", # not a big difference between bilinear and ngb, but the latter is still a bit closer to initial gaez, and the former yields negative values.  
         filename = gaez_resampled_output_name, 
         overwrite = TRUE)

# # takes the first layer
# resbil <- raster(resampled_output_name)
# resbil %>% values%>% summary()
# rm(resbil)
# 
# # aligne
# resample(x = gaez, 
#          y = drivers, 
#          method = "ngb", # bilinear or ngb changes nothing 
#          filename = resampled_output_name, 
#          overwrite = TRUE)
# 
# resngb <- raster(resampled_output_name)
# 
# resngb %>% values%>% summary()
# gaez$high_input_all.1 %>% values %>% summary()


#### COMMODITY DRIVEN DEFORESTATION DATA ####
# it is already aggregated at drivers resolution. Just crop it to the brazil aoi, and aline it to the previous rasters
commodity_lossdrivers <- brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", "agri_lossdrivers_commodity.tif"))

commodity_lossdrivers <- crop(commodity_lossdrivers, brazil_aoi)

lossdrivers_resampled_output_name <- here("temp_data", "processed_lossdrivers", "brazil_aoi", "agri_lossdrivers_commodity.tif")
resample(x = commodity_lossdrivers, 
         y = drivers, 
         method = "ngb", # no difference between bilinear and ngb
         filename = lossdrivers_resampled_output_name, 
         overwrite = TRUE)


#### PREPARE ALWAYS ZERO DEFORESTATION MASK #### 
commodity_lossdrivers <- brick(lossdrivers_resampled_output_name)

always_zero <- function(y){if_else(condition = (sum(y)==0), true = 0, false = 1)}
mask_path <- here("temp_data", "processed_lossdrivers", "brazil_aoi", "always_zero_mask_lossdrivers_commodity.tif")

overlay(x = commodity_lossdrivers,
        fun = always_zero,
        filename = mask_path,
        na.rm = TRUE,
        overwrite = TRUE)

#### STACK AND MASK MAPBIOMASS, GAEZ, AND COMMODITY DEFORESTATION #### 

# First, it is important to explicitly rename layers that are going to be stacked and then called to reshape the data frame 
# Rename deforestation layers 
names(commodity_lossdrivers) <- paste0("driven_loss_commodity.",seq(2001, 2019, 1)) 

# read and rename pasture
pastures <- brick(pastures_resampled_output_name)
names(pastures) <- paste0("pastures",seq(2001, 2019, 1)) # this is indeed what has been selected in this land use preparation above

# read and rename sugarcane 
sugarcane <- brick(sugarcane_resampled_output_name)
names(sugarcane) <- paste0("sugarcane",seq(2001, 2019, 1)) # this is indeed what has been selected in this land use preparation above

# read and rename gaez 
gaez <- brick(gaez_resampled_output_name)
gaez_crops <- list.files(path = here(gaez_dir, "High-input"), 
                         pattern = "", 
                         full.names = FALSE)
gaez_crops <- gsub(pattern = ".tif", replacement = "", x = gaez_crops)
names(gaez) <- gaez_crops

# add the 2000 forest cover
# it is not already masked
fc2k <- raster(here("temp_data", "processed_lossdrivers", "brazil_aoi", "fc_2000.tif"))
names(fc2k) <- "fc_2000"

# stack 
brazil_stack <- stack(commodity_lossdrivers, 
                      pastures, 
                      sugarcane,
                      gaez, 
                      fc2k)
names(brazil_stack)

### MASK TO REMOVE ALWAYS ZERO PIXELS AND LIGHTEN THE DATA FRAMES ### 
mask <- raster(mask_path)

brazil_stack <- mask(x = brazil_stack, 
                         mask = mask,
                         maskvalue = 0, # necessary here, because the there is no NA in the mask, only 0 and 1 (see the prepare_loss_drivers.R script)
                         updatevalue = NA)


# (note that masking changes the summary values)


### RASTER TO DATAFRAME ### 

# na.rm = TRUE is key here, as it removes previously masked pixels (NA) and ensures the output is not too large (memory intensive)
# We also set long to false because we reshape with a proper function for more control
wide_df <- raster::as.data.frame(driverloss_gaez, na.rm = TRUE, xy = TRUE, centroids = TRUE, long = FALSE) # ~700s. 

# Rename coordinate variables
names(wide_df)
head(wide_df[,c("x", "y")])
wide_df <- dplyr::rename(wide_df, lon = x, lat = y)


### WIDE TO LONG ### 

# Since we merged datasets in the raster format, we wont need a lonlat format id for each grid cell. 
# So we can simply create an ID that's a sequence. 
wide_df$grid_id <- seq(1, nrow(wide_df), 1) 

# the dot is, by construction of all variable names, only in the names of time varying variables. 
# fixed = TRUE is necessary (otherwise the dot is read as a regexp I guess)
# Note also that it is important that it is structured in a LIST when there are several varying variables in the *long* format
# Because: "Notice that the order of variables in varying is like x.1,y.1,x.2,y.2."
varying_vars <- list(names(driverloss_commodity),
                     names(driverloss_shifting),
                     names(driverloss_forestry),
                     names(driverloss_fire))
#varying_vars <- names(driverloss_gaez)[grep(".", names(driverloss_gaez), fixed = TRUE)]

# reshape to long.
long_df <- stats::reshape(wide_df,
                          varying = varying_vars,
                          v.names = c("driven_loss_commodity", "driven_loss_shifting", "driven_loss_forestry", "driven_loss_fire"),
                          sep = ".",
                          timevar = "year",
                          idvar = "grid_id", # don't put "lon" and "lat" in there, otherwise memory issue (see https://r.789695.n4.nabble.com/reshape-makes-R-run-out-of-memory-PR-14121-td955889.html)
                          ids = "grid_id", # lonlat is our cross-sectional identifier.
                          direction = "long",
                          new.row.names = NULL)#seq(from = 1, to = nrow(ibs_msk_df)*length(years), by = 1)
rm(wide_df)
names(long_df)
# replace the indices from the raster::as.data.frame with actual years.

years <- seq(2001, 2019, 1) # notice here again that it is not the same years as for phtfloss
long_df <- mutate(long_df, year = years[year])

long_df <- dplyr::arrange(long_df, grid_id, year)

# d <- long_df[long_df$driven_loss != long_df$driven_loss_commodity,]
# d <- mutate(d, diff = driven_loss - driven_loss_commodity)
# summary(d$diff)
# d[d$diff>0 , c("driven_loss", "driven_loss_commodity")]

