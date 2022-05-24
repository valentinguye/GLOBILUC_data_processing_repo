##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

neededPackages <- c("data.table", "plyr", "tidyr", "dplyr",  "Hmisc", "sjmisc", "stringr",
                    "here", "readstata13", "foreign", "readxl", "writexl",
                    "raster", "rgdal", "sp", "spdep", "sf","stars", "gfcanalysis",  "nngeo", # "osrm", "osrmr",
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
dir.create(here("temp_data", "processed_fc2000", "brazil_aoi"), recursive = TRUE)
dir.create(here("temp_data", "processed_pasture2000", "brazil_aoi"), recursive = TRUE)
dir.create(here("temp_data", "merged_datasets", "brazil_aoi"), recursive = TRUE)

### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")

### GLOBAL CRS USED throughout the study ### 
mercator_world_crs <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "


# This is to get the bounding box of Brazil as a region for Mapbiomass data aggregation in GEE 
# countries <- st_read(here("input_data", "Global_LSIB_Polygons_Detailed"))
# brazil <- countries[countries$COUNTRY_NA == "Brazil", "geometry"]
# plot(brazil)
# brazil_bb <- st_bbox(brazil) %>% st_as_sfc
# plot(brazil_bb, add = T)
# brazil_bb
# rm(countries, brazil)
# However, I manually changed the -33.75099 latitude to -30 to match the tropical aoi. 
# Moreover, with the aggregation, the precise extent of the 5km output from GEE is now slightly different. 
#  For these reasons, it is more appropriate to use the extent of the output from GEE in this script.  
mapbio <- raster( here("input_data", "MAPBIOMASS", "MapBiomass60_3km_unidir_pasture.tif"))
brazil_aoi <- extent(mapbio)
rm(mapbio)


### GAEZ OBJECTS
# in this script, GAEZ is the target raster of all aggregations / resamplings
gaez_dir <- here("temp_data", "GAEZ", "v4", "AEAY_out_density",  "Rain-fed")
gaez_crops <- list.files(path = here(gaez_dir, "High-input"), 
                         pattern = "", 
                         full.names = FALSE)
gaez_crops <- gsub(pattern = ".tif", replacement = "", x = gaez_crops)

## THIS IS GAEZ IN FULL TROPICAL AOI 
# a priori no issue if it's gaez in global aoi, i.e. not croped in prepare_gaez.R, it only needs to be a larger aoi than continental ones given above. 
# besides, note that we brick the file that was already saved as a single brick of raster layers. 
# Otherwise, calling brick on multiple layers takes some time, and calling stack on multiple layers implies that the object is kept in R memory, and it's ~.06Gb
gaez <- brick(here(gaez_dir, "high_input_all.tif"))

# first crop it to brazil aoi 
gaez_brazil <- crop(gaez, brazil_aoi)
rm(gaez)


### SOME PATHS THAT ARE CALLED THROUGHOUT THE SCRIPT 
pasture_extent_resampled_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "pasture_extent_resampledgaez_0120.tif")
pasture_unidir_resampled_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "pasture_unidir_resampledgaez_0120.tif")

sugarcane_extent_resampled_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "sugarcane_extent_resampledgaez_0120.tif")
sugarcane_unidir_resampled_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "sugarcane_unidir_resampledgaez_0120.tif")

rice_extent_resampled_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "rice_extent_resampledgaez_0120.tif")
rice_unidir_resampled_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "rice_unidir_resampledgaez_0120.tif")

soy_extent_resampled_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "soy_extent_resampledgaez_0120.tif")
soy_unidir_resampled_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "soy_unidir_resampledgaez_0120.tif")

fc2k_resampled_output_name <- here("temp_data", "processed_fc2000", "brazil_aoi", "fc_resampledgaez_2000.tif")
pst2k_resampled_output_name <- here("temp_data", "processed_pasture2000", "brazil_aoi", "pasture_resampledgaez_2000.tif")

commodity_resampled_output_name <- here("temp_data", "processed_lossdrivers", "brazil_aoi", "loss_commo_resampledgaez_0119.tif")

mask_path <- here("temp_data", "processed_lossdrivers", "brazil_aoi", "always_zero_mask_loss_commo_resampledgaez_0119.tif")

#### AGGREGATE AND ALIGNE PASTURE EXTENT ####
## Brick layers 

# import annual layers of forest loss (in hectares) as computed in GEE (and downloaded from Google Drive to input_data/)
# rasterlist <- list.files(path = here("input_data", "MAPBIOMASS"), 
#                          pattern = "^MapBiomass60_5km_pasture", 
#                          full.names = TRUE) %>% as.list()
# writeRaster(pasture[[1:19]], here("temp_data", "processed_mapbiomass", "brazil_aoi", "pasture0119.tif"), 
#             overwrite = TRUE)
# pasture <- brick(here("temp_data", "processed_mapbiomass", "brazil_aoi", "pasture0119.tif"))

pasture <- brick(here("input_data", "MAPBIOMASS", "MapBiomass60_3km_extent_pasture.tif"))

# values are in hectares of pasture in 5km x 5km grid cells - and THERE ARE NAs 
# pasture2001 <- pasture$MapBiomass60_5km_pasture2001
# plot(pasture2001)
# valpast01 <- values(pasture2001)
# anyNA(valpast01)

# define output file name
aggr_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "pasture_extent_9km_0120.tif")

# aggregate it from the ~5km cells to ~10km
raster::aggregate(pasture, fact = 3,
                  expand = FALSE,
                  fun = sum,
                  na.rm = FALSE, # NA values are on the eastern band only. Thus they can contaminate aggregation safely. 
                  filename = aggr_output_name,
                  overwrite = TRUE)

# align to GAEZ exactly 
aggregated <- brick(aggr_output_name)

resample(x = aggregated, 
         y = gaez_brazil, 
         method = "bilinear", # bilinear or ngb changes nothing 
         filename = pasture_extent_resampled_output_name, 
         overwrite = TRUE)

#### AGGREGATE AND ALIGNE PASTURE EXPANSION ####
pasture <- brick(here("input_data", "MAPBIOMASS", "MapBiomass60_3km_unidir_pasture.tif"))

# define output file name
aggr_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "pasture_unidir_9km_0120.tif")

# aggregate it from the ~5km cells to ~10km
raster::aggregate(pasture, fact = 3,
                  expand = FALSE,
                  fun = sum,
                  na.rm = FALSE, # NA values are on the eastern band only. Thus they can contaminate aggregation safely. 
                  filename = aggr_output_name,
                  overwrite = TRUE)

# align to GAEZ exactly 
aggregated <- brick(aggr_output_name)

resample(x = aggregated, 
         y = gaez_brazil, 
         method = "bilinear", # bilinear or ngb changes nothing 
         filename = pasture_unidir_resampled_output_name, 
         overwrite = TRUE)

#### AGGREGATE AND ALIGNE SUGARCANE EXTENT ####
## Bricked layers 
sugarcane <- brick(here("input_data", "MAPBIOMASS", "MapBiomass60_3km_extent_sugarcane.tif"))

# define output file name
aggr_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "sugarcane_extent_9km_0120.tif")

# aggregate it from the ~5km cells to ~10km
raster::aggregate(sugarcane, fact = 3,
                  expand = FALSE,
                  fun = sum,
                  na.rm = FALSE, # NA values are on the eastern band only. Thus they can contaminate aggregation safely. 
                  filename = aggr_output_name,
                  overwrite = TRUE)

# align to GAEZ exactly 
aggregated <- brick(aggr_output_name)

resample(x = aggregated, 
         y = gaez_brazil, 
         method = "bilinear", # bilinear or ngb changes nothing 
         filename = sugarcane_extent_resampled_output_name, 
         overwrite = TRUE)

#### AGGREGATE AND ALIGNE SUGARCANE EXPANSION ####
sugarcane <- brick(here("input_data", "MAPBIOMASS", "MapBiomass60_3km_unidir_sugarcane.tif"))

# define output file name
aggr_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "sugarcane_unidir_9km_0120.tif")

# aggregate it from the ~5km cells to ~10km
raster::aggregate(sugarcane, fact = 3,
                  expand = FALSE,
                  fun = sum,
                  na.rm = FALSE, # NA values are on the eastern band only. Thus they can contaminate aggregation safely. 
                  filename = aggr_output_name,
                  overwrite = TRUE)

# align to GAEZ exactly 
aggregated <- brick(aggr_output_name)

resample(x = aggregated, 
         y = gaez_brazil, 
         method = "bilinear", # bilinear or ngb changes nothing 
         filename = sugarcane_unidir_resampled_output_name, 
         overwrite = TRUE)


#### AGGREGATE AND ALIGNE RICE EXTENT ####
## Bricked layers 
rice <- brick(here("input_data", "MAPBIOMASS", "MapBiomass60_3km_extent_rice.tif"))

# define output file name
aggr_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "rice_extent_9km_0120.tif")

# aggregate it from the ~5km cells to ~10km
raster::aggregate(rice, fact = 3,
                  expand = FALSE,
                  fun = sum,
                  na.rm = FALSE, # NA values are on the eastern band only. Thus they can contaminate aggregation safely. 
                  filename = aggr_output_name,
                  overwrite = TRUE)

# align to GAEZ exactly 
aggregated <- brick(aggr_output_name)

resample(x = aggregated, 
         y = gaez_brazil, 
         method = "bilinear", # bilinear or ngb changes nothing 
         filename = rice_extent_resampled_output_name, 
         overwrite = TRUE)

#### AGGREGATE AND ALIGNE RICE EXPANSION ####
rice <- brick(here("input_data", "MAPBIOMASS", "MapBiomass60_3km_unidir_rice.tif"))

# define output file name
aggr_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "rice_unidir_9km_0120.tif")

# aggregate it from the ~5km cells to ~10km
raster::aggregate(rice, fact = 3,
                  expand = FALSE,
                  fun = sum,
                  na.rm = FALSE, # NA values are on the eastern band only. Thus they can contaminate aggregation safely. 
                  filename = aggr_output_name,
                  overwrite = TRUE)

# align to GAEZ exactly 
aggregated <- brick(aggr_output_name)

resample(x = aggregated, 
         y = gaez_brazil, 
         method = "bilinear", # bilinear or ngb changes nothing 
         filename = rice_unidir_resampled_output_name, 
         overwrite = TRUE)


#### AGGREGATE AND ALIGNE SOY EXTENT ####
## Bricked layers 
soy <- brick(here("input_data", "MAPBIOMASS", "MapBiomass60_3km_extent_soy.tif"))

# define output file name
aggr_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "soy_extent_9km_0120.tif")

# aggregate it from the ~5km cells to ~10km
raster::aggregate(soy, fact = 3,
                  expand = FALSE,
                  fun = sum,
                  na.rm = FALSE, # NA values are on the eastern band only. Thus they can contaminate aggregation safely. 
                  filename = aggr_output_name,
                  overwrite = TRUE)

# align to GAEZ exactly 
aggregated <- brick(aggr_output_name)

resample(x = aggregated, 
         y = gaez_brazil, 
         method = "bilinear", # bilinear or ngb changes nothing 
         filename = soy_extent_resampled_output_name, 
         overwrite = TRUE)

#### AGGREGATE AND ALIGNE SOY EXPANSION ####
soy <- brick(here("input_data", "MAPBIOMASS", "MapBiomass60_3km_unidir_soy.tif"))

# define output file name
aggr_output_name <- here("temp_data", "processed_mapbiomass", "brazil_aoi", "soy_unidir_9km_0120.tif")

# aggregate it from the ~5km cells to ~10km
raster::aggregate(soy, fact = 3,
                  expand = FALSE,
                  fun = sum,
                  na.rm = FALSE, # NA values are on the eastern band only. Thus they can contaminate aggregation safely. 
                  filename = aggr_output_name,
                  overwrite = TRUE)

# align to GAEZ exactly 
aggregated <- brick(aggr_output_name)

resample(x = aggregated, 
         y = gaez_brazil, 
         method = "bilinear", # bilinear or ngb changes nothing 
         filename = soy_unidir_resampled_output_name, 
         overwrite = TRUE)



#### 2000 FOREST COVER #### 
fc2k <- raster(here("input_data", "fc_2000_3km_10th.tif"))

fc2k <- crop(fc2k, brazil_aoi)

fc2k_aggr_output_name <- here("temp_data", "processed_fc2000", "brazil_aoi", "fc_9km_2000.tif")

aggregate(fc2k, fact = 3,
          expand = FALSE,
          fun = sum,
          na.rm = TRUE,  # no NA values a priori mais bon 
          filename = fc2k_aggr_output_name,
          overwrite = TRUE)

# align to GAEZ exactly 
aggregated_fc2000 <- raster(fc2k_aggr_output_name)

resample(x = aggregated_fc2000, 
         y = gaez_brazil, 
         method = "ngb", # ngb because bilinear yields negatie values
         filename = fc2k_resampled_output_name, 
         overwrite = TRUE)

#### 2000 PASTURE SHARE ####
pst2k <- raster(here("input_data", "CroplandPastureArea2000_Geotiff", "Pasture2000_5m.tif"))
pst2k <- crop(pst2k, brazil_aoi)
# resample directly (without aggregating first) from the ~9km cells to ~10km (aggregate does not work bc resolutions are to close)

# transform NAs to 0, such that the more numerous NAs in pasture data do not force losing information when the stack is turned to a data frame. 
pst2k <- reclassify(pst2k, cbind(NA, 0))

resample(x = pst2k, 
         y = gaez_brazil, 
         method = "ngb", # we use ngb and not bilinear because the output values' summary better fits that of the aggregated layer 
         # and the bilinear interpolation arguably smoothes the reprojection more than necessary given that from and to are already very similar.  
         filename = pst2k_resampled_output_name, 
         overwrite = TRUE)


#### DRIVEN DEFORESTATION DATA ####

### COMMODITY ###
# it is already aggregated at gaez_brazil resolution. Just crop it to the brazil aoi, and aline it to the previous rasters
losscommo <- brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", paste0("loss_commo_resampledgaez_0119.tif")))

losscommo <- crop(losscommo, brazil_aoi)

resample(x = losscommo, 
         y = gaez_brazil, 
         method = "ngb", # no difference between bilinear and ngb
         filename = commodity_resampled_output_name, 
         overwrite = TRUE)


#### PREPARE ALWAYS ZERO DEFORESTATION MASK #### 
losscommo <- brick(commodity_resampled_output_name)

always_zero <- function(y){if_else(condition = (sum(y)==0), true = 0, false = 1)}

overlay(x = losscommo,
        fun = always_zero,
        filename = mask_path,
        na.rm = TRUE,
        overwrite = TRUE)



#### STACK AND MASK MAPBIOMASS, GAEZ, AND COMMODITY DEFORESTATION #### 
 # Read layers to be stacked
losscommo <- brick(commodity_resampled_output_name)

pasture_extent <- brick(pasture_extent_resampled_output_name)
sugarcane_extent <- brick(sugarcane_extent_resampled_output_name)
rice_extent <- brick(rice_extent_resampled_output_name)
soy_extent <- brick(soy_extent_resampled_output_name)

pasture_unidir <- brick(pasture_unidir_resampled_output_name)
sugarcane_unidir <- brick(sugarcane_unidir_resampled_output_name)
rice_unidir <- brick(rice_unidir_resampled_output_name)
soy_unidir <- brick(soy_unidir_resampled_output_name)

fc2k <- raster(fc2k_resampled_output_name)
pst2k <- raster(pst2k_resampled_output_name)

# It is important to explicitly rename layers that are going to be stacked and then called to reshape the data frame 
# for time varying variables, the dot is important. 
names(losscommo) <- paste0("loss_commodity.",seq(2001, 2019, 1)) 

names(pasture_extent) <- paste0("extent_pasture.",seq(2001, 2020, 1)) # this is indeed what has been selected in this land use preparation above
names(sugarcane_extent) <- paste0("extent_sugarcane.",seq(2001, 2020, 1)) # this is indeed what has been selected in this land use preparation above
names(rice_extent) <- paste0("extent_rice.",seq(2001, 2020, 1)) # this is indeed what has been selected in this land use preparation above
names(soy_extent) <- paste0("extent_soy.",seq(2001, 2020, 1)) # this is indeed what has been selected in this land use preparation above

names(pasture_unidir) <- paste0("unidir_pasture.",seq(2001, 2020, 1)) # this is indeed what has been selected in this land use preparation above
names(sugarcane_unidir) <- paste0("unidir_sugarcane.",seq(2001, 2020, 1)) # this is indeed what has been selected in this land use preparation above
names(rice_unidir) <- paste0("unidir_rice.",seq(2001, 2020, 1)) # this is indeed what has been selected in this land use preparation above
names(soy_unidir) <- paste0("unidir_soy.",seq(2001, 2020, 1)) # this is indeed what has been selected in this land use preparation above

names(gaez_brazil) <- gaez_crops

names(fc2k) <- "fc_2000"
names(pst2k) <- "pasture_share_2000"

## STACK 
brazil_stack <- stack(losscommo,
                      pasture_extent, sugarcane_extent, rice_extent, soy_extent, 
                      pasture_unidir, sugarcane_unidir, rice_unidir, soy_unidir, 
                      gaez_brazil, fc2k, pst2k)
names(brazil_stack)

# save the stack at this point, as a clean export 
# writeRaster(brazil_stack,
#             filename = here("temp_data", "merged_datasets", "brazil_aoi", "drivenloss_mapbiomass_gaez"),
#             overwrite = TRUE, 
#             format = "raster")

### MASK TO REMOVE ALWAYS ZERO PIXELS AND LIGHTEN THE DATA FRAMES ### 
mask <- raster(mask_path)

# ~500 seconds
brazil_stack <- mask(x = brazil_stack, 
                         mask = mask,
                         maskvalue = 0, # necessary here, because the there is no NA in the mask, only 0 and 1 (see the prepare_loss_drivers.R script)
                         updatevalue = NA)


# (note that masking changes the summary values)


### RASTER TO DATAFRAME ### 

# na.rm = TRUE is key here, as it removes previously masked pixels (NA) and ensures the output is not too large (memory intensive)
# We also set long to false because we reshape with a proper function for more control
wide_df <- raster::as.data.frame(brazil_stack, na.rm = TRUE, xy = TRUE, centroids = TRUE, long = FALSE) # ~700s. 

# Rename coordinate variables
names(wide_df)
head(wide_df[,c("x", "y")])
wide_df <- dplyr::rename(wide_df, lon = x, lat = y)


# SPECIAL TO SOY SCRIPT HERE: ADD 2001 COLUMN FOR THE SOY EXPANSION  
# because we want 2001 extent to be in the final data set, but expansion (unidir) is only avaiable from 2002 and number of years must be 
# equal for all varying vars

# unidir 2001 layers are available, since Mapbiomass data was available from 1985

wide_df$loss_commodity.2020 <- NA

### WIDE TO LONG ### 

# Since we merged datasets in the raster format, we wont need a lonlat format id for each grid cell. 
# So we can simply create an ID that's a sequence. 
wide_df$grid_id <- seq(1, nrow(wide_df), 1) 

# the dot is, by construction of all variable names, only in the names of time varying variables. 
# fixed = TRUE is necessary (otherwise the dot is read as a regexp I guess)
# Note also that it is important that it is structured in a LIST when there are several varying variables in the *long* format
# Because: "Notice that the order of variables in varying is like x.1,y.1,x.2,y.2."
varying_vars <- list(paste0("loss_commodity.", seq(2001, 2020, 1)),
                     paste0("extent_pasture.",seq(2001, 2020, 1)), 
                     paste0("unidir_pasture.",seq(2001, 2020, 1)), 
                     paste0("extent_sugarcane.",seq(2001, 2020, 1)), 
                     paste0("unidir_sugarcane.",seq(2001, 2020, 1)), 
                     paste0("extent_rice.",seq(2001, 2020, 1)), 
                     paste0("unidir_rice.",seq(2001, 2020, 1)), 
                     paste0("extent_soy.",seq(2001, 2020, 1)), 
                     paste0("unidir_soy.",seq(2001, 2020, 1)))
#varying_vars <- names(driverloss_gaez)[grep(".", names(driverloss_gaez), fixed = TRUE)]

# reshape to long.
long_df <- stats::reshape(wide_df,
                          varying = varying_vars,
                          v.names = c("loss_commodity", 
                                      "extent_pasture", "extent_sugarcane", "extent_rice", "extent_soy", 
                                      "unidir_pasture", "unidir_sugarcane", "unidir_rice", "unidir_soy"),
                          sep = ".",
                          timevar = "year",
                          idvar = "grid_id", # don't put "lon" and "lat" in there, otherwise memory issue (see https://r.789695.n4.nabble.com/reshape-makes-R-run-out-of-memory-PR-14121-td955889.html)
                          ids = "grid_id", # lonlat is our cross-sectional identifier.
                          direction = "long",
                          new.row.names = NULL)#seq(from = 1, to = nrow(ibs_msk_df)*length(years), by = 1)
rm(wide_df)
names(long_df)
# replace the indices from the raster::as.data.frame with actual years.

years <- seq(2001, 2020, 1) # notice here again that it is not the same years as for phtfloss
long_df <- mutate(long_df, year = years[year])

long_df <- dplyr::arrange(long_df, grid_id, year)

saveRDS(long_df, here("temp_data", "merged_datasets", "brazil_aoi", "loss_commodity_aeay_long.Rdata"))

rm(long_df, varying_vars, gaez_brazil, mask, fc2k, pst2k)


#### BIGGER CELL VARIABLES #### 
## Prepare base data
path <- here("temp_data", "merged_datasets", "brazil_aoi", "loss_commodity_aeay_long.Rdata")
df <- readRDS(path)

# Remove gaez variables
df <- dplyr::select(df,-all_of(gaez_crops))
# Use cross section only
df_cs <- df[!duplicated(df$grid_id),]
# Spatial
df_cs <- st_as_sf(df_cs, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

rm(df)

## Prepare bigger square grids
grid_base <- brazil_stack[[1]]

# for ~45km grid cells (5 times larger grid cells in both dimensions, hence 25 times larger)
bigger_5 <- aggregate(grid_base, fact = 5, expand = TRUE, fun = sum)
bigger_5_stars <- st_as_stars(bigger_5)
bigger_5_sf <- st_as_sf(bigger_5_stars, as_points = FALSE, merge = FALSE, na.rm = FALSE, long = FALSE)

# Make the grid cell index. It is continent-specific, otherwise it won't be unique across continents. 
bigger_5_sf$grid_id_5 <- seq(from = 1, to = nrow(bigger_5_sf))
bigger_5_sf <- bigger_5_sf[,c("grid_id_5", "geometry")]

# repeat for ~90km grid cells (10 times larger grid cells in both dimensions, hence 10 times larger)
bigger_10 <- aggregate(grid_base, fact = 10, expand = TRUE, fun = sum)
bigger_10_stars <- st_as_stars(bigger_10)
bigger_10_sf <- st_as_sf(bigger_10_stars, as_points = FALSE, merge = FALSE, na.rm = FALSE, long = FALSE)

# Make the grid cell index. It is continent-specific, otherwise it won't be unique across continents. 
bigger_10_sf$grid_id_10 <- seq(from = 1, to = nrow(bigger_10_sf))
bigger_10_sf <- bigger_10_sf[,c("grid_id_10", "geometry")]

# repeat for ~180km grid cells (20 times larger grid cells in both dimensions, hence 20 times larger)
bigger_20 <- aggregate(grid_base, fact = 20, expand = TRUE, fun = sum)
bigger_20_stars <- st_as_stars(bigger_20)
bigger_20_sf <- st_as_sf(bigger_20_stars, as_points = FALSE, merge = FALSE, na.rm = FALSE, long = FALSE)

# Make the grid cell index. It is continent-specific, otherwise it won't be unique across continents. 
bigger_20_sf$grid_id_20 <- seq(from = 1, to = nrow(bigger_20_sf))
bigger_20_sf <- bigger_20_sf[,c("grid_id_20", "geometry")]

# DO NOT TRANSFORM, because it makes the inner_join associate two bigger squares for each smaller one, 
# I don't really know why, but it does not do so with geographic coordinates, and there is no mismatch, and no problem of imprecision due to 
# geographic coordinates being used as planer because aoi is not near the pole. 
# df_cs <- st_transform(df_cs, crs = mercator_world_crs)
# bigger_50km_sf <- st_transform(bigger_50km_sf, crs = mercator_world_crs)
# bigger_100km_sf <- st_transform(bigger_100km_sf, crs = mercator_world_crs)

# join the variables
df_cs <- st_join(x = df_cs,
                 y = bigger_5_sf,
                 join = st_within,
                 prepared = FALSE, # tests shows that this changes nothing, whether shapes are transformed or not
                 left = TRUE)

df_cs <- st_join(x = df_cs,
                 y = bigger_10_sf,
                 join = st_within,
                 prepared = FALSE, # tests shows that this changes nothing, whether shapes are transformed or not
                 left = TRUE)

df_cs <- st_join(x = df_cs,
                 y = bigger_20_sf,
                 join = st_within,
                 prepared = FALSE, # tests shows that this changes nothing, whether shapes are transformed or not
                 left = TRUE)

length(unique(df_cs$grid_id)) == nrow(df_cs)

df_cs <- st_drop_geometry(df_cs)

# Keep only new variable and id
df_cs <- df_cs[,c("grid_id", "grid_id_5", "grid_id_10", "grid_id_20")]

saveRDS(df_cs, here("temp_data", "merged_datasets", "brazil_aoi", "loss_commodity_aeay_cs_biggercells.Rdata"))


#### ADD REMAINING FOREST VARIABLE #### 
# neither country nor continent information is relevant here. 

df <- readRDS(here("temp_data", "merged_datasets", "brazil_aoi", "loss_commodity_aeay_long.Rdata"))

# Remove gaez variables
df <- dplyr::select(df,-all_of(gaez_crops))

year_list <- list()

# for some grid cells, deforestation from different drivers is positive in the same year. 
# This probably comes from resampling/reprojecting operations, and can be observed here only when forest loss was positive in those places
# dplyr::filter(df, driven_loss_commodity>0 & driven_loss_any > driven_loss_commodity)
# df[df$grid_id == 110,]

# just ensure that driven_loss_any is not counting 4 times too much deforestation in those places. 
# ensuring driven_loss_commodity is not missing, as is the case for the whole year 2020
# df[!(is.na(df$driven_loss_commodity)) & df$driven_loss_commodity>0 & df$driven_loss_any > df$driven_loss_commodity, ] <- dplyr::filter(df, 
#                                                                                                   driven_loss_commodity>0 & driven_loss_any > driven_loss_commodity) %>% 
#   mutate(driven_loss_any = driven_loss_commodity)


# in the first year (2001), the past year accumulated deforestation is null. 
year_list[["2001"]] <- df[df$year == 2001, c("grid_id", "year")] 
year_list[["2001"]][,"accu_defo_since2k"] <- 0

# then, each year's deforestation accumulated in the past is the sum of *past years'* deforestation
years <- 2008:2011 # we need it only 2007-2010, because we possibly use only those years as starting years not 2002:max(df$year)
for(y in years){
  sub_ <- df[df$year < y,]
  year_list[[as.character(y)]] <- ddply(sub_, "grid_id", summarise,
                                        accu_defo_since2k = sum(loss_commodity, na.rm = TRUE))
  year_list[[as.character(y)]][,"year"] <- y
}

accu_defo_df <- bind_rows(year_list)

df <- left_join(df, accu_defo_df, by = c("grid_id", "year"))


# summary(df$accu_lucpfp_since2k)
df <- dplyr::mutate(df, 
                    remaining_fc = fc_2000 - accu_defo_since2k)

fc_2009 <- df[df$year == 2009, c("grid_id", "remaining_fc")]
names(fc_2009) <- c("grid_id", "fc_2009")
df <- left_join(df, fc_2009, by = "grid_id")


# df[df$grid_id == 1267,c("grid_id", "year", "lucpfap_pixelcount", "accu_lucpfp_since2k", "remain_pf_pixelcount")] 

# put keep only new variables in remaining
remaining <- df[,c("grid_id", "year", "remaining_fc", "accu_defo_since2k", "fc_2009")] # fc_2000 is added as a raster layer in merge_* scripts

saveRDS(remaining, here("temp_data", "merged_datasets", "brazil_aoi", "loss_commodity_aeay_long_remaining.Rdata"))

rm(year_list, sub_, accu_defo_df, df)


#### GROUP AND STANDARDIZE AEAY VARIABLES #### 
# all groupings in this section are motivated on the GAEZ v4 model documentation, and in particular Table A4-1.3

df <- readRDS(here("temp_data", "merged_datasets", "brazil_aoi", "loss_commodity_aeay_long.Rdata"))

# Use cross section only
df_cs <- df[!duplicated(df$grid_id),]

rm(df)

# NOW, prices are needed only for SOME crops, those that are grouped but GAEZ yield in dry weights are not comparable. 
prices <- readRDS(here("temp_data", "prepared_international_prices.Rdata"))

### This matrix is used for maping crops from GAEZ with commodities from price data sets
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

# this vector is used to convert potential production in dry weight from GAEZ into weight of traded commodities. 
# Olive, Oil palm, and Sugar crops are already provided in GAEZ in units of traded goods (oil or sugar)
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


# no vlaue missing for any of these commodities in the broad period 1991-2019
# working_prices <- prices[,c("year", mapmat[,"Prices"])]%>%filter(year>1990 & year < 2020)

# get prices for 2000
price_avg <- prices %>% 
  filter(year>=1995 & year <= 2004) %>% 
  #filter(year==2000) %>% 
  summarise(across(.cols = any_of(mapmat[,"Prices"]), 
                   .fns = mean, na.rm = TRUE))


## Aggregate suitability indexes to match price data
# makes sense for SUgar fodder and rice subcrops because yields expressed in comaprable units  
df_cs <- df_cs %>% rowwise() %>% mutate(Biomass = max(c(Miscanthus, Reedcanarygrass, Sorghumbiomass*10, Switchgrass)), # sorghum biomass is expressed in kg/ha, and not 10kg/ha as it is the case for the three other crops
                                        #  don't include Jatropha because it is expressed in seeds and not above ground biomass in GAEZ. 
                                        Fodder = max(c(Alfalfa, Napiergrass)),   
                                        Rice = max(c(Drylandrice, Wetlandrice)),
                                        Sugar = max(c(Sugarbeet, Sugarcane)) # Especially necessary to match the international price of sugar
                                        # We surely wont need these in the AEAY workflow, as we need to interact these variables with a price and there is no price for those. 
                                        # bioenergy_crops = max(c(Jatropha, Miscanthus, Reedcanarygrass, Sorghumbiomass, Switchgrass)),
                                        # cereal_crops = max(c(Buckwheat, Foxtailmillet, Pearlmillet, Rye)),
                                        # pulses_crops = max(c(Chickpea, Cowpea, Drypea, Gram, Phaseolousbean, Pigeonpea)),
                                        # roots_crops = max(c(Cassava, Sweetpotato, Whitepotato, Yam)),
                                        # # oil_crops = max(c(Groundnut, Jatropha, Olive, Rapeseed, Sunflower)), # we have prices for all of them
                                        # vegetables_crops = max(c(Cabbage, Carrot, Onion, Tomato)),
                                        # # fruits_crops = max(c(Banana, Citrus, Coconut)), # For coconut we have price only for the oil while it is not expressed in oil in GAEZ. 
                                        # industrial_crops = max(c(Flax)) # Rubber, directly responsible for deforestation, are not mixed with these.
                                        # # narcotics_crops = max(c(Cocoa, Coffee, Tea, Tobacco)), # We have prices for all of them
) %>% as.data.frame()

# keep only crops of interest 
df_cs <- df_cs[, names(df_cs) %in% c("grid_id", "lon", "lat", mapmat[,"Crops"])]
## TONS OF WHAT? For some crops, the AEAY quantity does not necessarily match the market price unit  
# "For most crops the agro-climatic potential yield is given as kg dry weight per hectare. 
# For alfalfa, miscanthus, switchgrass, reed canary grass, napier grass, pasture legumes and grasses the yield is given in 10 kg dry weight per hectare. 
# For sugar beet and sugarcane (and hence Sugar, the max of them) yields are in kg sugar per hectare, 
# and for oil palm and olives in kg oil per hectare. Cotton yield is given as kg lint per hectare." 
# https://gaez.fao.org/pages/theme-details-theme-3

# BANANA. 
# Just convert from dry matter weight
df_cs <- dplyr::mutate(df_cs, Banana = Banana / conv_fac["Banana"])

# BARLEY. 
# Just convert from dry matter weight
df_cs <- dplyr::mutate(df_cs, Barley = Barley / conv_fac["Barley"])

# COCONUT. 
# Coconut cannot be converted to coconut oil, as it is made of only a by product of coconut; 
# We would not correctly estimate the value of coconut if we counted the whole coconut yield as the copra byproduct;   

# COTTON.
# From FAO documentation and v4 data visualization, the unit of cotton AEAY is kg lint/ha. 
# cotton lint is "raw ginned cotton which is ready for baling" (baling is 'packing' cotton in standard volumes), see there for definitions http://agropedia.iitk.ac.in/content/glossary-useful-terms-related-cotton
# the price in pink sheet is in $/kg and comes from Cotton Outlook A Index, which is expressed for 'raw cotton', see there https://www.cotlook.com/information-2/the-cotlook-indices-an-explanation/
# which I understand as being raw ginned cotton, as it is given for a particular grade, and ginned (baled) cotton is graded, not pure raw cotton (before ginning)
# THUS, the price is expressed in the same unit as the AEAY
# NOTE that the conversion factor of 0.33 used in GAEZ should not be applied, because it converts from seed cotton to cotton lint to match with FAOSTAT data, which is expressed in seed cotton weight. 

# FODDER. 
# Convert fodder crop yield into beef by a feed conversion ratio. 
# Following Galloway et al. 2007 who get a ratio of 20 (feed to meat conversion rate of 0.05) for ruminants (beef and sheep and goats) on non arable land (i.e. for fodder, not feed from crops which is more efficient)
# This is quite similar to Alexander et al. 2016, 25 tons of feed transform into 1 ton of beef meat (in edible weight). 
# Galloway being more specific about FCR of non-arable land feed, we retain this. 
# Lower (more efficient) FCR found in the literature, typically below 10, represent feed, not fodder (often necessary to compare with non-ruminants)
# Thus, every ton of fodder (coming from agro-climatically achievable yields in ton/ha) is scaled to 1/20 ton of beef meat  
# BUT needs to account for a dry matter content conversion to actually grazed feed. 
# We apply an dry matter content of 0.271 as reported for summer green chop alfalfa in https://www.ccof.org/sites/default/files/Feed%20Type%20DMI%20Table%20Final.pdf
# (which has a similar dry matter content as napier grass)
df_cs <- dplyr::mutate(df_cs, Fodder = Fodder * 0.05 / 0.271)

# GROUNDNUT. 
# Oil content is 45-56% in https://link.springer.com/chapter/10.1007%2F978-94-011-0733-4_6 as reporte by https://link.springer.com/article/10.1007/s11746-017-2981-3
# it is 31% in table 32 in https://www.ers.usda.gov/webdocs/publications/41880/33132_ah697_002.pdf
# Because it's a bit unclear what the value is in the literature, and price is available in PS for groudnuts (not oil), we don't convert to oil.
# Nevertheless, we convert from DM weight kernel (in GAEZ) to shelled weight (Pink Sheet), using GAEZ conversion factor
df_cs <- dplyr::mutate(df_cs, Groundnut = Groundnut / conv_fac["Groundnut"])

# MAIZE. 
# Just convert from dry matter weight
df_cs <- dplyr::mutate(df_cs, Maizegrain = Maizegrain / conv_fac["Maizegrain"])

# OAT
df_cs <- dplyr::mutate(df_cs, Oat = Oat / conv_fac["Oat"])

# OLIVE AND OIL PALM. 
# They are expressed in oil already, as in the price data. Hence nothing to do. 

# RAPESEED.
# Seeds contain around 41% oil see Table 4 in Yasar 2018, and https://www.agmrc.org/commodities-products/grains-oilseeds/rapeseed 
df_cs <- dplyr::mutate(df_cs, Rapeseed = Rapeseed * 0.41 / conv_fac["Rapeseed"])

# RUBBER. 
# "In general, latex contains about 30-40 % of rubber particles and 55-65 % of water. However, fresh latex shows 15-45 % of rubber hydrocarbon and about
# 2-4 % of non-rubber ingredients [2]. Latex is usually sold either in the form of dry rubber sheet or concentrated rubber solution."
# https://www.measurement.sk/2014/Kerdtongmee.pdf
# So the RSS3 product the PS price is for, is Ribbed Smoked Sheet (dry rubber sheet), i.e. a raw form of latex, that is said to have a dry rubber content of 30-40%. 
# Dry rubber content from this page (60%) is for concentrated rubber/latex, not RSS http://www.unistarglobal.com/natural_rubber.php 
# So here 1 ton of dry rubber is diluted to produce a higher weight of RSS. 
df_cs <- dplyr::mutate(df_cs, Rubber = Rubber / 0.35)

# SORGHUM. 
# Just convert from dry matter weight
df_cs <- dplyr::mutate(df_cs, Sorghum = Sorghum / conv_fac["Sorghum"])

# SUNFLOWER
# The oil content is set at 42%, https://www.sciencedirect.com/science/article/pii/B9780123849472006747
# it is 40.5% in table 32 in https://www.ers.usda.gov/webdocs/publications/41880/33132_ah697_002.pdf and https://www.agmrc.org/commodities-products/grains-oilseeds/sunflower-profile 
# in Table 4 in Yasar 2018 it's 40-50%, 
# 22–55% oil content (Flagella et al. 2002; Gonzalez-Martin et al. 2013) as reported in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5976617/
df_cs <- dplyr::mutate(df_cs, Sunflower = Sunflower * 0.42 / conv_fac["Sunflower"])

# SOY
# 18% of the seed is extracted as oil and 82% as meal, so we make two distinct yields 
df_cs <- dplyr::mutate(df_cs, Soybean = Soybean / conv_fac["Soybean"])
df_cs <- dplyr::mutate(df_cs, Soybean_oil = Soybean * 0.18 / conv_fac["Soybean"])
df_cs <- dplyr::mutate(df_cs, Soybean_meal = Soybean * 0.82 / conv_fac["Soybean"])
# Compute also the value for mere Soy beans

# RICE 
# Just convert from dry matter weight
df_cs <- dplyr::mutate(df_cs, Rice = Rice / conv_fac["Rice"])

# TOBACCO. 
# prices are for unmanufactured tobacco 
# Just convert from dry matter weight
df_cs <- dplyr::mutate(df_cs, Tobacco = Tobacco / conv_fac["Tobacco"])

# WHEAT. 
# Just convert from dry matter weight
df_cs <- dplyr::mutate(df_cs, Wheat = Wheat / conv_fac["Wheat"])


## CONVERT FROM kg/ha to ton/ha
# All prices have been converted to $/ton in prepare_prices.R but all yields are expressed in kg/ha
df_cs <- dplyr::mutate(df_cs, across(.cols = all_of(mapmat[,"Crops"]),
                                     .fns = ~./1000)) 
# Convert Nappier grass and alfalfa (components of Fodder), and Biomass from now 10 tons to tons. 
df_cs <- dplyr::mutate(df_cs, across(.cols = all_of(c("Fodder", "Biomass")),
                                     .fns = ~.*10)) 

## Interact with average prices to get Expected Agro-Ecological Attainable Revenue (EAEAR)

# Prices have been converted to $/t in prepare_prices.R
# we actually need to convert yields to revenues only for crops that we want to group
# yet, for sthe sake of generality, we convert all crops to their revenues
# for non-grouped crops, it's not a big deal if the conversion to $/ha is not accurate, as comparisons will only be within crops, between grid cells 
for(aeay_i in mapmat[,"Crops"]){
  price_i <- price_avg[mapmat[mapmat[,"Crops"]==aeay_i,"Prices"]]%>%as.numeric()
  eaear_i <- paste0("eaear_", aeay_i)
  df_cs <- dplyr::mutate(df_cs, 
                         !!as.symbol(eaear_i) := !!as.symbol(aeay_i) * price_i)
}

df_cs <- df_cs %>% rowwise() %>% mutate(eaear_Cereals = max(c(eaear_Barley, eaear_Sorghum, eaear_Wheat)), 
                                        eaear_Oilfeed_crops = max(c(eaear_Groundnut, eaear_Rapeseed, eaear_Sunflower)), 
                                        eaear_Cocoa_Coffee = max(c(eaear_Cocoa, eaear_Coffee))) %>% as.data.frame()

# and for Soy commodities: 
df_cs <- dplyr::mutate(df_cs, eaear_Soy_compo =  eaear_Soybean_meal + eaear_Soybean_oil)

# the highest values for Soy_compo represent the value added from processing into oil/meals
summary(df_cs$eaear_Soybean)
summary(df_cs$eaear_Soy_compo)
# Retain processed soy commodities, because this is more traded, and more comparable to the other oil seeds that are expressed in oil too. 
df_cs <- dplyr::select(df_cs, -eaear_Soybean, -eaear_Soybean_meal, -eaear_Soybean_oil)


# Select variables to save: all the eaear variables, standardized or not
var_names <- grep(pattern = "eaear_", names(df_cs), value = TRUE) 
df_cs <- df_cs[,c("grid_id", var_names)]

saveRDS(df_cs, here("temp_data", "merged_datasets", "brazil_aoi",  "loss_commodity_aeay_cs_stdeaear.Rdata"))

rm(df_cs)

#### MERGE ADDITIONAL VARIABLES ####  
# Base dataset (including outcome variable(s))
df_base <- readRDS(here("temp_data", "merged_datasets", "brazil_aoi", "loss_commodity_aeay_long.Rdata"))

# just compute country and continent variables, even if invariant, so they can be called in generic function
df_base$country_name <- "Brazil"
df_base$continent_name <- "America"
# Create country year fixed effect
df_base <- mutate(df_base, country_year = paste0(country_name, "_", year))


## EAEAR
df_stdeaear <- readRDS(here("temp_data", "merged_datasets", "brazil_aoi",  "loss_commodity_aeay_cs_stdeaear.Rdata"))  

final <- left_join(df_base, df_stdeaear, by = "grid_id") # no issue with using grid_id as a key here, bc df_remain was computed just above from the df_base data
rm(df_base, df_stdeaear)

## BIGGER CELLS
df_biggercells <- readRDS(here("temp_data", "merged_datasets", "brazil_aoi", "loss_commodity_aeay_cs_biggercells.Rdata"))

final <- left_join(final, df_biggercells, by = "grid_id")
rm(df_biggercells)

# Create bigger cell-year identifier
final <- mutate(final, grid_id_5_year = paste0(grid_id_5, "_", year))
final <- mutate(final, grid_id_10_year = paste0(grid_id_10, "_", year))
final <- mutate(final, grid_id_20_year = paste0(grid_id_20, "_", year))
# length(unique(final$grid_id_50km_year))==length(unique(final$grid_id_50km))*length(unique(final$year))


## REMAINING
# df_remain <- readRDS(here("temp_data", "merged_datasets", "brazil_aoi", "loss_commodity_aeay_long_remaining.Rdata"))
# 
# final <- left_join(final, df_remain, by = c("grid_id", "year"))  # no issue with using grid_id as a key here, bc df_remain was computed just above from the df_base data
# rm(df_remain)

saveRDS(final, here("temp_data", "merged_datasets", "brazil_aoi", "loss_commodity_aeay_long_final.Rdata"))

rm(final)




