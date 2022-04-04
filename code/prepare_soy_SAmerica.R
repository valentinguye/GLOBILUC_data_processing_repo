##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

neededPackages <- c("data.table", "plyr", "tidyr", "dplyr",  "Hmisc", "sjmisc", "stringr",
                    "here", "readstata13", "foreign", "readxl", "writexl",
                    "raster", "rgdal", "sp", "spdep", "sf", "stars", "gfcanalysis",  "nngeo", # "osrm", "osrmr",
                    "lubridate","exactextractr",
                    "doParallel", "foreach", "snow", 
                    "knitr", "kableExtra",
                    "DataCombine", 
                    "fixest", 
                    "boot", "fwildclusterboot", "sandwich",
                    "ggplot2", "leaflet", "tmap", "dotwhisker")

# Install them in their project-specific versions
renv::restore(packages = neededPackages)

### USER HAS TO CONFIRM BY TYPING < y > IN CONSOLE AT THIS POINT 

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
dir.create(here("temp_data", "processed_soy_SAmerica", "southam_aoi"), recursive = TRUE)
dir.create(here("temp_data", "processed_lossdrivers", "southam_aoi"), recursive = TRUE)
dir.create(here("temp_data", "processed_fc2000", "southam_aoi"), recursive = TRUE)
dir.create(here("temp_data", "processed_pasture2000", "southam_aoi"), recursive = TRUE)
dir.create(here("temp_data", "merged_datasets", "southam_aoi"), recursive = TRUE)

### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")

### GLOBAL CRS USED throughout the study ### 
mercator_world_crs <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "

### GAEZ OBJECTS
gaez_dir <- here("temp_data", "GAEZ", "v4", "AEAY_out_density",  "Rain-fed")
gaez_crops <- list.files(path = here(gaez_dir, "High-input"), 
                         pattern = "", 
                         full.names = FALSE)
gaez_crops <- gsub(pattern = ".tif", replacement = "", x = gaez_crops)

### SOUTH AMERICA AOI 
sam <- raster(here("input_data", "SouthAmerica_Soybean", "soy_SAmerica_3km_unidir.tif"))
southam_aoi <- extent(sam)
rm(sam)


## All aggregations are performed to the curtis drivers data, because this is the lowest resolution we will need
drivers <- raster(here("input_data", "curtis", "Goode_FinalClassification_19_05pcnt_prj", "Goode_FinalClassification_19_05pcnt_prj.tif"))
# first crop it to brazil aoi 
drivers <- crop(drivers, southam_aoi)
# plot(drivers)

### SOME PATHS THAT ARE CALLED THROUGHOUT THE SCRIPT 
extentsoy_resampled_output_name <- here("temp_data", "processed_soy_SAmerica", "southam_aoi", "resampled_soy_extent_0120.tif")
unidirsoy_resampled_output_name <- here("temp_data", "processed_soy_SAmerica", "southam_aoi", "resampled_soy_unidir_0220.tif")
gaez_resampled_output_name <- here(gaez_dir, "southam_resampleddrivers_high_input.tif")
fc2k_resampled_output_name <- here("temp_data", "processed_fc2000", "southam_aoi", "resampled_fc_2000.tif")
pst2k_resampled_output_name <- here("temp_data", "processed_pasture2000", "southam_aoi", "resampled_drivers_pasture_2000.tif")

any_resampled_output_name <- here("temp_data", "processed_lossdrivers", "southam_aoi", "loss_drivers_any.tif")
commodity_resampled_output_name <- here("temp_data", "processed_lossdrivers", "southam_aoi", "loss_drivers_commodity.tif")
shifting_resampled_output_name <- here("temp_data", "processed_lossdrivers", "southam_aoi", "loss_drivers_shifting.tif")
forestry_resampled_output_name <- here("temp_data", "processed_lossdrivers", "southam_aoi", "loss_drivers_forestry.tif")
fire_resampled_output_name <- here("temp_data", "processed_lossdrivers", "southam_aoi", "loss_drivers_fire.tif")
urba_resampled_output_name <- here("temp_data", "processed_lossdrivers", "southam_aoi", "loss_drivers_urba.tif")

mask_path <- here("temp_data", "processed_lossdrivers", "southam_aoi", "always_zero_mask_lossdrivers_any.tif")


#### AGGREGATE AND ALIGNE SOY EXTENTS ####
## Brick layers 

# import annual layers of forest loss (in hectares) as computed in GEE (and downloaded from Google Drive to input_data/)
soy <- brick(here("input_data", "SouthAmerica_Soybean", "soy_SAmerica_3km_extent.tif"))
names(soy)
# rename layers here (was not done in GEE because no loop there to make it conveniently)
names(soy) <- paste0("extent_soy_", seq(2001, 2020, 1))

# rasterlist <- list.files(path = here("input_data", "SouthAmerica_Soybean"), 
#                          pattern = "^soy_SAmerica_5km_", 
#                          full.names = TRUE) %>% as.list()
# soy <- brick(rasterlist)
# 
# writeRaster(soy[[1:20]], here("temp_data", "processed_soy_SAmerica", "southam_aoi", "soy0119.tif"), 
#             overwrite = TRUE)
# soy <- brick(here("temp_data", "processed_soy_SAmerica", "southam_aoi", "soy0119.tif"))


# define output file name
aggr_output_name <- here("temp_data", "processed_soy_SAmerica", "southam_aoi", "aggr_soy_extent_0120.tif")

# aggregate it from the ~5km cells to ~10km
raster::aggregate(soy, fact = c(res(drivers)[1]/res(soy)[1], res(drivers)[2]/res(soy)[2]),
                  expand = FALSE,
                  fun = sum,
                  na.rm = FALSE, # NA values are on the eastern band only. Thus they can contaminate aggregation safely. 
                  filename = aggr_output_name,
                  overwrite = TRUE)

# align to DRIVERS exactly 
aggregated_ssa <- brick(aggr_output_name)

resample(x = aggregated_ssa, 
         y = drivers, 
         method = "bilinear", # bilinear or ngb changes nothing 
         filename = extentsoy_resampled_output_name, 
         overwrite = TRUE)




#### AGGREGATE AND ALIGNE SOY EXPANSION (UNIDIR) #### 
soy <- brick(here("input_data", "SouthAmerica_Soybean", "soy_SAmerica_3km_unidir.tif"))
names(soy)
# keeps years from 2002 to 2020 (i.e. all those in the brick)
# 2002 because expansion data is available only as of this year, and 2020 because although 2020 is not available in deforestation data, we want max information in final dataset 


# values are in hectares of soy in 5km x 5km grid cells - and THERE ARE NAs 
ssa2002 <- soy$unidir_soy_2002
plot(ssa2002)
valpast01 <- values(ssa2002)
summary(valpast01)
anyNA(valpast01)

# define output file name
aggr_output_name <- here("temp_data", "processed_soy_SAmerica", "southam_aoi", "aggr_soy_unidir_0220.tif")

# aggregate it from the ~5km cells to ~10km
raster::aggregate(soy, fact = c(res(drivers)[1]/res(soy)[1], res(drivers)[2]/res(soy)[2]),
                  expand = FALSE,
                  fun = sum,
                  na.rm = FALSE, # NA values are on the eastern band only. Thus they can contaminate aggregation safely. 
                  filename = aggr_output_name,
                  overwrite = TRUE)

# align to DRIVERS exactly 
aggregated_ssa <- brick(aggr_output_name)

resample(x = aggregated_ssa, 
         y = drivers, 
         method = "bilinear", # bilinear or ngb changes nothing 
         filename = unidirsoy_resampled_output_name, 
         overwrite = TRUE)




#### ALIGNE GAEZ ####
gaez <- brick(here(gaez_dir, "high_input_all.tif"))
gaez <- crop(gaez, southam_aoi)

# gaez$high_input_all.1 %>% values() %>% summary()

# resample directly (without aggregating first) from the ~9km cells to ~10km (aggregate does not work bc resolutions are to close)

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

#### 2000 FOREST COVER #### 
fc2k <- raster(here("input_data", "fc_2000_3km_10th.tif"))

fc2k <- crop(fc2k, southam_aoi)

fc2k_aggr_output_name <- here("temp_data", "processed_fc2000", "southam_aoi", "aggr_fc_2000.tif")

aggregate(fc2k, fact = c(res(drivers)[1]/res(fc2k)[1], res(drivers)[2]/res(fc2k)[2]),
          expand = FALSE,
          fun = sum,
          na.rm = TRUE,  # no NA values a priori mais bon 
          filename = fc2k_aggr_output_name,
          overwrite = TRUE)

# align to DRIVERS exactly 
aggregated_fc2000 <- raster(fc2k_aggr_output_name)

resample(x = aggregated_fc2000, 
         y = drivers, 
         method = "ngb", # ngb because bilinear yields negatie values
         filename = fc2k_resampled_output_name, 
         overwrite = TRUE)

#### 2000 PASTURE SHARE ####
pst2k <- raster(here("input_data", "CroplandPastureArea2000_Geotiff", "Pasture2000_5m.tif"))
pst2k <- crop(pst2k, southam_aoi)
# resample directly (without aggregating first) from the ~9km cells to ~10km (aggregate does not work bc resolutions are to close)

# transform NAs to 0, such that the more numerous NAs in pasture data do not force losing information when the stack is turned to a data frame. 
pst2k <- reclassify(pst2k, cbind(NA, 0))

resample(x = pst2k, 
         y = drivers, 
         method = "ngb", # we use ngb and not bilinear because the output values' summary better fits that of the aggregated layer 
         # and the bilinear interpolation arguably smoothes the reprojection more than necessary given that from and to are already very similar.  
         filename = pst2k_resampled_output_name, 
         overwrite = TRUE)


#### DRIVEN DEFORESTATION DATA ####
### ANY ### 
# it is already aggregated at drivers resolution. Just crop it to the brazil aoi, and aline it to the previous rasters
any_lossdrivers <- brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", "loss_drivers_any.tif"))

any_lossdrivers <- crop(any_lossdrivers, southam_aoi)

resample(x = any_lossdrivers, 
         y = drivers, 
         method = "ngb", # no difference between bilinear and ngb
         filename = any_resampled_output_name, 
         overwrite = TRUE)

### COMMODITY ###
# it is already aggregated at drivers resolution. Just crop it to the brazil aoi, and aline it to the previous rasters
commodity_lossdrivers <- brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", "loss_drivers_commodity.tif"))

commodity_lossdrivers <- crop(commodity_lossdrivers, southam_aoi)

resample(x = commodity_lossdrivers, 
         y = drivers, 
         method = "ngb", # no difference between bilinear and ngb
         filename = commodity_resampled_output_name, 
         overwrite = TRUE)

### SHIFTING ###
# it is already aggregated at drivers resolution. Just crop it to the brazil aoi, and aline it to the previous rasters
shifting_lossdrivers <- brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", "loss_drivers_shifting.tif"))

shifting_lossdrivers <- crop(shifting_lossdrivers, southam_aoi)

resample(x = shifting_lossdrivers, 
         y = drivers, 
         method = "ngb", # no difference between bilinear and ngb
         filename = shifting_resampled_output_name, 
         overwrite = TRUE)

### FORESTRY ###
# it is already aggregated at drivers resolution. Just crop it to the brazil aoi, and aline it to the previous rasters
forestry_lossdrivers <- brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", "loss_drivers_forestry.tif"))

forestry_lossdrivers <- crop(forestry_lossdrivers, southam_aoi)

resample(x = forestry_lossdrivers, 
         y = drivers, 
         method = "ngb", # no difference between bilinear and ngb
         filename = forestry_resampled_output_name, 
         overwrite = TRUE)

### FIRE ###
# it is already aggregated at drivers resolution. Just crop it to the brazil aoi, and aline it to the previous rasters
fire_lossdrivers <- brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", "loss_drivers_fire.tif"))

fire_lossdrivers <- crop(fire_lossdrivers, southam_aoi)

resample(x = fire_lossdrivers, 
         y = drivers, 
         method = "ngb", # no difference between bilinear and ngb
         filename = fire_resampled_output_name, 
         overwrite = TRUE)

### URBANIZATION ###
# it is already aggregated at drivers resolution. Just crop it to the brazil aoi, and aline it to the previous rasters
urba_lossdrivers <- brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", "loss_drivers_urba.tif"))

urba_lossdrivers <- crop(urba_lossdrivers, southam_aoi)

resample(x = urba_lossdrivers, 
         y = drivers, 
         method = "ngb", # no difference between bilinear and ngb
         filename = urba_resampled_output_name, 
         overwrite = TRUE)


#### PREPARE ALWAYS ZERO DEFORESTATION MASK #### 
any_lossdrivers <- brick(any_resampled_output_name)

always_zero <- function(y){if_else(condition = (sum(y)==0), true = 0, false = 1)}

overlay(x = any_lossdrivers,
        fun = always_zero,
        filename = mask_path,
        na.rm = TRUE,
        overwrite = TRUE)



#### STACK AND MASK SOY MAPS, GAEZ, AND COMMODITY DEFORESTATION #### 
# Read layers to be stacked
drivenloss_commodity <- brick(commodity_resampled_output_name)
drivenloss_shifting <- brick(shifting_resampled_output_name)
drivenloss_forestry <- brick(forestry_resampled_output_name)
drivenloss_fire <- brick(fire_resampled_output_name)
drivenloss_urba <- brick(urba_resampled_output_name)

extent_soy <- brick(extentsoy_resampled_output_name)
unidir_soy <- brick(unidirsoy_resampled_output_name)

gaez <- brick(gaez_resampled_output_name)
fc2k <- raster(fc2k_resampled_output_name)
pst2k <- raster(pst2k_resampled_output_name)

# It is important to explicitly rename layers that are going to be stacked and then called to reshape the data frame 
# for time varying variables, the dot is important. 
names(drivenloss_commodity) <- paste0("driven_loss_commodity.",seq(2001, 2019, 1)) 
names(drivenloss_shifting) <- paste0("driven_loss_shifting.",seq(2001, 2019, 1)) # note the difference with the names of phtfloss (not the same years)
names(drivenloss_forestry) <- paste0("driven_loss_forestry.",seq(2001, 2019, 1)) # note the difference with the names of phtfloss (not the same years)
names(drivenloss_fire) <- paste0("driven_loss_fire.",seq(2001, 2019, 1)) # note the difference with the names of phtfloss (not the same years)
names(drivenloss_urba) <- paste0("driven_loss_urba.",seq(2001, 2019, 1)) # note the difference with the names of phtfloss (not the same years)

names(extent_soy) <- paste0("extent_soy.",seq(2001, 2020, 1)) # note the difference with the names of phtfloss (not the same years)
names(unidir_soy) <- paste0("unidir_soy.",seq(2002, 2020, 1)) # note the start in 2002

names(gaez) <- gaez_crops

names(fc2k) <- "fc_2000"
names(pst2k) <- "pasture_share_2000"

## STACK 
southam_stack <- stack(drivenloss_commodity,
                              drivenloss_shifting, 
                              drivenloss_forestry, 
                              drivenloss_fire,
                              drivenloss_urba,
                              extent_soy, 
                              unidir_soy,
                              gaez, fc2k, pst2k)
names(southam_stack)

# # save the stack at this point, as a clean export 
# writeRaster(southam_stack,
#             filename = here("temp_data", "merged_datasets", "southam_aoi", "drivenloss_soy_gaez"),
#             overwrite = TRUE, 
#             format = "raster")


### MASK TO REMOVE ALWAYS ZERO PIXELS AND LIGHTEN THE DATA FRAMES ### 
mask <- raster(mask_path)

southam_stack <- mask(x = southam_stack, 
                     mask = mask,
                     maskvalue = 0, # necessary here, because the there is no NA in the mask, only 0 and 1 (see the prepare_loss_drivers.R script)
                     updatevalue = NA)


# (note that masking changes the summary values)


### RASTER TO DATAFRAME ### 

# na.rm = TRUE is key here, as it removes previously masked pixels (NA) and ensures the output is not too large (memory intensive)
# We also set long to false because we reshape with a proper function for more control
wide_df <- raster::as.data.frame(southam_stack, na.rm = TRUE, xy = TRUE, centroids = TRUE, long = FALSE) # ~700s. 
# set na.rm = FALSE and run below until saveRDS to extract driverloss_all_aeay_long_allcells.Rdata

# Rename coordinate variables
names(wide_df)
head(wide_df[,c("x", "y")])
wide_df <- dplyr::rename(wide_df, lon = x, lat = y)


# SPECIAL TO SOY SCRIPT HERE: ADD 2001 COLUMN FOR THE SOY EXPANSION  
# because we want 2001 extent to be in the final data set, but expansion (unidir) is only avaiable from 2002 and number of years must be 
# equal for all varying vars
wide_df$unidir_soy.2001 <- NA
wide_df$driven_loss_commodity.2020 <- NA
wide_df$driven_loss_shifting.2020 <- NA
wide_df$driven_loss_forestry.2020 <- NA
wide_df$driven_loss_fire.2020 <- NA
wide_df$driven_loss_urba.2020 <- NA

### WIDE TO LONG ### 

# Since we merged datasets in the raster format, we wont need a lonlat format id for each grid cell. 
# So we can simply create an ID that's a sequence. 
wide_df$grid_id <- seq(1, nrow(wide_df), 1) 

# the dot is, by construction of all variable names, only in the names of time varying variables. 
# fixed = TRUE is necessary (otherwise the dot is read as a regexp I guess)
# Note also that it is important that it is structured in a LIST when there are several varying variables in the *long* format
# Because: "Notice that the order of variables in varying is like x.1,y.1,x.2,y.2."
varying_vars <- list(paste0("driven_loss_commodity.", seq(2001, 2020, 1)),
                     paste0("driven_loss_shifting.", seq(2001, 2020, 1)),
                     paste0("driven_loss_forestry.", seq(2001, 2020, 1)),
                     paste0("driven_loss_fire.", seq(2001, 2020, 1)), 
                     paste0("driven_loss_urba.", seq(2001, 2020, 1)),
                     paste0("extent_soy.",seq(2001, 2020, 1)), 
                     paste0("unidir_soy.",seq(2001, 2020, 1)))
#varying_vars <- names(driverloss_gaez)[grep(".", names(driverloss_gaez), fixed = TRUE)]

# reshape to long.
long_df <- stats::reshape(wide_df,
                          varying = varying_vars,
                          v.names = c("driven_loss_commodity", "driven_loss_shifting", "driven_loss_forestry", "driven_loss_fire", "driven_loss_urba", 
                                      "extent_soy", "unidir_soy"),
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

saveRDS(long_df, here("temp_data", "merged_datasets", "southam_aoi", "driverloss_all_aeay_long.Rdata"))

rm(long_df, varying_vars, southam_stack, gaez, mask, commodity_lossdrivers, fc2k, pst2k)


#### COUNTRY VARIABLE #### 
countries <- st_read(here("input_data", "Global_LSIB_Polygons_Detailed"))
length(unique(countries$COUNTRY_NA)) == nrow(countries)

path <- here("temp_data", "merged_datasets", "southam_aoi", "driverloss_all_aeay_long.Rdata")
df <- readRDS(path)

# Remove gaez variables
df <- dplyr::select(df,-all_of(gaez_crops))

# Use cross section only
df_cs <- df[!duplicated(df$grid_id),]

rm(df)

# Spatial
df_cs <- st_as_sf(df_cs, coords = c("lon", "lat"), crs = 4326, remove = FALSE)


# transform shapes such that st_nearest_feature works with projected data
df_cs <- st_transform(df_cs, crs = mercator_world_crs)
countries <- st_transform(countries, crs = mercator_world_crs)

# This is much much faster (like 5 minutes vs. 6h).
# df_cs <- st_join(x = countries[,c("OBJECTID", "COUNTRY_NA")],
#                  y = df_cs,
#                  join = st_contains,
#                  prepared = TRUE,
#                  left = FALSE)# performs inner join so returns only records that spatially match.


# However, we use st_nearest_feature so that all points match a country
df_cs <- st_join(x = df_cs,
                 y = countries,
                 join = st_nearest_feature,
                 left = TRUE)

rm(countries)

# names(df_cs)[names(df_cs) == "OBJECTID"] <- "country_id"
names(df_cs)[names(df_cs) == "COUNTRY_NA"] <- "country_name"

df_cs <- st_drop_geometry(df_cs)

# Keep only new variable and id
df_cs <- df_cs[,c("grid_id", "country_name")]

# saveRDS(df_cs, path)
saveRDS(df_cs, here("temp_data", "merged_datasets", "southam_aoi", "driverloss_all_aeay_cs_country_nf.Rdata"))
rm(df_cs)


#### BIGGER CELL VARIABLES #### 
## Prepare base data
path <- here("temp_data", "merged_datasets", "southam_aoi", "driverloss_all_aeay_long.Rdata")
df <- readRDS(path)

# Remove gaez variables
df <- dplyr::select(df,-all_of(gaez_crops))
# Use cross section only
df_cs <- df[!duplicated(df$grid_id),]
# Spatial
df_cs <- st_as_sf(df_cs, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

rm(df)

## Prepare bigger square grids
# for ~50km grid cells (5 times larger grid cells in both dimensions, hence 25 times larger)
bigger_50km <- aggregate(drivers, fact = 5, expand = TRUE, fun = sum)
bigger_50km_stars <- st_as_stars(bigger_50km)
bigger_50km_sf <- st_as_sf(bigger_50km_stars, as_points = FALSE, merge = FALSE, na.rm = FALSE, long = FALSE)

bigger_50km_sf$grid_id_50km <- seq(from = 1, to = nrow(bigger_50km_sf))
bigger_50km_sf <- bigger_50km_sf[,c("grid_id_50km", "geometry")]

# repeat for ~100km grid cells (10 times larger grid cells in both dimensions, hence 100 times larger)
bigger_100km <- aggregate(drivers, fact = 10, expand = TRUE, fun = sum)
bigger_100km_stars <- st_as_stars(bigger_100km)
bigger_100km_sf <- st_as_sf(bigger_100km_stars, as_points = FALSE, merge = FALSE, na.rm = FALSE, long = FALSE)

bigger_100km_sf$grid_id_100km <- seq(from = 1, to = nrow(bigger_100km_sf))
bigger_100km_sf <- bigger_100km_sf[,c("grid_id_100km", "geometry")]

# repeat for ~500km grid cells (50 times larger grid cells in both dimensions, hence 2500 times larger)
bigger_500km <- aggregate(drivers, fact = 50, expand = TRUE, fun = sum)
bigger_500km_stars <- st_as_stars(bigger_500km)
bigger_500km_sf <- st_as_sf(bigger_500km_stars, as_points = FALSE, merge = FALSE, na.rm = FALSE, long = FALSE)

bigger_500km_sf$grid_id_500km <- seq(from = 1, to = nrow(bigger_500km_sf))
bigger_500km_sf <- bigger_500km_sf[,c("grid_id_500km", "geometry")]

# DO NOT TRANSFORM, because it makes the inner_join associate two bigger squares for each smaller one, 
# I don't really know why, but it does not do so with geographic coordinates, and there is no mismatch, and no problem of imprecision due to 
# geographic coordinates being used as planer because aoi is not near the pole. 
# df_cs <- st_transform(df_cs, crs = mercator_world_crs)
# bigger_50km_sf <- st_transform(bigger_50km_sf, crs = mercator_world_crs)
# bigger_100km_sf <- st_transform(bigger_100km_sf, crs = mercator_world_crs)

# join the variables
df_cs <- st_join(x = bigger_50km_sf,
                 y = df_cs,
                 join = st_contains,
                 prepared = FALSE, # tests shows that this changes nothing, whether shapes are transformed or not
                 left = FALSE)# performs inner join so returns only records that spatially match.

df_cs <- st_join(x = bigger_100km_sf,
                 y = df_cs,
                 join = st_contains,
                 prepared = FALSE, # tests shows that this changes nothing, whether shapes are transformed or not
                 left = FALSE)# performs inner join so returns only records that spatially match.

df_cs <- st_join(x = bigger_500km_sf,
                 y = df_cs,
                 join = st_contains,
                 prepared = FALSE, # tests shows that this changes nothing, whether shapes are transformed or not
                 left = FALSE)# performs inner join so returns only records that spatially match.


df_cs <- st_drop_geometry(df_cs)

length(unique(df_cs$grid_id)) == nrow(df_cs)

# Keep only new variable and id
df_cs <- df_cs[,c("grid_id", "grid_id_50km", "grid_id_100km", "grid_id_500km")]

saveRDS(df_cs, here("temp_data", "merged_datasets", "southam_aoi", "driverloss_all_aeay_cs_biggercells.Rdata"))

rm(df_cs, bigger_50km_sf, bigger_100km_sf, bigger_500km_sf)

#### ADD REMAINING FOREST VARIABLE #### 
# neither country nor continent information is relevant here. 

df <- readRDS(here("temp_data", "merged_datasets", "southam_aoi", "driverloss_all_aeay_long.Rdata"))

# Remove gaez variables
df <- dplyr::select(df,-all_of(gaez_crops))

year_list <- list()

df <- mutate(df, 
             driven_loss_any = driven_loss_commodity + driven_loss_shifting + driven_loss_forestry + driven_loss_fire + driven_loss_urba)

# for some grid cells, deforestation from different drivers is positive in the same year. 
# This probably comes from resampling/reprojecting operations, and can be observed here only when forest loss was positive in those places
# dplyr::filter(df, driven_loss_commodity>0 & driven_loss_any > driven_loss_commodity)
# df[df$grid_id == 110,]

# just ensure that driven_loss_any is not counting 4 times too much deforestation in those places. 
df[df$driven_loss_commodity>0 & df$driven_loss_any > df$driven_loss_commodity, ] <- dplyr::filter(df, 
                                                                                                  driven_loss_commodity>0 & driven_loss_any > driven_loss_commodity) %>% 
  mutate(driven_loss_any = driven_loss_commodity)


# in the first year (2001), the past year accumulated deforestation is null. 
year_list[["2001"]] <- df[df$year == 2001, c("grid_id", "year")] 
year_list[["2001"]][,"accu_defo_since2k"] <- 0

# then, each year's deforestation accumulated in the past is the sum of *past years'* deforestation
years <- 2008:2011 # we need it only 2007-2010, because we possibly use only those years as starting years not 2002:max(df$year)
for(y in years){
  sub_ <- df[df$year < y,]
  year_list[[as.character(y)]] <- ddply(sub_, "grid_id", summarise,
                                        accu_defo_since2k = sum(driven_loss_any, na.rm = TRUE))
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

saveRDS(remaining, here("temp_data", "merged_datasets", "southam_aoi", "driverloss_all_aeay_long_remaining.Rdata"))

rm(year_list, sub_, accu_defo_df, df)


#### GROUP AND STANDARDIZE AEAY VARIABLES #### 
# all groupings in this section are motivated on the GAEZ v4 model documentation, and in particular Table A4-1.3

df <- readRDS(here("temp_data", "merged_datasets", "southam_aoi", "driverloss_all_aeay_long.Rdata"))

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

### STANDARDIZE ### 

# divide each revenue variable by the sum of them, to standardize.
# To understand this line, see https://dplyr.tidyverse.org/articles/rowwise.html#row-wise-summary-functions
df_cs <- dplyr::mutate(df_cs, eaear_sum = rowSums(across(.cols = any_of(eaear2std))))

df_cs <- dplyr::mutate(df_cs, across(.cols = any_of(eaear2std),
                                     .fns = ~./eaear_sum, 
                                     .names = paste0("{.col}", "_std"))) 

df_cs <- dplyr::select(df_cs, -eaear_sum)

## Second way to standardize: for each crop that can be matched with a price, 
# standardize by dividing by the sum of the suitability indexes of the N (N = 1,2) crops with the highest suitability (among all, not only among the six drivers), 
# and give a 0 value to the crops that are not in the top N suitability index. 
# if N = 1, this procedure is equivalent to sj = 1[Sj = max(Si)]

# this code is intricate but it handles potential but unlikely cases where two crops have equal EAEAR

# identify the highest suitability index values (in every grid cell)
df_cs <- df_cs %>% rowwise(grid_id) %>% dplyr::mutate(max_eaear = max(c_across(cols = any_of(eaear2std)))) %>% as.data.frame()

# and for the alternative set of crops
df_cs <- df_cs %>% rowwise(grid_id) %>% dplyr::mutate(max_eaearbis = max(c_across(cols = any_of(eaear2std_bis)))) %>% as.data.frame()

## N = 1
# if N = 1, this procedure is equivalent to sj = 1[Sj = max(Si)]
df_cs <- dplyr::mutate(df_cs, across(.cols = any_of(eaear2std),
                                     .fns = ~if_else(.==max_eaear, true = 1, false = 0), 
                                     .names = paste0("{.col}", "_ismax")))

# and then standardize by the number of different crops being the highest
all_crops_ismax <- paste0(eaear2std,"_ismax")
df_cs <- dplyr::mutate(df_cs, n_max = rowSums(across(.cols = (any_of(all_crops_ismax)))))

unique(df_cs$n_max) # 16 is when GAEZ is NA
# df_cs[df_cs$n_max==16,]

df_cs <- dplyr::mutate(df_cs, across(.cols = (any_of(all_crops_ismax)),
                                     .fns = ~./n_max, 
                                     .names = paste0("{.col}", "_std1")))

# remove those columns
df_cs <- dplyr::select(df_cs, !ends_with("_ismax"))  

# rename new ones 
names(df_cs)[grepl("_std1", names(df_cs))] <- paste0(eaear2std, "_std1")

# _std do sum up to 1. 
# df_cs[87687,paste0(eaear2std, "_std2")]%>%sum()

## N = 2: 2nd highest:
# helper function to apply within dplyr rowwise framework
max2nd <- function(x){max(x[x!=max(x)])}

df_cs <- df_cs %>% rowwise() %>% mutate(max_eaear_2nd = max2nd(c_across(cols = (any_of(eaear2std))))) %>% as.data.frame()
# this returns a warning about aucun argument pour max ; -Inf est renvoyé" when there are only zero values for every crop in the grid cell
# not a problem

# new column for each crop, telling whether it's in the top 2 or not
df_cs <- dplyr::mutate(df_cs, across(.cols = any_of(eaear2std),
                                     .fns = ~if_else(.>=max_eaear_2nd, true = 1, false = 0), 
                                     .names = paste0("{.col}", "_istop2")))

all_crops_istop2 <- paste0(eaear2std,"_istop2")

df_cs <- dplyr::mutate(df_cs, n_top2 = rowSums(across(.cols = (any_of(all_crops_istop2)))))
unique(df_cs$n_top2) # 16 is when GAEZ is NA
# df_cs[df_cs$n_max==16,]
# in other words there is always only 2 crops in the top 2     

# multiply istop2 columns with their corresponding EAEAR columns
for(crop in eaear2std){
  df_cs <- mutate(df_cs, !!as.symbol(paste0(crop, "_istop2")) := !!as.symbol(paste0(crop, "_istop2")) * !!as.symbol(crop) )
}

# sum them, and standardize with this sum
df_cs <- dplyr::mutate(df_cs, sum_top2 = rowSums(across(.cols = (any_of(all_crops_istop2)))))

df_cs <- dplyr::mutate(df_cs, across(.cols = (any_of(all_crops_istop2)),
                                     .fns = ~./sum_top2, 
                                     .names = paste0("{.col}", "_std2")))

# remove those columns
df_cs <- dplyr::select(df_cs, !ends_with("_istop2"))  
# rename columns
names(df_cs)[grepl("_istop2_std2", names(df_cs))] <- paste0(eaear2std, "_std2")  


### ### ### 
## Repeat std2 for the broader set of crops 
df_cs <- df_cs %>% rowwise() %>% mutate(max_eaear_2nd_bis = max2nd(c_across(cols = (any_of(eaear2std_bis))))) %>% as.data.frame()

# new column for each crop, telling whether it's in the top 2 or not
df_cs <- dplyr::mutate(df_cs, across(.cols = any_of(eaear2std_bis),
                                     .fns = ~if_else(.>=max_eaear_2nd_bis, true = 1, false = 0), 
                                     .names = paste0("{.col}", "_istop2")))

all_crops_istop2_bis <- paste0(eaear2std_bis,"_istop2")

# multiply istop2 columns with their corresponding EAEAR columns
for(crop in eaear2std_bis){
  df_cs <- mutate(df_cs, !!as.symbol(paste0(crop, "_istop2")) := !!as.symbol(paste0(crop, "_istop2")) * !!as.symbol(crop) )
}

# sum them, and standardize with this sum
df_cs <- dplyr::mutate(df_cs, sum_top2_bis = rowSums(across(.cols = (any_of(all_crops_istop2_bis)))))

df_cs <- dplyr::mutate(df_cs, across(.cols = (any_of(all_crops_istop2_bis)),
                                     .fns = ~./sum_top2_bis, 
                                     .names = paste0("{.col}", "_std2bis")))

# remove those columns
df_cs <- dplyr::select(df_cs, !ends_with("_istop2"))  
# rename columns
names(df_cs)[grepl("_istop2_std2bis", names(df_cs))] <- paste0(eaear2std_bis, "_std2bis")  

### ### ### ### ### 

# Select variables to save: all the eaear variables, standardized or not
df_cs <- dplyr::select(df_cs, -max_eaear_2nd, - max_eaear_2nd_bis)
var_names <- grep(pattern = "eaear_", names(df_cs), value = TRUE) 
df_cs <- df_cs[,c("grid_id", var_names)]

saveRDS(df_cs, here("temp_data", "merged_datasets", "southam_aoi",  "driverloss_all_aeay_cs_stdeaear.Rdata"))

rm(df_cs)

#### MERGE ADDITIONAL VARIABLES ####  
# Base dataset (including outcome variable(s))
df_base <- readRDS(here("temp_data", "merged_datasets", "southam_aoi", "driverloss_all_aeay_long.Rdata"))


## COUNTRY
# just compute country and continent variables, even if invariant, so they can be called in generic function
df_country <- readRDS(here("temp_data", "merged_datasets", "southam_aoi", "driverloss_all_aeay_cs_country_nf.Rdata"))

# Merge them and remove to save memory 
final <- left_join(df_base, df_country, by = "grid_id")
rm(df_base, df_country)
final <- mutate(final, country_year = paste0(country_name, "_", year))

# just compute country and continent variables, even if invariant, so they can be called in generic function
final$continent_name <- "America"
# Create country year fixed effect

## BIGGER CELLS
df_biggercells <- readRDS(here("temp_data", "merged_datasets", "southam_aoi", "driverloss_all_aeay_cs_biggercells.Rdata"))

final <- left_join(final, df_biggercells, by = "grid_id")
rm(df_biggercells)

# Create bigger cell-year identifier
final <- mutate(final, grid_id_50km_year = paste0(grid_id_50km, "_", year))
final <- mutate(final, grid_id_100km_year = paste0(grid_id_100km, "_", year))
final <- mutate(final, grid_id_500km_year = paste0(grid_id_500km, "_", year))
# length(unique(final$grid_id_50km_year))==length(unique(final$grid_id_50km))*length(unique(final$year))


## EAEAR
df_stdeaear <- readRDS(here("temp_data", "merged_datasets", "southam_aoi",  "driverloss_all_aeay_cs_stdeaear.Rdata"))  

final <- left_join(final, df_stdeaear, by = "grid_id") # no issue with using grid_id as a key here, bc df_remain was computed just above from the df_base data
rm(df_stdeaear)

## REMAINING
df_remain <- readRDS(here("temp_data", "merged_datasets", "southam_aoi", "driverloss_all_aeay_long_remaining.Rdata"))

final <- left_join(final, df_remain, by = c("grid_id", "year"))  # no issue with using grid_id as a key here, bc df_remain was computed just above from the df_base data
rm(df_remain)

saveRDS(final, here("temp_data", "merged_datasets", "southam_aoi", "driverloss_all_aeay_long_final.Rdata"))

rm(final)














