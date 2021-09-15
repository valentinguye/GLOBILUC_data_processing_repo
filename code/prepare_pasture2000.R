##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages = c("dplyr", "here", "readstata13", "DescTools",
                   "raster", "rgdal", "sp", "sf",
                   "doParallel", "foreach", "parallel")

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
dir.create(here("temp_data", "processed_pasture2000", "tropical_aoi"), recursive = TRUE)

### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")


# Import the pasture data
pasture <- raster(here("input_data", "CroplandPastureArea2000_Geotiff", "Pasture2000_5m.tif"))

croped_path <- here("temp_data", "processed_pasture2000", "tropical_aoi", "tropical_aoi_pasture2000.tif")
# crop to ***TROPICAL*** AOI 
ext <- extent(c(-180, 179.9167, -30, 30))
pasture <- crop(x = pasture, 
                y = ext, 
                filename = croped_path, 
                overwrite = TRUE)

croped <- raster(croped_path)

# Align to GAEZ exactly
# resolution is already similar, if not exactly equal

# Read in a GAEZ grid, because we want to align to it. 
gaez <- raster(here("temp_data", "GAEZ", "AES_index_value", "Rain-fed", "High-input", "Alfalfa.tif"))

resampled_path <- here("temp_data", "processed_pasture2000", "tropical_aoi", "resampled_pasture2000.tif")

resample(x = croped, 
         y = gaez, 
         method = "ngb", # we use ngb and not bilinear because the output values' summary better fits that of the aggregated layer 
         # and the bilinear interpolation arguably smoothes the reprojection more than necessary given that from and to are already very similar.  
         filename = resampled_path, 
         overwrite = TRUE)




# r <- values(resampled)
# r1 <- r[!is.na(r)]
# summary(r1)
# 
# it's important to compare with croped raster and not the whole pasture (global) raster that includes out of aoi cells
# v <- values(croped)
# v1 <- v[!is.na(v)]
# summary(v1)


### MASK PASTURE TO REMOVE ALWAYS ZERO PIXELS AND LIGHTEN THE DATA FRAMES ### 
resampled <- raster(resampled_path)

mask <- raster(here("temp_data", "processed_lossdrivers", "tropical_aoi", "always_zero_mask_lossdrivers.tif"))

masked_pasture_output_name <- here("temp_data", "processed_pasture2000", "tropical_aoi", "driverloss_masked_pasture.tif")

mask(x = resampled, 
     mask = mask, 
     maskvalue = 0, # necessary here, because there is no NA in the mask, only 0 and 1 (see the prepare_loss_drivers.R script)
     updatevalue = NA, 
     filename = masked_pasture_output_name,
     overwrite = TRUE)


# The pasture data have a bit more NAs (~30000) than the gaez data. 
# This might be a question of what they respectively counted as the shore. 
# Anyways, in order to merge the pasture data, we need either to trim the other 
# data sets to the data available (i.e. not NA) in pasture (and loose information), 
# or to merge based on the centroids

### STACK RASTERS TO MERGE ###
# we do not stack rasters to merge them here. It is just for one variable. 

### RASTER TO DATAFRAME ### 
pasture_m <- raster(masked_pasture_output_name)

# na.rm = TRUE is key here, as it removes previously masked pixels (NA) and ensures the output is not too large (memory intensive)
# We also set long to false because we reshape with a proper function for more control
wide_df <- raster::as.data.frame(pasture_m, na.rm = TRUE, xy = TRUE, centroids = TRUE, long = FALSE) # ~700s. 

# Rename coordinate variables
names(wide_df)
head(wide_df[,c("x", "y")])
wide_df <- dplyr::rename(wide_df, lon = x, lat = y)

# We do not creat a sequantial ID, as this would be confounding (see explanation above)

saveRDS(wide_df, here("temp_data", "processed_pasture2000", "tropical_aoi", "pasture_4_driverloss_df.Rdata"))

rm(pasture_m, mask, wide_df)