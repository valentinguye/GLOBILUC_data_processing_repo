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
                    "ggplot2", "leaflet", "dotwhisker")

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


### THIS IS THE IMPORTANT LINE DEFINING THIS WHOLE SCRIPT: 
ext <- extent(c(-95, -25, -57, 15))




# Renaming layers with the following name mapping matrix
mapmat_data <- c("alf", "Alfalfa", 
                 "ban", "Banana",
                 "bck", "Buckwheat", 
                 "brl", "Barley",
                 "bsg", "Sorghumbiomass", 
                 "cab", "Cabbage",
                 "car", "Carrot", 
                 "chk", "Chickpea",
                 "cit", "Citrus", 
                 "coc", "Cocoa",
                 "cof", "Coffee", 
                 "con", "Coconut",
                 "cot", "Cotton", 
                 "cow", "Cowpea",
                 "csv", "Cassava", 
                 "flx", "Flax",
                 "fml", "Foxtailmillet", 
                 "grd", "Groundnut",
                 "grm", "Gram", 
                 "jtr", "Jatropha",
                 "mis", "Miscanthus", 
                 "mlt", "Millet",
                 "mze", "Maizegrain",
                 "mzs", "Maizesilage",
                 "nap", "Napiergrass" ,
                 "oat", "Oat",
                 "olp", "Oilpalm" ,
                 "olv", "Olive",
                 "oni", "Onion" ,
                 "pea", "Drypea",
                 "phb", "Phaseolousbean" ,
                 "pig", "Pigeonpea",
                 "pml", "Pearlmillet" ,
                 "pst", "Pasture",
                 "rcd", "Drylandrice",
                 "rcg", "Reedcanarygrass",
                 "rcw", "Wetlandrice" ,
                 "rsd", "Rapeseed",
                 "rub", "Rubber" ,
                 "rye", "Rye",
                 "sfl", "Sunflower" ,
                 "soy", "Soybean",
                 "spo", "Sweetpotato" ,
                 "srg", "Sorghum",
                 "sub", "Sugarbeet" ,
                 "suc", "Sugarcane",
                 "swg", "Switchgrass" ,
                 "tea", "Tea",
                 "tob", "Tobacco" ,
                 "tom", "Tomato",
                 "whe", "Wheat" ,
                 "wpo", "Whitepotato",
                 "yam", "Yam")

# 3 of these are not in the data downloaded (but that does matter): Millet, Maizesilage, and Pasture
mapmat <- matrix(data = mapmat_data, 
                 nrow = length(mapmat_data)/2,
                 ncol = 2, 
                 byrow = TRUE)

colnames(mapmat) <- c("abrev", "Names")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### PREPARE GAEZ #### 

### PREPARE SUITABILITY INDICES
# dir.create(here("temp_data", "GAEZ", "South_America", "AES_index_value", "Rain-fed", "High-input"), recursive = TRUE) 
# dir.create(here("temp_data", "GAEZ", "South_America", "AES_index_value_current", "Rain-fed", "High-input"), recursive = TRUE) 

datadir <- here("input_data", "GAEZ", "v4", "AES_index_value", "Rain-fed", "High-input")
targetdir <- here("temp_data", "GAEZ", "South_America", "v4", "AES_index_value", "Rain-fed", "High-input")
tmpdir <- here("temp_data", "tmp")


if (!dir.exists(tmpdir)) dir.create(tmpdir, recursive = TRUE)
if (dir.exists(targetdir)) {
  file.remove(list.files(path = targetdir,
                         pattern = ".tif", full.names = TRUE))
} else dir.create(targetdir, recursive = TRUE)



files <- list.files(path = datadir, pattern = ".tif")
crops <- unlist(strsplit(files, split = ".tif"))
for(crop in crops){
  crops[grepl(crop, crops)] <- mapmat[grepl(crop, paste0("sxHr_", mapmat[,"abrev"])),"Names"]
}


## Import most crops
for (j in 1:length(files)) {
  #if (any(crops[j] == cropsToAggregate)) next
  print(files[j])
  unzip(zipfile = here(datadir, files[j]), exdir = tmpdir)
  dt <- raster(paste0(tmpdir, "/data.asc"))
  
  # crop to south america AOI
  dt_trop <- crop(dt, ext)
  
  # A few points are -0.09 with no apparent reason. 
  dt_trop[dt_trop<0] <- NA 
  
  names(dt_trop) <- crops[j]
  writeRaster(dt_trop,
              filename = here(targetdir, paste0(crops[j], ".tif")),
              overwrite = TRUE)
}
rm(dt, dt_trop)

### CHECK GRASS ### 
# dt <- raster(here("input_data", "GAEZ", "Agro_climatically_attainable_yield", "Rain-fed", "High-input", "Grass", "data.asc"))
# # crop to tropical AOI
# dt_trop <- crop(dt, ext)
# # A few points are -0.09 with no apparent reason. 
# dt_trop[dt_trop<0] <- NA 
# names(dt_trop) <- "Grass"
# 
# summary(values(dt_trop))

## Create a brick for convenience. 
rasterlist_gaez <- list.files(path = targetdir, 
                              pattern = "", 
                              full.names = TRUE) %>% as.list()
gaez_all <- brick(rasterlist_gaez)

writeRaster(gaez_all, here("temp_data", "GAEZ", "South_America", "v4", "AES_index_value", "Rain-fed", "high_input_all.tif"), 
            overwrite = TRUE)

rm(rasterlist_gaez, gaez_all)


### PREPARE ATTAINABLE YIELDS
# dir.create(here("temp_data", "GAEZ", "South_America", "AES_index_value", "Rain-fed", "High-input"), recursive = TRUE) 
# dir.create(here("temp_data", "GAEZ", "South_America", "AES_index_value_current", "Rain-fed", "High-input"), recursive = TRUE) 

datadir <- here("input_data", "GAEZ", "v4", "AEAY_out_density", "Rain-fed", "High-input")
targetdir <- here("temp_data", "GAEZ", "South_America", "v4", "AEAY_out_density", "Rain-fed", "High-input")
tmpdir <- here("temp_data", "tmp")

if (!dir.exists(tmpdir)) dir.create(tmpdir, recursive = TRUE)
if (dir.exists(targetdir)) {
  file.remove(list.files(path = targetdir,
                         pattern = ".tif", full.names = TRUE))
} else dir.create(targetdir, recursive = TRUE)



files <- list.files(path = datadir, pattern = ".tif")
crops <- unlist(strsplit(files, split = ".tif"))
for(crop in crops){
  # note the difference here, we match *yl*Hr and not sxHr
  crops[grepl(crop, crops)] <- mapmat[grepl(crop, paste0("ylHr_", mapmat[,"abrev"])),"Names"]
}

## Import most crops
for (j in 1:length(files)) {
  #if (any(crops[j] == cropsToAggregate)) next
  print(files[j])
  unzip(zipfile = here(datadir, files[j]), exdir = tmpdir)
  dt <- raster(paste0(tmpdir, "/data.asc"))
  
  # crop to tropical AOI
  dt_trop <- crop(dt, ext)
  
  # A few points are -0.09 with no apparent reason. 
  dt_trop[dt_trop<0] <- NA 
  
  names(dt_trop) <- crops[j]
  writeRaster(dt_trop,
              filename = here(targetdir, paste0(crops[j], ".tif")),
              overwrite = TRUE)
  
}
rm(dt, dt_trop)


## Create a brick for convenience. 
rasterlist_gaez <- list.files(path = targetdir, 
                              pattern = "", 
                              full.names = TRUE) %>% as.list()
gaez_all <- brick(rasterlist_gaez)

writeRaster(gaez_all, here("temp_data", "GAEZ", "South_America", "v4", "AEAY_out_density", "Rain-fed", "high_input_all.tif"), 
            overwrite = TRUE)


rm(rasterlist_gaez, gaez_all)


#### PREPARE GLASS GLC DATA #### 
# We reproduce the data preparation workflow applied to years 1983-2015, even if only 2001-2015 is necessary, it is safer. 
# Read data in
rasterlist <- list.files(path = "input_data/GLASS-GLC", 
                         pattern = paste0("GLASS-GLC_7classes_"), 
                         full.names = TRUE) %>% as.list()
parcels_brick <- brick(rasterlist)

# crop to ***SOUTH AMERICA*** AOI 
southam_aoi <- crop(parcels_brick, ext)

# Write
writeRaster(southam_aoi, here("temp_data", "processed_glass-glc", "southam_aoi", "brick_southam_aoi.tif"), 
            overwrite = TRUE)


southam_aoi <- brick( here("temp_data", "processed_glass-glc", "southam_aoi", "brick_southam_aoi.tif"))

# Create annual layers of forest loss defined as: class is not 20 in a given year while it was 20 in *all the previous year*.
# This restricts annual loss to that occurring for the first time

# construct previous: a collection of annual layers, each giving the mean of the class value in the previous years. 
previous_years <- seq(1982, 2014, 1) 
for(t in 1:length(previous_years)){ # goes only up to 2014, as we don't need the average up to 2015.
  calc(southam_aoi[[1:t]], fun = mean, 
       filename = here("temp_data", "processed_glass-glc", "southam_aoi", paste0("past_mean_lu_",previous_years[t], ".tif")), 
       datatype = "FLT4S", # necessary so that a 19.9 mean is not counted as a 20 (i.e. so far undisturbed forest pixel)
       overwrite = TRUE)
}

# this is a raster of 33 layers, giving the mean of GLC class value in 1982, 1982-83, 1982-84, ..., 1982-2014. 
rasterlist <- list.files(path = here("temp_data", "processed_glass-glc", "southam_aoi"), 
                         pattern = "past_mean_lu_", 
                         full.names = TRUE) %>% as.list()
previous <- brick(rasterlist)

#unique(values(previous))

make_first_loss <- function(previous, current){if_else(condition = (previous == 20 & current != 20), 
                                                       true = 1, false = 0)}

years <- seq(1982, 2015, 1) 
for(t in 2:length(years)){ # starts from 1983 as we need t-1 and thus t starts from 2
  overlay(previous[[t-1]], southam_aoi[[t]], fun = make_first_loss, 
          filename = here("temp_data", "processed_glass-glc", "southam_aoi", paste0("first_loss_",years[t], ".tif")), 
          datatype = "INT1U", 
          overwrite = TRUE) 
}




rasterlist <- list.files(path = here("temp_data", "processed_glass-glc", "southam_aoi"), 
                         pattern = "first_loss_", 
                         full.names = TRUE) %>% as.list()
first_loss <- brick(rasterlist)

aggregate(first_loss, 
          fact = 2, 
          fun = sum, 
          expand = FALSE, 
          na.rm = FALSE, 
          filename = here("temp_data", "processed_glass-glc", "southam_aoi", "aggr_first_loss.tif"),
          datatype = "INT1U",
          overwrite = TRUE)


aggr_first_loss <- brick(here("temp_data", "processed_glass-glc", "southam_aoi", "aggr_first_loss.tif"))

gaez <- raster(here("temp_data", "GAEZ", "South_America", "v4", "AES_index_value", "Rain-fed", "High-input", "Banana.tif"))

beginCluster() # this uses by default detectCores() - 1

resample(x = aggr_first_loss, y = gaez,
         method = "ngb",
         filename = here("temp_data", "processed_glass-glc", "southam_aoi", "resampled_first_loss.tif"),
         datatype = "INT1U",
         overwrite = TRUE )

endCluster()


## Create the mask layer
# Create a layer that has values either : NA if first_loss always 0 across all years, 1 otherwise
res_first_loss <- brick(here("temp_data", "processed_glass-glc", "southam_aoi", "resampled_first_loss.tif"))
# not using if_else here to allow NA as an output... 
always_zero <- function(y){
  if(sum(y, na.rm = TRUE) == 0){d <- NA}else{d <- 1}
  return(d)}

mask_path <- here("temp_data", "processed_glass-glc", "southam_aoi", "always_zero_mask.tif")

overlay(x = res_first_loss, 
        fun = always_zero, 
        filename = mask_path,
        na.rm = TRUE, # but there is no NA anyway
        datatype = "INT2U", # INT2U to allow have NAs
        overwrite = TRUE)  

mask <- raster(mask_path)
# plot(mask)

# then use it to mask the brick 

mask(x = res_first_loss, 
     mask = mask, 
     filename = here("temp_data", "processed_glass-glc", "southam_aoi", "masked_first_loss.tif"), 
     datatype = "INT2U", 
     overwrite = TRUE)

first_loss <- brick(here("temp_data", "processed_glass-glc", "southam_aoi", "masked_first_loss.tif"))

cell_area <- area(first_loss)

make_area <- function(values, areas){
  # *0.25 to transform 0:4 scale into proportion of cell. *100 to convert km2 (returned by raster::area) to hectares, to match phtfl 
  return(values*0.25*areas*100) 
}

overlay(first_loss, cell_area, 
        fun = make_area, 
        filename = here("temp_data", "processed_glass-glc", "southam_aoi", "ha_first_loss.tif"), 
        overwrite = TRUE)


#### MERGE AESI #### 

### NEW FOLDERS USED IN THIS SCRIPT 
dir.create(here("temp_data", "merged_datasets", "southam_aoi"), recursive = TRUE)

### This script's target dir
targetdir <- here("temp_data", "merged_datasets", "southam_aoi")


### READ IN AND RENAME GAEZ DATA ### 

## SUITABILITY INDICES 
gaez_dir <- here("temp_data", "GAEZ", "South_America", "v4", "AES_index_value", "Rain-fed")
gaez <- brick(here(gaez_dir, "high_input_all.tif"))

# Rename layers (will be lost when writing the masked_gaez in the current code, so useless here and we rename later)
gaez_crops <- list.files(path = here(gaez_dir, "High-input"), 
                         pattern = "", 
                         full.names = FALSE)
gaez_crops <- gsub(pattern = ".tif", replacement = "", x = gaez_crops)
names(gaez) <- gaez_crops



### 1. MERGE GLASS-GLC AND SUITABILITY INDICES

## Read in GLASS-GLC data 
first_loss <- brick(here("temp_data", "processed_glass-glc", "southam_aoi", "ha_first_loss.tif"))

# Rename layers 
glc_sbqt_years <- seq(1983, 2015, 1)
names(first_loss) <- paste0("first_loss.",glc_sbqt_years)

glass <- first_loss


### MASK GAEZ TO REMOVE ALWAYS ZERO PIXELS AND LIGHTEN THE DATA FRAMES ### 

mask <- raster(here("temp_data", "processed_glass-glc", "southam_aoi", "always_zero_mask.tif"))

mask(x = gaez, 
     mask = mask, 
     filename = here("temp_data", "GAEZ", "South_America", "v4", "AES_index_value", "Rain-fed", "glass_masked_high_input_all.tif"), 
     overwrite = TRUE)

gaez_m <- brick(here("temp_data", "GAEZ", "South_America", "v4", "AES_index_value", "Rain-fed", "glass_masked_high_input_all.tif"))
# Rename layers (important, as writing the masked gaez lost the layer names)
names(gaez_m) <- gaez_crops

# (note that masking changes the summary values of gaez)


# Stack together the annual layers of GLASS-GLC data and GAEZ crop cross sections 
glass_gaez <- stack(glass, gaez_m)
names(glass_gaez)


### RASTER TO DATAFRAME ### 

# na.rm = TRUE is key here, as it removes previously masked pixels (NA) and ensures the output is not too large (memory intensive)
# We also set long to false because we reshape with a proper function for more control
wide_df <- raster::as.data.frame(glass_gaez, na.rm = TRUE, xy = TRUE, centroids = TRUE, long = FALSE) # ///!!!\\\ ~700s (rather 2-3 hours last time) 

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
varying_vars <- list(names(glass_gaez)[grep("first_loss.", names(glass_gaez), fixed = TRUE)])


# reshape to long.
long_df <- stats::reshape(wide_df,
                          varying = varying_vars,
                          v.names = c("first_loss"),
                          sep = ".",
                          timevar = "year",
                          idvar = "grid_id", # don't put "lon" and "lat" in there, otherwise memory issue (see https://r.789695.n4.nabble.com/reshape-makes-R-run-out-of-memory-PR-14121-td955889.html)
                          ids = "grid_id", # lonlat is our cross-sectional identifier.
                          direction = "long",
                          new.row.names = NULL)#seq(from = 1, to = nrow(ibs_msk_df)*length(years), by = 1)
rm(wide_df)
names(long_df)
# replace the indices from the raster::as.data.frame with actual years.

long_df <- mutate(long_df, year = glc_sbqt_years[year])

long_df <- dplyr::arrange(long_df, grid_id, year)


saveRDS(long_df, here(targetdir, "glass_aesi_long.Rdata"))

rm(long_df, varying_vars, glass_gaez, gaez_m, mask, glass, glc_sbqt_years, first_loss)



#### MERGE AEAY #### 

### NEW FOLDERS USED IN THIS SCRIPT 
dir.create(here("temp_data", "merged_datasets", "southam_aoi"), recursive = TRUE)

### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")


### This script's target dir
targetdir <- here("temp_data", "merged_datasets", "southam_aoi")


### READ IN AND RENAME GAEZ DATA ### 

## SUITABILITY INDICES 
gaez_dir <- here("temp_data", "GAEZ", "South_America", "v4", "AEAY_out_density", "Rain-fed")
gaez <- brick(here(gaez_dir, "high_input_all.tif"))

# Rename layers (will be lost when writing the masked_gaez in the current code, so useless here and we rename later)
gaez_crops <- list.files(path = here(gaez_dir, "High-input"), 
                         pattern = "", 
                         full.names = FALSE)
gaez_crops <- gsub(pattern = ".tif", replacement = "", x = gaez_crops)
names(gaez) <- gaez_crops




### 1. MERGE GLASS-GLC AND SUITABILITY INDICES

## Read in GLASS-GLC data 
first_loss <- brick(here("temp_data", "processed_glass-glc", "southam_aoi", "ha_first_loss.tif"))

# Rename layers 
glc_sbqt_years <- seq(1983, 2015, 1)
names(first_loss) <- paste0("first_loss.",glc_sbqt_years)


glass <- first_loss


### MASK GAEZ TO REMOVE ALWAYS ZERO PIXELS AND LIGHTEN THE DATA FRAMES ### 

mask <- raster(here("temp_data", "processed_glass-glc", "southam_aoi", "always_zero_mask.tif"))

mask(x = gaez, 
     mask = mask, 
     filename = here("temp_data", "GAEZ", "South_America", "v4", "AEAY_out_density", "Rain-fed", "glass_masked_high_input_all.tif"), 
     overwrite = TRUE)

gaez_m <- brick(here("temp_data", "GAEZ", "South_America", "v4", "AEAY_out_density", "Rain-fed", "glass_masked_high_input_all.tif"))
# Rename layers (important, as writing the masked gaez lost the layer names)
names(gaez_m) <- gaez_crops

# (note that masking changes the summary values of gaez)


### STACK RASTERS TO MERGE ###

# # Mask all the layers with ocean mask from GAEZ (more masked pixels than phtfloss that has only some rectangles between continents masked. plot both to see this)
# # take one layer from gaez
# gaez_mask <- gaez[[1]]
# mask(x = glass, mask = gaez_mask, 
#      filename = here("temp_data", "processed_glass-glc", "southam_aoi", "glass_masked.tif"),
#      overwrite = TRUE)
# 
# glass <- brick(here("temp_data", "processed_glass-glc", "southam_aoi", "glass_masked.tif"))
# names(glass) <- glass_names

# Stack together the annual layers of GLASS-GLC data and GAEZ crop cross sections 
glass_gaez <- stack(glass, gaez_m)
names(glass_gaez)


### RASTER TO DATAFRAME ### 

# na.rm = TRUE is key here, as it removes previously masked pixels (NA) and ensures the output is not too large (memory intensive)
# We also set long to false because we reshape with a proper function for more control
wide_df <- raster::as.data.frame(glass_gaez, na.rm = TRUE, xy = TRUE, centroids = TRUE, long = FALSE) # ///!!!\\\ ~700s (rather 2-3 hours last time) 

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
varying_vars <- list(names(glass_gaez)[grep("first_loss.", names(glass_gaez), fixed = TRUE)])


# reshape to long.
long_df <- stats::reshape(wide_df,
                          varying = varying_vars,
                          v.names = c("first_loss"),
                          sep = ".",
                          timevar = "year",
                          idvar = "grid_id", # don't put "lon" and "lat" in there, otherwise memory issue (see https://r.789695.n4.nabble.com/reshape-makes-R-run-out-of-memory-PR-14121-td955889.html)
                          ids = "grid_id", # lonlat is our cross-sectional identifier.
                          direction = "long",
                          new.row.names = NULL)#seq(from = 1, to = nrow(ibs_msk_df)*length(years), by = 1)
rm(wide_df)
names(long_df)
# replace the indices from the raster::as.data.frame with actual years.

long_df <- mutate(long_df, year = glc_sbqt_years[year])

long_df <- dplyr::arrange(long_df, grid_id, year)


saveRDS(long_df, here(targetdir, "glass_aeay_long.Rdata"))

rm(long_df, varying_vars, glass_gaez, gaez_m, mask, glass, glc_sbqt_years, first_loss)


# glass <- readRDS(here(targetdir, "glass_aeay_long.Rdata"))


#### ADD VARIABLES AESI #### 

# Define GAEZ AESI variables 
gaez_crops <- list.files(path = here("temp_data", "GAEZ", "South_America", "v4", "AES_index_value", "Rain-fed", "High-input"), 
                         pattern = "", 
                         full.names = FALSE)
gaez_crops <- gsub(pattern = ".tif", replacement = "", x = gaez_crops)

# Contrary to what is done in the main add_variables.R script, here we execute code only for first_loss.
origindir <- here("temp_data", "merged_datasets", "southam_aoi")
name <- "glass_aesi_long"

### COUNTRIES ### 

countries <- st_read(here("input_data", "Global_LSIB_Polygons_Detailed"))

# Data need to be projected 
crs_southam <- 31970

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
                       countries$COUNTRY_NA=="Venezuela"  ,] %>% st_transform(crs_southam)
southam <- southam[,c("COUNTRY_NA", "geometry")]


path <- paste0(here(origindir, name), ".Rdata")
df <- readRDS(path)

# Remove gaez variables
df <- dplyr::select(df,-all_of(gaez_crops))

# Use cross section only
df_cs <- df[!duplicated(df$grid_id),]

rm(df)

# Spatial
df_cs <- st_as_sf(df_cs, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
df_cs <- st_transform(df_cs, crs= crs_southam)

# This is much much faster (like 5 minutes vs. 6h).
# df_cs <- st_join(x = countries[,c("OBJECTID", "COUNTRY_NA")],
#                  y = df_cs,
#                  join = st_contains,
#                  prepared = TRUE,
#                  left = FALSE)# performs inner join so returns only records that spatially match.
#

# However, we use st_nearest_feature so that all points match a country
df_cs <- st_join(x = df_cs,
                 y = southam,
                 join = st_nearest_feature,
                 left = TRUE)

# names(df_cs)[names(df_cs) == "OBJECTID"] <- "country_id"
names(df_cs)[names(df_cs) == "COUNTRY_NA"] <- "country_name"

df_cs <- st_drop_geometry(df_cs)

# Keep only new variable and id
df_cs <- df_cs[,c("grid_id", "country_name")]

df_cs$country_name[df_cs$country_name=="Isla Brasilera (disp)"] <- "Brazil"
df_cs$country_name[df_cs$country_name=="Bolivia"] <- "Bolivia (Plurinational State of)"
df_cs$country_name[df_cs$country_name=="Venezuela"] <- "Venezuela (Bolivarian Republic of)"
df_cs$country_name[df_cs$country_name=="French Guiana (Fr)"] <- "France"
df_cs$country_name[df_cs$country_name=="Trinidad & Tobago"] <- "Trinidad and Tobago"

# # We save the cross section, not the panel, as it is not necessary 
# saveRDS(df_cs, path)
saveRDS(df_cs, paste0(here(origindir, name), "_country_nf.Rdata"))
rm(df_cs)


### STANDARDIZE AND AGGREGATE ### 

path <- paste0(here(origindir, name), ".Rdata")
df <- readRDS(path)
# Use cross section only
df_cs <- df[!duplicated(df$grid_id),]

rm(df)

## Aggregate suitability indexes to match price data
df_cs <- df_cs %>% rowwise() %>% mutate(Sugar = max(c(Sugarbeet, Sugarcane)),# Especially necessary to match the international price of sugar
                                        Fodder = max(c(Alfalfa, Napiergrass)), # To adjust once we have Grass as a GAEZ crop  
                                        Rice = max(c(Drylandrice, Wetlandrice)),
                                        Sorghum2 = max(c(Sorghum, Sorghumbiomass)),
                                        bioenergy_crops = max(c(Jatropha, Miscanthus, Reedcanarygrass, Switchgrass)),
                                        cereal_crops = max(c(Buckwheat, Foxtailmillet, Pearlmillet, Rye)),
                                        pulses_crops = max(c(Chickpea, Cowpea, Drypea, Gram, Phaseolousbean, Pigeonpea)),
                                        roots_crops = max(c(Cassava, Sweetpotato, Whitepotato, Yam)),
                                        # oil_crops = max(c(Groundnut, Jatropha, Olive, Rapeseed, Sunflower)), # we have prices for all of them
                                        vegetables_crops = max(c(Cabbage, Carrot, Onion, Tomato)),
                                        # fruits_crops = max(c(Banana, Citrus, Coconut)), # For coconut we have price only for the oil but it does not matter because it's the AESI workstream anyway 
                                        industrial_crops = max(c(Flax)) # Rubber, directly responsible for deforestation, are not mixed with these.
                                        # narcotics_crops = max(c(Cocoa, Coffee, Tea, Tobacco)), # We have prices for all of them
                                        
) %>% as.data.frame()

# sugar crops and oil crops could alternatively be categorized as bioenergy feedstock, and Miscanthus etc. as fodder crops (according to Wikipedia).

## Standardize crops that can be matched with a price and crop groups
indivcrops_to_std <- c("Banana", "Barley", "Citrus", "Cocoa", "Coconut", "Coffee", "Cotton", "Fodder", 
                       "Groundnut", "Maizegrain", "Oat", "Oilpalm", "Olive", "Rapeseed", "Rice", "Rubber", 
                       "Sorghum2", "Soybean", "Sugar", "Sunflower", "Tea", "Tobacco", "Wheat")

# To understand these lines, see https://dplyr.tidyverse.org/articles/rowwise.html#row-wise-summary-functions
df_cs <- dplyr::mutate(df_cs, si_sum = rowSums(across(.cols = (contains("_crops") | any_of(indivcrops_to_std)))))

df_cs <- dplyr::mutate(df_cs, across(.cols = (contains("_crops") | any_of(indivcrops_to_std)),
                                     .fns = ~./si_sum, 
                                     .names = paste0("{.col}", "_std")))

# Select only newly constructed variables, and id
var_names <- names(df_cs)[grepl(pattern = "_std", x = names(df_cs))] # this is to add if we want non stded grouped crops too: | grepl(pattern = "_crops", x = names(df_cs))
# and this is if we want the new variables in non std format too | names(df_cs) %in% c("Fodder", "Rice", "Sugar")
df_cs <- df_cs[,c("grid_id", var_names)]

saveRDS(df_cs, paste0(here(origindir, name), "_stdsi.Rdata"))  
rm(df_cs, path)


### MERGE AESI DATASETS WITH ADDED VARIABLES ### 

# Base dataset (including outcome variable(s))
base_path <- paste0(here(origindir, name), ".Rdata")
df_base <- readRDS(base_path)

# Remove non-standardized suitability indexes
df_base <- dplyr::select(df_base,-all_of(gaez_crops))

# Country variable
country_path <- paste0(here(origindir, name), "_country_nf.Rdata")
df_country <- readRDS(country_path)

# Merge them and remove to save memory 
final <- left_join(df_base, df_country, by = "grid_id")
rm(df_base, df_country)

# Standardized and aggregated suitability indexes
stdsi_path <- paste0(here(origindir, name), "_stdsi.Rdata")
df_stdsi <- readRDS(stdsi_path)  

final <- left_join(final, df_stdsi, by = "grid_id")
rm(df_stdsi)

# Create country trends variable
final <- mutate(final, country_year = paste0(country_name, "_", year))

saveRDS(final, paste0(here(origindir, name), "_final.Rdata"))

rm(final)


#### ADD VARIABLES AEAY #### 

# Define GAEZ AEAY variables 
gaez_crops <- list.files(path = here("temp_data", "GAEZ", "South_America", "v4", "AEAY_out_density", "Rain-fed", "High-input"), 
                         pattern = "", 
                         full.names = FALSE)
gaez_crops <- gsub(pattern = ".tif", replacement = "", x = gaez_crops)


origindir <- here("temp_data", "merged_datasets", "southam_aoi")
name <- "glass_aeay_long"

### GROUP AEAY CROPS ### 
path <- paste0(here(origindir, name), ".Rdata")
df <- readRDS(path)
# Use cross section only
df_cs <- df[!duplicated(df$grid_id),]

rm(df)

## Aggregate suitability indexes to match price data
df_cs <- df_cs %>% rowwise() %>% mutate(Sugar = max(c(Sugarbeet, Sugarcane)),# Especially necessary to match the international price of sugar
                                        Fodder = max(c(Alfalfa, Napiergrass)),   
                                        Rice = max(c(Drylandrice, Wetlandrice)),
                                        Sorghum2 = max(c(Sorghum, Sorghumbiomass))
                                        # We surely wont need these in the AEAY workflow, as we need to interact these variables with a price and there is no price for those. 
                                        # bioenergy_crops = max(c(Jatropha, Miscanthus, Reedcanarygrass, Sorghumbiomass, Switchgrass)),
                                        # cereal_crops = max(c(Buckwheat, Foxtailmillet, Pearlmillet, Rye)),
                                        # pulses_crops = max(c(Chickpea, Cowpea, Drypea, Gram, Phaseolousbean, Pigeonpea)),
                                        # roots_crops = max(c(Cassava, Sweetpotato, Whitepotato, Yam)),
                                        # # oil_crops = max(c(Groundnut, Jatropha, Olive, Rapeseed, Sunflower)), # we have prices for all of them
                                        # vegetables_crops = max(c(Cabbage, Carrot, Onion, Tomato)),
                                        # # fruits_crops = max(c(Banana, Citrus, Coconut)), # For coconut we have price only for the oil but it does not matter because it's the AESI workstream anyway 
                                        # industrial_crops = max(c(Flax)) # Rubber, directly responsible for deforestation, are not mixed with these.
                                        # # narcotics_crops = max(c(Cocoa, Coffee, Tea, Tobacco)), # We have prices for all of them
) %>% as.data.frame()

# Select only newly constructed variables, and id
df_cs <- df_cs[,c("grid_id", "Fodder", "Rice", "Sugar", "Sorghum2")]

# # merge back to panel - NOPE, not anymore
# df <- left_join(df, df_cs, by = "grid_id")

# Keep only new variables and id
# df_cs <- df_cs[,(names(df_cs)=="grid_id" | grepl(pattern = "_std", x = names(df_cs)))]

saveRDS(df_cs, paste0(here(origindir, name), "_groupedcrops.Rdata"))  
rm(df_cs, path)


### MERGE AEAY DATASETS WITH ADDED VARIABLES

# Base dataset (including outcome variable(s))
base_path <- paste0(here(origindir, name), ".Rdata")
df_base <- readRDS(base_path)

# Country variable
# do that because we did not run the spatial join with countries for aeay dataset
country_path <- paste0(here(origindir, "glass_aesi_long"), "_country_nf.Rdata")
df_country <- readRDS(country_path)

# Merge them and remove to save memory 
final <- left_join(df_base, df_country, by = "grid_id")
rm(df_base, df_country)

# Grouped crops
grouped_path <- paste0(here(origindir, name), "_groupedcrops.Rdata")
df_grouped <- readRDS(grouped_path)  

final <- left_join(final, df_grouped, by = "grid_id")
rm(df_grouped)

# Create country trends variable
final <- mutate(final, country_year = paste0(country_name, "_", year))

saveRDS(final, paste0(here(origindir, name), "_final.Rdata"))

rm(final)