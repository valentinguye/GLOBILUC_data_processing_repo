# In this script, we merge in 3 steps: 
# 1. First, we merge GLASS-GLC data with GAEZ, 
# 2. Second, we merge the first loss variable computed for 1983 to 2020 with GAEZ 
# 3. Third, we merge phtfl loss data with GAEZ
# We don't merge everything in a single data frame because it is memory intensive.  
# Also because we mask gaez with masks based on outcome variables being "always zero" and thus the 
# masks are different for first loss 1983-2015 (three variables of GLASS-GLC), first loss 1983-2020, and for phtfloss. 


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
dir.create(here("temp_data", "merged_datasets", "tropical_aoi"), recursive = TRUE)

### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")


### This script's target dir
targetdir <- here("temp_data", "merged_datasets", "tropical_aoi")


# About layers renaming
# The index given to each layer in the bricks we read here, is representative of the increasing alphabetical order of the names of 
# the individual layers originally read in in the prepare_ scripts (this is ensured by the use of list.files() to get these names.)
# Thus we can rename with name vectors that reproduce this order.  

### READ IN AND RENAME GAEZ DATA ### 

## SUITABILITY INDICES 
gaez_dir <- here("temp_data", "GAEZ", "v4", "AES_index_value", "Rain-fed")
gaez <- brick(here(gaez_dir, "high_input_all.tif"))

# Rename layers (will be lost when writing the masked_gaez in the current code, so useless here and we rename later)
gaez_crops <- list.files(path = here(gaez_dir, "High-input"), 
                         pattern = "", 
                         full.names = FALSE)
gaez_crops <- gsub(pattern = ".tif", replacement = "", x = gaez_crops)
names(gaez) <- gaez_crops

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 



#### 1. MERGE GLASS-GLC AND GAEZ #### 

## Read in GLASS-GLC data 
first_loss <- brick(here("temp_data", "processed_glass-glc", "tropical_aoi", "ha_first_loss.tif"))
nd_first_loss <- brick(here("temp_data", "processed_glass-glc", "tropical_aoi", "ha_nd_first_loss.tif"))
sbqt_direct_lu <- brick(here("temp_data", "processed_glass-glc", "tropical_aoi", "masked_sbqt_direct_lu.tif"))
sbqt_mode_lu <- brick(here("temp_data", "processed_glass-glc", "tropical_aoi", "masked_sbqt_mode_lu.tif"))
# Rename layers 
glc_sbqt_years <- seq(1983, 2015, 1)
names(first_loss) <- paste0("first_loss.",glc_sbqt_years)
names(nd_first_loss) <- paste0("nd_first_loss.", seq(1983, 2014, 1))
names(sbqt_direct_lu) <- paste0("sbqt_direct_lu.",glc_sbqt_years)
names(sbqt_mode_lu) <- paste0("sbqt_mode_lu.",glc_sbqt_years)

glass <- stack(first_loss, nd_first_loss, sbqt_direct_lu, sbqt_mode_lu)  
# # save names as they will be lost when writing the masked data 
# glass_names <- names(glass)

# banana1 <- gaez_all[[2]] 
# banana2 <- gaez[[2]] 
# 
# all.equal(values(banana1), values(banana2))

# 
# banana <- gaez_all[[2]]
# pl <- phtfloss[[1]]
# 
# plot(banana)
# vbanana <- values(banana)
# 
# sum(is.na(vbanana))
# 
# sum(is.na(values(pl)))
# plot(pl)


### MASK GAEZ TO REMOVE ALWAYS ZERO PIXELS AND LIGHTEN THE DATA FRAMES ### 

mask <- raster(here("temp_data", "processed_glass-glc", "tropical_aoi", "always_zero_mask.tif"))

mask(x = gaez, 
     mask = mask, 
     filename = here("temp_data", "GAEZ", "v4", "AES_index_value", "Rain-fed", "glass_masked_high_input_all.tif"), 
     overwrite = TRUE)

gaez_m <- brick(here("temp_data", "GAEZ", "v4", "AES_index_value", "Rain-fed", "glass_masked_high_input_all.tif"))
# Rename layers (important, as writing the masked gaez lost the layer names)
names(gaez_m) <- gaez_crops

# (note that masking changes the summary values of gaez)


### STACK RASTERS TO MERGE ###

# # Mask all the layers with ocean mask from GAEZ (more masked pixels than phtfloss that has only some rectangles between continents masked. plot both to see this)
# # take one layer from gaez
# gaez_mask <- gaez[[1]]
# mask(x = glass, mask = gaez_mask, 
#      filename = here("temp_data", "processed_glass-glc", "tropical_aoi", "glass_masked.tif"),
#      overwrite = TRUE)
# 
# glass <- brick(here("temp_data", "processed_glass-glc", "tropical_aoi", "glass_masked.tif"))
# names(glass) <- glass_names

# add a layer for year 2015 for nd_first_loss (necessary for reshape to work)
# set its value to -1, not NA, otherwise all obs. get removed in as.data.frame below
nd_first_loss.2015 <- raster(glass, layer = 0) 
values(nd_first_loss.2015) <- -1
names(nd_first_loss.2015) <- "nd_first_loss.2015"

# Stack together the annual layers of GLASS-GLC data and GAEZ crop cross sections 
glass_gaez <- stack(glass, nd_first_loss.2015, gaez_m) # 
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
varying_vars <- list(paste0("first_loss.", seq(1983, 2015, 1)),
                     paste0("nd_first_loss.", seq(1983, 2015, 1)),
                     paste0("sbqt_direct_lu.", seq(1983, 2015, 1)),
                     paste0("sbqt_mode_lu.", seq(1983, 2015, 1)))
                    
# reshape to long.
long_df <- stats::reshape(wide_df,
                       varying = varying_vars,
                       v.names = c("first_loss", "nd_first_loss", "sbqt_direct_lu", "sbqt_mode_lu"),
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

# replace the -1 values in the vector of nd_first_loss in 2015 by NAs
long_df[long_df$year == 2015, "nd_first_loss"] <- NA

# rearrange
long_df <- dplyr::arrange(long_df, grid_id, year)

saveRDS(long_df, here(targetdir, "glass_aesi_long.Rdata"))

rm(long_df, varying_vars, glass_gaez, gaez_m, mask, glass, glc_sbqt_years, first_loss, sbqt_direct_lu, sbqt_mode_lu)


# glass <- readRDS(here(targetdir, "glass_aesi_long.Rdata"))


#### 2. MERGE 1983-2020 FIRST LOSS AND GAEZ #### 

## Read in FIRST LOSS 1983 2020 
firstloss8320 <- brick(here("temp_data", "processed_glass-glc", "tropical_aoi", "ha_firstloss8320.tif")) 
# Rename layers 
names(firstloss8320) <- paste0("firstloss_glassgfc.",seq(1983, 2020, 1))


### MASK GAEZ TO REMOVE ALWAYS ZERO PIXELS AND LIGHTEN THE DATA FRAMES ### 

mask <- raster(here("temp_data", "processed_glass-glc", "tropical_aoi", "always_zero_mask8320.tif"))

mask(x = gaez, 
     mask = mask, 
     filename = here("temp_data", "GAEZ", "v4", "AES_index_value", "Rain-fed", "firstloss8320_masked_high_input_all.tif"), 
     overwrite = TRUE)

gaez_m <- brick(here("temp_data", "GAEZ", "v4", "AES_index_value", "Rain-fed", "firstloss8320_masked_high_input_all.tif"))
# Rename layers (important, as writing the masked gaez lost the layer names)
names(gaez_m) <- gaez_crops

# (note that masking changes the summary values of gaez)


### STACK RASTERS TO MERGE ###

# Stack together the annual layers of firstloss8320-GLC data and GAEZ crop cross sections 
firstloss_gaez <- stack(firstloss8320, gaez_m)
names(firstloss_gaez)


### RASTER TO DATAFRAME ### 

# na.rm = TRUE is key here, as it removes previously masked pixels (NA) and ensures the output is not too large (memory intensive)
# We also set long to false because we reshape with a proper function for more control
wide_df <- raster::as.data.frame(firstloss_gaez, na.rm = TRUE, xy = TRUE, centroids = TRUE, long = FALSE) # 

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
varying_vars <- names(firstloss_gaez)[grep(".", names(firstloss_gaez), fixed = TRUE)]

# reshape to long.
long_df <- stats::reshape(wide_df,
                          varying = varying_vars,
                          v.names = "firstloss_glassgfc",
                          sep = ".",
                          timevar = "year",
                          idvar = "grid_id", # don't put "lon" and "lat" in there, otherwise memory issue (see https://r.789695.n4.nabble.com/reshape-makes-R-run-out-of-memory-PR-14121-td955889.html)
                          ids = "grid_id", # lonlat is our cross-sectional identifier.
                          direction = "long",
                          new.row.names = NULL)#seq(from = 1, to = nrow(ibs_msk_df)*length(years), by = 1)
rm(wide_df)
names(long_df)
# replace the indices from the raster::as.data.frame with actual years.

years <- seq(1983, 2020, 1)
long_df <- mutate(long_df, year = years[year])

long_df <- dplyr::arrange(long_df, grid_id, year)


saveRDS(long_df, here(targetdir, "firstloss8320_aesi_long.Rdata"))

rm(years, long_df, varying_vars, firstloss_gaez, gaez_m, mask, firstloss8320)



#### 3. MERGE PHTF LOSS AND GAEZ #### 

## Read in PHTF LOSS 
phtfloss <- brick(here("temp_data", "processed_phtfloss", "tropical_aoi", "masked_phtfloss.tif")) 
# Rename layers 
names(phtfloss) <- paste0("phtf_loss.",seq(2002, 2020, 1))


### MASK GAEZ TO REMOVE ALWAYS ZERO PIXELS AND LIGHTEN THE DATA FRAMES ### 

mask <- raster(here("temp_data", "processed_phtfloss", "tropical_aoi", "always_zero_mask_phtfloss.tif"))

mask(x = gaez, 
     mask = mask, 
     filename = here("temp_data", "GAEZ", "v4", "AES_index_value", "Rain-fed", "phtfloss_masked_high_input_all.tif"), 
     overwrite = TRUE)

gaez_m <- brick(here("temp_data", "GAEZ", "v4", "AES_index_value", "Rain-fed", "phtfloss_masked_high_input_all.tif"))
# Rename layers (important, as writing the masked gaez lost the layer names)
names(gaez_m) <- gaez_crops

# (note that masking changes the summary values of gaez)


### STACK RASTERS TO MERGE ###

# Stack together the annual layers of phtfloss-GLC data and GAEZ crop cross sections 
phtfloss_gaez <- stack(phtfloss, gaez_m)
names(phtfloss_gaez)


### RASTER TO DATAFRAME ### 

# na.rm = TRUE is key here, as it removes previously masked pixels (NA) and ensures the output is not too large (memory intensive)
# We also set long to false because we reshape with a proper function for more control
wide_df <- raster::as.data.frame(phtfloss_gaez, na.rm = TRUE, xy = TRUE, centroids = TRUE, long = FALSE) # ~700s. 

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
varying_vars <- names(phtfloss_gaez)[grep(".", names(phtfloss_gaez), fixed = TRUE)]

# reshape to long.
long_df <- stats::reshape(wide_df,
                          varying = varying_vars,
                          v.names = "phtf_loss",
                          sep = ".",
                          timevar = "year",
                          idvar = "grid_id", # don't put "lon" and "lat" in there, otherwise memory issue (see https://r.789695.n4.nabble.com/reshape-makes-R-run-out-of-memory-PR-14121-td955889.html)
                          ids = "grid_id", # lonlat is our cross-sectional identifier.
                          direction = "long",
                          new.row.names = NULL)#seq(from = 1, to = nrow(ibs_msk_df)*length(years), by = 1)
rm(wide_df)
names(long_df)
# replace the indices from the raster::as.data.frame with actual years.

years <- seq(2002, 2020, 1)
long_df <- mutate(long_df, year = years[year])

long_df <- dplyr::arrange(long_df, grid_id, year)


saveRDS(long_df, here(targetdir, "phtfloss_aesi_long.Rdata"))

rm(long_df, varying_vars, phtfloss_gaez, gaez_m, mask, phtfloss)


#### 4. MERGE DRIVER LOSS AND GAEZ #### 

## Read in DRIVER LOSS 
driverloss <- brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", "masked_lossdrivers.tif")) 
# Rename layers 
names(driverloss) <- paste0("driven_loss.",seq(2001, 2019, 1)) # note the difference with the names of phtfloss (not the same years)


driverloss_commodity <- brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", "masked_lossdrivers_commodity.tif")) 
# Rename layers 
names(driverloss_commodity) <- paste0("driven_loss_commodity.",seq(2001, 2019, 1)) # note the difference with the names of phtfloss (not the same years)

### MASK GAEZ TO REMOVE ALWAYS ZERO PIXELS AND LIGHTEN THE DATA FRAMES ### 

mask <- raster(here("temp_data", "processed_lossdrivers", "tropical_aoi", "always_zero_mask_lossdrivers.tif"))

mask(x = gaez, 
     mask = mask,
     maskvalue = 0, # necessary here, because the there is no NA in the mask, only 0 and 1 (see the prepare_loss_drivers.R script)
     updatevalue = NA, 
     filename = here("temp_data", "GAEZ", "v4", "AES_index_value", "Rain-fed", "driverloss_masked_high_input_all.tif"), 
     overwrite = TRUE)

gaez_m <- brick(here("temp_data", "GAEZ", "v4", "AES_index_value", "Rain-fed", "driverloss_masked_high_input_all.tif"))
# Rename layers (important, as writing the masked gaez lost the layer names)
names(gaez_m) <- gaez_crops

# (note that masking changes the summary values of gaez)

### ADD THE 2000 FOREST COVER ### 
# (it is already masked with always zero driven loss, in prepare_fc2000)
fc2k <- raster(here("temp_data", "processed_lossdrivers", "tropical_aoi", "masked_fc_2000.tif"))
names(fc2k) <- "fc_2000"

### STACK RASTERS TO MERGE ###

# Stack together the annual layers of driverloss-GLC data and GAEZ crop cross sections 
driverloss_gaez <- stack(driverloss, driverloss_commodity, gaez_m, fc2k)
names(driverloss_gaez)


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
varying_vars <- list(names(driverloss), 
                     names(driverloss_commodity))
#varying_vars <- names(driverloss_gaez)[grep(".", names(driverloss_gaez), fixed = TRUE)]

# reshape to long.
long_df <- stats::reshape(wide_df,
                          varying = varying_vars,
                          v.names = c("driven_loss","driven_loss_commodity"),
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

saveRDS(long_df, here(targetdir, "driverloss_aesi_long.Rdata"))

rm(long_df, varying_vars, driverloss_gaez, gaez_m, mask, driverloss, fc2k)

#### 5. MERGE ALL DRIVERS - AND GAEZ #### 

## Read in DRIVER LOSS 
driverloss_commodity <- brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", "masked_lossdrivers_commodity.tif")) 
# Rename layers 
names(driverloss_commodity) <- paste0("driven_loss_commodity.",seq(2001, 2019, 1)) # note the difference with the names of phtfloss (not the same years)

driverloss_shifting <- brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", "masked_lossdrivers_shifting.tif")) 
# Rename layers 
names(driverloss_shifting) <- paste0("driven_loss_shifting.",seq(2001, 2019, 1)) # note the difference with the names of phtfloss (not the same years)

driverloss_forestry <- brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", "masked_lossdrivers_forestry.tif")) 
# Rename layers 
names(driverloss_forestry) <- paste0("driven_loss_forestry.",seq(2001, 2019, 1)) # note the difference with the names of phtfloss (not the same years)

driverloss_fire <- brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", "masked_lossdrivers_fire.tif")) 
# Rename layers 
names(driverloss_fire) <- paste0("driven_loss_fire.",seq(2001, 2019, 1)) # note the difference with the names of phtfloss (not the same years)


### MASK GAEZ TO REMOVE ALWAYS ZERO PIXELS AND LIGHTEN THE DATA FRAMES ### 

mask <- raster(here("temp_data", "processed_lossdrivers", "tropical_aoi", "always_zero_mask_lossdrivers_any.tif"))

mask(x = gaez, 
     mask = mask,
     maskvalue = 0, # necessary here, because the there is no NA in the mask, only 0 and 1 (see the prepare_loss_drivers.R script)
     updatevalue = NA, 
     filename = here("temp_data", "GAEZ", "v4", "AES_index_value", "Rain-fed", "any_driverloss_masked_high_input_all.tif"), 
     overwrite = TRUE)

gaez_m <- brick(here("temp_data", "GAEZ", "v4", "AES_index_value", "Rain-fed", "any_driverloss_masked_high_input_all.tif"))
# Rename layers (important, as writing the masked gaez lost the layer names)
names(gaez_m) <- gaez_crops

# (note that masking changes the summary values of gaez)

### ADD THE 2000 FOREST COVER ### 
# (it is already masked with always zero driven loss, in prepare_fc2000)
fc2k <- raster(here("temp_data", "processed_lossdrivers", "tropical_aoi", "masked_fc_2000_any.tif"))
names(fc2k) <- "fc_2000"

### STACK RASTERS TO MERGE ###

# Stack together the annual layers of driverloss-GLC data and GAEZ crop cross sections 
driverloss_gaez <- stack(driverloss_commodity,
                         driverloss_shifting, 
                         driverloss_forestry, 
                         driverloss_fire, 
                         gaez_m, fc2k)
names(driverloss_gaez)


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


saveRDS(long_df, here(targetdir, "driverloss_all_aesi_long.Rdata"))

rm(long_df, varying_vars, driverloss_gaez, gaez_m, mask, driverloss, fc2k)


# Note: much more cells have no some deforestation in the driverloss data than in the first_loss. This is due to the starting resolution of each original data set. 



rm(gaez, gaez_crops, gaez_dir, targetdir)



# We rename layers in a specific way: 
# First we get names in the order from which they are read (but brick() seems to preserve the increasing alphabetical order of file names)
# Then we use a regex replacement method in case that layer_names was not saved in the order we think (increasing alphabetical)
# layer_names <- names(phtfloss)
# layer_names2 <- gsub(pattern = ".+?(resampled_phtf_loss.)", replacement = "phtf_loss_", layer_names)
# names(phtfloss) <- layer_names2



