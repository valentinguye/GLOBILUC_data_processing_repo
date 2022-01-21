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
dir.create(here("temp_data", "processed_fc2000", "tropical_aoi"), recursive = TRUE)

### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")

### TROPICAL AOI 
ext <- extent(c(-180, 179.9167, -30, 30))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


fc2k <- raster(here("input_data", "fc_2000_3km_10th.tif"))

# crop it to the tropical aoi
fc2k <- crop(fc2k, ext)

# Currently, the data are hectares of forest cover in 5x5km grid cells
# We aggregate this to the GAEZ resolution by adding up the hectares

gaez <- raster(here("temp_data", "GAEZ", "AES_index_value", "Rain-fed", "High-input", "Alfalfa.tif"))

# define output file name
aggr_output_name <- here("temp_data", "processed_fc2000", "tropical_aoi", "aggr_fc_2000.tif")

# aggregate it from the ~5km cells to ~10km
raster::aggregate(fc2k, fact = c(res(gaez)[1]/res(fc2k)[1], res(gaez)[2]/res(fc2k)[2]),
                  expand = FALSE,
                  fun = sum,
                  na.rm = TRUE, # NA values are only in the sea. Where there is no forest loss, like in a city in Brazil, the value is 0 (see with plot())
                  filename = aggr_output_name,
                  # datatype = "INT2U", # let the data be float, as we have decimals in the amount of hectares. 
                  overwrite = TRUE)

#### ALIGN to GAEZ ####
aggregated <- raster(aggr_output_name)

resampled_output_name <- here("temp_data", "processed_fc2000", "tropical_aoi", "resampled_fc_2000.tif")

resample(x = aggregated, 
         y = gaez, 
         method = "ngb", # we use ngb and not bilinear because the output values' summary better fits that of the aggregated layer 
         # and the bilinear interpolation arguably smoothes the reprojection more than necessary given that from and to are already very similar.  
         filename = resampled_output_name, 
         overwrite = TRUE)

### MASK WITH ALWAYS 0 PIXELS
resampled <- raster(resampled_output_name)

mask <- raster(here("temp_data", "processed_lossdrivers", "tropical_aoi", "always_zero_mask_lossdrivers.tif"))
# plot(mask)

mask(resampled, 
     mask = mask, 
     maskvalue = 0, # necessary here, because the always_zero function used converted to 0 and not to NA
     updatevalue = NA, 
     filename = here("temp_data", "processed_lossdrivers", "tropical_aoi", "masked_fc_2000.tif"), 
     overwrite = TRUE)

# AND WITH ANY DRIVER ALAYS-ZERO MASK 
any_mask <- raster(here("temp_data", "processed_lossdrivers", "tropical_aoi", "always_zero_mask_lossdrivers_any.tif"))
# plot(mask)

mask(resampled, 
     mask = any_mask, 
     maskvalue = 0, # necessary here, because the always_zero function used converted to 0 and not to NA
     updatevalue = NA, 
     filename = here("temp_data", "processed_lossdrivers", "tropical_aoi", "masked_fc_2000_any.tif"), 
     overwrite = TRUE)

