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
dir.create(here("temp_data", "processed_phtfloss", "tropical_aoi"), recursive = TRUE)

### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 



#### PREPARE PHTF LOSS #### 

### Brick layers ### 

# import annual layers of phtf loss (in hectares) as computed in GEE (and downloaded from Google Drive to input_data/)
rasterlist <- list.files(path = here("input_data", "phtfLossSumGlass_maxP"), 
                         pattern = "phtfLossSumGlass_maxP_", 
                         full.names = TRUE) %>% as.list()
phtfloss <- brick(rasterlist)

# test <- phtfloss[[1]]
# test
# plot(test)
# summary(values(test))

writeRaster(phtfloss, here("temp_data", "processed_phtfloss", "tropical_aoi", "phtf_loss.tif"), 
            overwrite = TRUE)


### Aggregate to the GAEZ resolution ### 

phtfloss <- brick(here("temp_data", "processed_phtfloss", "tropical_aoi", "phtf_loss.tif"))

# Currently, the data are hectares of phtf loss in 5x5km grid cells annually. 
# We aggregate this to 10km grid cells by adding up the hectares of phtf loss

# Read in a GAEZ grid, because we want to align to it. 
gaez <- raster(here("temp_data", "GAEZ", "AES_index_value", "Rain-fed", "High-input", "Alfalfa.tif"))

# define output file name
aggr_output_name <- here("temp_data", "processed_phtfloss", "tropical_aoi", "aggr_phtf_loss.tif")

# aggregate it from the ~5km cells to ~10km
raster::aggregate(phtfloss, fact = c(res(gaez)[1]/res(phtfloss)[1], res(gaez)[2]/res(phtfloss)[2]),
                  expand = FALSE,
                  fun = sum,
                  na.rm = TRUE, # NA values are only in the sea. Where there is no phtf loss, like in a city in Brazil, the value is 0 (see with plot())
                  filename = aggr_output_name,
                  # datatype = "INT2U", # let the data be float, as we have decimals in the amount of hectares. 
                  overwrite = TRUE)


### Align to GAEZ exactly ###
aggregated <- brick(aggr_output_name)

resample(x = aggregated, 
         y = gaez, 
         method = "ngb", # we use ngb and not bilinear because the output values' summary better fits that of the aggregated layer 
         # and the bilinear interpolation arguably smoothes the reprojection more than necessary given that from and to are already very similar.  
         filename = here("temp_data", "processed_phtfloss", "tropical_aoi", "resampled_phtf_loss.tif"), 
         overwrite = TRUE)

# aligned <- brick(aligned_output_name)
# 
# aggregated_2002 <- aggregated[[1]]
# aligned_2002 <- aligned[[1]]
# 
# aggregated_values <- values(aggregated_2002)
# aligned_values <- values(aligned_2002)
# 
# summary(aggregated_values[aggregated_values != 0 & !is.na(aggregated_values)])
# summary(aligned_values[aligned_values != 0 & !is.na(aligned_values)])



### MASK ALWAYS 0 PIXELS

## Create the mask layer
# Create a layer that has values either : NA if phtfloss always 0 across all years, 1 otherwise
phtfloss <- brick(here("temp_data", "processed_phtfloss", "tropical_aoi", "resampled_phtf_loss.tif"))
# # not using if_else here to allow NA as an output...
# always_zero <- function(y){
#   if(sum(y) == 0){d <- NA}else{d <- 1}
#   return(d)}
# don't know why but it will work only with the function like this, converting to 0 and not to NA 
# (which we handle in the masking function next)
always_zero <- function(y){if_else(condition = (sum(y)==0), true = 0, false = 1)}
mask_path <- here("temp_data", "processed_phtfloss", "tropical_aoi", "always_zero_mask_phtfloss.tif")

overlay(x = phtfloss, 
        fun = always_zero, 
        filename = mask_path,
        na.rm = TRUE, # but there is no NA anyway
        overwrite = TRUE)  

mask <- raster(mask_path)
plot(mask)

mask(phtfloss, 
     mask = mask, 
     maskvalue = 0, # necessary here, because the always_zero function used converted to 0 and not to NA
     updatevalue = NA, 
     filename = here("temp_data", "processed_phtfloss", "tropical_aoi", "masked_phtfloss.tif"), 
     overwrite = TRUE)


# masked <- brick(here("temp_data", "processed_phtfloss", "tropical_aoi", "masked_phtfloss.tif"))






