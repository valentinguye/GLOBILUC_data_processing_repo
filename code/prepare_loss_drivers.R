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
dir.create(here("temp_data", "processed_lossdrivers", "tropical_aoi"), recursive = TRUE)

### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")

### TROPICAL AOI 
ext <- extent(c(-180, 179.9167, -30, 30))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 



#### PREPARE LOSS #### 

### Brick layers ### 

# import annual layers of forest loss (in hectares) as computed in GEE (and downloaded from Google Drive to input_data/)
rasterlist <- list.files(path = here("input_data", "10thLossSumGlass_maxP"), 
                         pattern = "^10thLossSumGlass_maxP_", 
                         full.names = TRUE) %>% as.list()
lossdrivers <- brick(rasterlist)

# test <- lossdrivers[[1]]
# test
# plot(test)
# summary(values(test))

writeRaster(lossdrivers, here("temp_data", "processed_lossdrivers", "tropical_aoi", "loss_drivers.tif"), 
            overwrite = TRUE)

#### PREPARE DRIVERS #### 

drivers <- raster(here("input_data", "curtis", "Goode_FinalClassification_19_05pcnt_prj", "Goode_FinalClassification_19_05pcnt_prj.tif"))

# crop it to the tropical aoi
drivers <- crop(drivers, ext)


# We need to map all values that correspond to other drivers than agricultural ones, to a single value
# drivers data set has NA, 1:5 values 
# according code in https://www.sustainabilityconsortium.org/tsc-downloads/forest-data/?ind=1536355391488&filename=Supplementary_Dataset_1.txt&wpdmdl=24673&refresh=60f2aa214d3801626516001
# (see "Create final classification using model output" section)
# 1 is commodity driven class
# 2 is shifting agriculture class
# 3 is forestry
# 4 is fires
# 5 is urbanization
# 0 is minor loss
# NA is mask for water and <.5% loss 

make_binary <- function(x){if_else(condition = (x==1 | x==2), true = 1, false = 0)}

drivers_b_output_name <- here("temp_data", "processed_lossdrivers", "tropical_aoi", "drivers_binary.tif")

calc(drivers, 
     fun = make_binary,
     filename = drivers_b_output_name, 
     overwrite = TRUE)


#### AGGREGATE AND alignE TO THE DRIVERS RESOLUTION ####

lossdrivers <- brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", "loss_drivers.tif"))

# Currently, the data are hectares of forest loss in 5x5km grid cells annually. 
# We aggregate this to 10km grid cells by adding up the hectares of loss

# define output file name
aggr_output_name <- here("temp_data", "processed_lossdrivers", "tropical_aoi", "aggr_loss_drivers.tif")

# aggregate it from the ~5km cells to ~10km
raster::aggregate(lossdrivers, fact = c(res(drivers)[1]/res(lossdrivers)[1], res(drivers)[2]/res(lossdrivers)[2]),
                  expand = FALSE,
                  fun = sum,
                  na.rm = TRUE, # NA values are only in the sea. Where there is no forest loss, like in a city in Brazil, the value is 0 (see with plot())
                  filename = aggr_output_name,
                  # datatype = "INT2U", # let the data be float, as we have decimals in the amount of hectares. 
                  overwrite = TRUE)

# align to DRIVERS exactly 
aggregated <- brick(aggr_output_name)

resampled_output_name <- here("temp_data", "processed_lossdrivers", "tropical_aoi", "drivers_resampled_loss.tif")

resample(x = aggregated, 
         y = drivers, 
         method = "ngb", # we use ngb and not bilinear because the output values' summary better fits that of the aggregated layer 
         # and the bilinear interpolation arguably smoothes the reprojection more than necessary given that from and to are already very similar.  
         filename = resampled_output_name, 
         overwrite = TRUE)


### MASK LOSS WITH DRIVERS ### 

drivers_b <- raster(drivers_b_output_name)
resampled <- brick(resampled_output_name)

agri_lossdrivers_output_name <- here("temp_data", "processed_lossdrivers", "tropical_aoi", "agri_lossdrivers.tif")

mask(resampled, 
     mask = drivers_b, 
     maskvalue = 0, 
     updatevalue = 0, 
     filename = agri_lossdrivers_output_name, 
     overwrite = TRUE)




#### ALIGN to GAEZ ####

# Read in a GAEZ grid, because we want to align to it. 
gaez <- raster(here("temp_data", "GAEZ", "AES_index_value", "Rain-fed", "High-input", "Alfalfa.tif"))

agri_lossdrivers <- brick(agri_lossdrivers_output_name)

gaez_resampled_output_name <- here("temp_data", "processed_lossdrivers", "tropical_aoi", "resampled_loss_drivers.tif")

resample(x = agri_lossdrivers, 
         y = gaez, 
         method = "ngb", # we use ngb and not bilinear because the output values' summary better fits that of the aggregated layer 
         # and the bilinear interpolation arguably smoothes the reprojection more than necessary given that from and to are already very similar.  
         filename = gaez_resampled_output_name, 
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
# Create a layer that has values either : NA if lossdrivers always 0 across all years, 1 otherwise
final_lossdrivers <- brick(gaez_resampled_output_name)
# # not using if_else here to allow NA as an output...
# always_zero <- function(y){
#   if(sum(y) == 0){d <- NA}else{d <- 1}
#   return(d)}
# don't know why but it will work only with the function like this, converting to 0 and not to NA 
# (which we handle in the masking function next)
always_zero <- function(y){if_else(condition = (sum(y)==0), true = 0, false = 1)}
mask_path <- here("temp_data", "processed_lossdrivers", "tropical_aoi", "always_zero_mask_lossdrivers.tif")

overlay(x = final_lossdrivers, 
        fun = always_zero, 
        filename = mask_path,
        na.rm = TRUE, # but there is no NA anyway
        overwrite = TRUE)  

mask <- raster(mask_path)
# plot(mask)

mask(final_lossdrivers, 
     mask = mask, 
     maskvalue = 0, # necessary here, because the always_zero function used converted to 0 and not to NA
     updatevalue = NA, 
     filename = here("temp_data", "processed_lossdrivers", "tropical_aoi", "masked_lossdrivers.tif"), 
     overwrite = TRUE)


# masked <- brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", "masked_lossdrivers.tif"))






