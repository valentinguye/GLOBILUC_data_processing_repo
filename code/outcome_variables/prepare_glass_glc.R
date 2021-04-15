
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
dir.create(here("temp_data", "processed_glass-glc", "tropical_aoi"), recursive = TRUE)


### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")

# GLASS GLC codebook 
# 0 = No data, 10 = cropland, 20 = forest, 30 = Grassland, 40 = Shrubland, 70 = Tundra, 90 = Barren land, 100 = Snow/ice

#### PREPARE GLASS DATA ####

# Read data in
rasterlist <- list.files(path = "input_data/GLASS-GLC", 
                         pattern = paste0("GLASS-GLC_7classes_"), 
                         full.names = TRUE) %>% as.list()
parcels_brick <- brick(rasterlist)

# crop to ***TROPICAL*** AOI 
ext <- extent(c(-180, 179.9167, -30, 30))
tropical_aoi <- crop(parcels_brick, ext)

# Write
writeRaster(tropical_aoi, here("temp_data", "processed_glass-glc", "tropical_aoi", "brick_tropical_aoi.tif"), 
            overwrite = TRUE)


##################
tropical_aoi <- brick( here("temp_data", "processed_glass-glc", "tropical_aoi", "brick_tropical_aoi.tif"))

#### FIRST FOREST LOSS ####
# Create annual layers of forest loss defined as: class is not 20 in a given year while it was 20 in *all the previous year*.
# This restricts annual loss to that occurring for the first time

# construct previous: a collection of annual layers, each giving the mean of the class value in the previous years. 
previous_years <- seq(1982, 2014, 1) 
for(t in 1:length(previous_years)){ # goes only up to 2014, as we don't need the average up to 2015.
  calc(tropical_aoi[[1:t]], fun = mean, 
       filename = here("temp_data", "processed_glass-glc", "tropical_aoi", paste0("past_mean_lu_",previous_years[t], ".tif")), 
       datatype = "FLT4S", # necessary so that a 19.9 mean is not counted as a 20 (i.e. so far undisturbed forest pixel)
       overwrite = TRUE)
}

# this is a raster of 33 layers, giving the mean of GLC class value in 1982, 1982-83, 1982-84, ..., 1982-2014. 
rasterlist <- list.files(path = here("temp_data", "processed_glass-glc", "tropical_aoi"), 
                         pattern = "past_mean_lu_", 
                         full.names = TRUE) %>% as.list()
previous <- brick(rasterlist)

#unique(values(previous))

make_first_loss <- function(previous, current){if_else(condition = (previous == 20 & current != 20), 
                                                       true = 1, false = 0)}

years <- seq(1982, 2015, 1) 
for(t in 2:length(years)){ # starts from 1983 as we need t-1 and thus t starts from 2
  overlay(previous[[t-1]], tropical_aoi[[t]], fun = make_first_loss, 
          filename = here("temp_data", "processed_glass-glc", "tropical_aoi", paste0("first_loss_",years[t], ".tif")), 
          datatype = "INT1U", 
          overwrite = TRUE) 
}

rasterlist <- list.files(path = here("temp_data", "processed_glass-glc", "tropical_aoi"), 
                         pattern = "first_loss_", 
                         full.names = TRUE) %>% as.list()
first_loss <- brick(rasterlist)

#### SUBSEQUENT LAND USE #### 
# Create annual layers indicating the sbqt LU class to forest
# i.e. the class in the year loss is detected. 

make_sbqt_direct_lu <- function(first_loss, current){
  return(first_loss*current)
}
# note that 0 now means "first loss never occurs", and not only "No Data"
sbqt_years <- seq(1983, 2015, 1) 
for(t in 1:length(sbqt_years)){ 
  overlay(first_loss[[t]], tropical_aoi[[t+1]], # since tropical_aoi starts in 1982
          fun = make_sbqt_direct_lu, 
          filename = here("temp_data", "processed_glass-glc", "tropical_aoi", paste0("sbqt_direct_lu_",sbqt_years[t], ".tif")), 
          datatype = "INT1U", 
          overwrite = TRUE) 
}

rasterlist <- list.files(path = here("temp_data", "processed_glass-glc", "tropical_aoi"), 
                         pattern = "sbqt_direct_lu", 
                         full.names = TRUE) %>% as.list()
sbqt_direct_lu <- brick(rasterlist)


### Identify (for each grid cell) the most frequently observed LU in years after a forest loss event. 

## Calculate the mode LU for all grid cells, from any year to the last year. 

# this keeps the first value instance if values occur equally frequently. 
getmode <- function(v, na.rm = TRUE) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

sbqt_years <- seq(1983, 2015, 1) 
for(t in 1:length(sbqt_years)){ # starts from 1983 as we need t-1 and thus t starts from 2
  calc(tropical_aoi[[t:length(sbqt_years)]], fun = getmode, 
       filename = here("temp_data", "processed_glass-glc", "tropical_aoi", paste0("future_mode_lu_",sbqt_years[t], ".tif")), 
       datatype = "INT1U", 
       overwrite = TRUE) 
}
# read it
rasterlist <- list.files(path = here("temp_data", "processed_glass-glc", "tropical_aoi"), 
                         pattern = "future_mode_lu_", 
                         full.names = TRUE) %>% as.list()
mode_lu <- brick(rasterlist)

## Filter the most frequently observed LU from the forest loss event onward (1), where forest loss occurred (2).
# (1) comes from definition of mode_lu, and (2) from the function make_sbqt_mode_lu below
make_sbqt_mode_lu <- function(first_loss, mode_lu){
  return(first_loss*mode_lu)
}
# note that 0 now means "first loss never occurs", and not only "No Data"
sbqt_years <- seq(1983, 2015, 1) 
for(t in 1:length(sbqt_years)){ 
  overlay(first_loss[[t]], mode_lu[[t]], 
          fun = make_sbqt_mode_lu, 
          filename = here("temp_data", "processed_glass-glc", "tropical_aoi", paste0("sbqt_mode_lu_",sbqt_years[t], ".tif")), 
          datatype = "INT1U", 
          overwrite = TRUE) 
}



rasterlist <- list.files(path = here("temp_data", "processed_glass-glc", "tropical_aoi"), 
                         pattern = "sbqt_mode_lu_", 
                         full.names = TRUE) %>% as.list()
sbqt_mode_lu <- brick(rasterlist)


values(tropical_aoi)[17,]
values(first_loss)[17,]
values(mode_lu)[17,]
values(sbqt_mode_lu)[17,]



#### AGGREGATION AND REPROJECTION TO GAEZ EXTENT #### 
# So, here we want to aggregate in turn three variables: 
# 1. The forest loss itself (0;1)
# 2. The subsequent direct LU (0, 10, 20... 100)
# 3. The subsequent mode LU (0. 10, 20... 100)

# For 1. we can just sum --> 0:4 

# For 2. and 3. we reclassify as follows: 
# (A,0,0,0) --> A
# (A, A, 0, 0) --> A 
# etc. 
# (A, B, 0, 0) --> LU mix (AB)
# (A, B, B, 0) --> B
# (A, A, B, B) --> LU mix (AB)
# (A, B, B, B) --> B
# Note: we choose 50 as the value for the "mix" class. 
# na.rm argument is necessary because raster::aggregate needs to find it. It has no effect as it does not appear in the function, 
# and no data is missing.
# also, the function is robust to cases where length(v) > 4, although this should  not happen with fact = 2 in aggregate
aggregate_sbqt <- function(v, na.rm = FALSE){
  uniqv <- unique(v)
  lu <- uniqv[uniqv != 0]
  
  if(length(lu) == 0){ # this is (0, 0, ... 0)
    return(uniqv)
  }
  
  if(length(lu) == 1){ # these are the (A, 0, 0, ... 0) to (A, A, ... A) cases 
    return(lu)
  }
  
  if(length(lu) > 1){
    mode <- DescTools::Mode(v[v != 0]) # we exclude 0 from the mode, because when it is, this is not what we are interested in. 
    # We only want to know what is the most represented LU subsequent to forest loss, i.e. not 0 in sbqt_direct_lu and sbqt_mode_lu. 
    
    # mode is NA if there is exactly one instance of each value (or only the same value, which is excluded here).
    # mode has length > 1 if several values have the same frequencies. 
    if(all(is.na(mode)) | length(mode) > 1){ # so this includes the (A, B, C, ... Z), (A, B, C, ... Z, 0, ... 0), and the (N_A, N_B) with N_A = N_B cases
      return(50)
    }else{# this includes the (A, B, ... B, 0, ...0), (A, B, B, ... B) and (A, B, C, ... C)cases
      return(mode[1])
    }
  }
}

### Input our three variables of interest
# 1. First loss
rasterlist <- list.files(path = here("temp_data", "processed_glass-glc", "tropical_aoi"), 
                         pattern = "first_loss_", 
                         full.names = TRUE) %>% as.list()
first_loss <- brick(rasterlist)

# 2. Subsequent direct LU
rasterlist <- list.files(path = here("temp_data", "processed_glass-glc", "tropical_aoi"), 
                         pattern = "sbqt_direct_lu", 
                         full.names = TRUE) %>% as.list()
sbqt_direct_lu <- brick(rasterlist)

# 3. Subsequent mode LU
rasterlist <- list.files(path = here("temp_data", "processed_glass-glc", "tropical_aoi"), 
                         pattern = "sbqt_mode_lu_", 
                         full.names = TRUE) %>% as.list()
sbqt_mode_lu <- brick(rasterlist)


### Input GAEZ template
gaez <- raster(here("temp_data", "GAEZ", "AES_index_value", "Rain-fed", "High-input", "Banana.tif"))

### 1. FIRST FOREST LOSS
aggregate(first_loss, 
          fact = 2, 
          fun = sum, 
          expand = FALSE, 
          na.rm = FALSE, 
          filename = here("temp_data", "processed_glass-glc", "tropical_aoi", "aggr_first_loss.tif"),
          datatype = "INT1U",
          overwrite = TRUE)

aggr_first_loss <- brick(here("temp_data", "processed_glass-glc", "tropical_aoi", "aggr_first_loss.tif"))

# and now align precisely to GAEZ 

beginCluster() # this uses by default detectCores() - 1

resample(x = aggr_first_loss, y = gaez,
         method = "ngb",
         filename = here("temp_data", "processed_glass-glc", "tropical_aoi", "resampled_first_loss.tif"),
         datatype = "INT1U",
         overwrite = TRUE )

endCluster()

### 2. SUBSEQUENT DIRECT LU 
aggregate(sbqt_direct_lu, 
          fact = 2, 
          fun = aggregate_sbqt, 
          expand = FALSE, 
          na.rm = FALSE, 
          filename = here("temp_data", "processed_glass-glc", "tropical_aoi", "aggr_sbqt_direct_lu.tif"),
          datatype = "INT1U",
          overwrite = TRUE)

aggr_sbqt_direct_lu <- brick(here("temp_data", "processed_glass-glc", "tropical_aoi", "aggr_sbqt_direct_lu.tif"))

# and now align precisely to GAEZ 

beginCluster() # this uses by default detectCores() - 1

resample(x = aggr_sbqt_direct_lu, y = gaez,
         method = "ngb",
         filename = here("temp_data", "processed_glass-glc", "tropical_aoi", "resampled_sbqt_direct_lu.tif"),
         datatype = "INT1U",
         overwrite = TRUE )

endCluster()

### 3. SUBSEQUENT MODE LU
aggregate(sbqt_mode_lu, 
          fact = 2, 
          fun = aggregate_sbqt, 
          expand = FALSE, 
          na.rm = FALSE, 
          filename = here("temp_data", "processed_glass-glc", "tropical_aoi", "aggr_sbqt_mode_lu.tif"),
          datatype = "INT1U",
          overwrite = TRUE)

aggr_sbqt_mode_lu <- brick(here("temp_data", "processed_glass-glc", "tropical_aoi", "aggr_sbqt_mode_lu.tif"))

# and now align precisely to GAEZ 

beginCluster() # this uses by default detectCores() - 1

resample(x = aggr_sbqt_mode_lu, y = gaez,
         method = "ngb",
         filename = here("temp_data", "processed_glass-glc", "tropical_aoi", "resampled_sbqt_mode_lu.tif"),
         datatype = "INT1U",
         overwrite = TRUE)

endCluster()


#### MASK ALWAYS 0 PIXELS #### 

### Create the mask layer
# Create a layer that has values either : NA if first_loss always 0 across all years, 1 otherwise
res_first_loss <- brick(here("temp_data", "processed_glass-glc", "tropical_aoi", "resampled_first_loss.tif"))
# not using if_else here to allow NA as an output... 
always_zero <- function(y){
  if(sum(y) == 0){d <- NA}else{d <- 1}
  return(d)}

mask_path <- here("temp_data", "processed_glass-glc", "tropical_aoi", "always_zero_mask.tif")

overlay(x = res_first_loss, 
         fun = always_zero, 
         filename = mask_path,
         na.rm = TRUE, # but there is no NA anyway
         datatype = "INT2U", # INT2U to allow have NAs
         overwrite = TRUE)  

mask <- raster(mask_path)
# plot(mask)

# then use it to mask the 3 brick variables

### 1. FIRST LOSS 

mask(x = res_first_loss, 
     mask = mask, 
     filename = here("temp_data", "processed_glass-glc", "tropical_aoi", "masked_first_loss.tif"), 
     datatype = "INT2U")

# masked <- brick( here("temp_data", "processed_glass-glc", "tropical_aoi", "masked_first_loss.tif"))
# plot(masked[[1]])

### 2. SUBSEQUENT DIRECT LU
res_sbqt_direct_lu <- brick(here("temp_data", "processed_glass-glc", "tropical_aoi", "resampled_sbqt_direct_lu.tif"))

mask(x = res_sbqt_direct_lu, 
     mask = mask, 
     filename = here("temp_data", "processed_glass-glc", "tropical_aoi", "masked_sbqt_direct_lu.tif"), 
     datatype = "INT2U")

### 3. SUBSEQUENT MODE LU
res_sbqt_mode_lu <- brick(here("temp_data", "processed_glass-glc", "tropical_aoi", "resampled_sbqt_mode_lu.tif"))

mask(x = res_sbqt_mode_lu, 
     mask = mask, 
     filename = here("temp_data", "processed_glass-glc", "tropical_aoi", "masked_sbqt_mode_lu.tif"), 
     datatype = "INT2U")





plot(sbqt_mode_lu[[20]])
plot(aggr_sbqt_mode_lu[[20]])
plot(res_sbqt_mode_lu[[20]])


# dis <- sbqt_mode_lu[[20]]
# hist(values(dis)[values(dis)!=0])
# 
# aggrbrick <- sbqt_mode_lu_aggr[[20]]
# plot(aggrbrick)
# hist(values(aggrbrick)[values(aggrbrick)!=0])
# 
# aggr <- aggregate(dis, fact = 2, fun = aggregate_sbqt, expand = FALSE, na.rm = FALSE)
# proj <- projectRaster(dis, gaez, method = "ngb") # much less precise than aggregate_sbqt it seems
# 
# aligned <- resample(aggr, gaez, method = "ngb") # resample equivalent to projectRaster here.
# 
# all.equal(values(aggr), values(aggrbrick))
# plot(dis)
# plot(aggr)
# plot(aggrbrick)
# plot(proj)
# plot(aligned)
# plot(proj_aggr)
# 
# plot(sbqt_mode_lu[[20]])
# plot(sbqt_mode_lu_aggr[[20]])






#### MASK ALWAYS 0 PIXELS IN GLASS RESOLUTION (FOR GEE) ####

# Input first loss in Glass resolution 
rasterlist <- list.files(path = here("temp_data", "processed_glass-glc", "tropical_aoi"), 
                         pattern = "first_loss_", 
                         full.names = TRUE) %>% as.list()
first_loss <- brick(rasterlist)

always_zero2 <- function(y){
  if(sum(y) == 0){d <- 0}else{d <- 1} # note that value for true is not NA but 0, for further convenience in GEE
  return(d)}

mask_path <- here("temp_data", "processed_glass-glc", "tropical_aoi", "always_zero_mask_glassres.tif")

overlay(x = first_loss, 
        fun = always_zero2, 
        filename = mask_path,
        na.rm = TRUE, # but there is no NA anyway
        datatype = "INT2U", # INT2U to allow have NAs
        overwrite = TRUE)  

mask <- raster(mask_path)
plot(mask)






