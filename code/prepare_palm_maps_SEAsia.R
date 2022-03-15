##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 

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
dir.create(here("temp_data", "processed_aop", "SEAsia_aoi"), recursive = TRUE)
dir.create(here("temp_data", "processed_aop", "SEAsia_aoi", "convenience_folder"), recursive = TRUE)
dir.create(here("temp_data", "processed_aop", "SEAsia_aoi", "unidir_making"), recursive = TRUE)

### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")


### SOUTH EAST ASIA AOI 
aop <- raster(here("input_data", "aop_lower", "lower", "2001_op_lower.tif"))
SEA_ext <- extent(aop)
SEA_ext <- as(SEA_ext, 'SpatialPolygons')
SEA_ext <- st_as_sfc(SEA_ext)
st_crs(SEA_ext) <- crs(aop)
SEA_ext <- st_transform(SEA_ext, crs = 4326)
SEA_ext
SEA_aoi <- extent(c(94.955908, 120.404282, -5.941733, 7.363556 ))
# r <- raster(here("input_data", "10thLossSumGlass_maxP_SEAsia_2001.tif"))


all_years <- c(2001:2016) %>% as.character()

midpoint_years <- c("2001", "2002", "2003", "2004", "2005", "2006", "2011", "2012", "2013", "2014")


### DRIVERS CROPED AT TROPICAL AOI
# this is necessary in other parts than prepare drivers, because it is the raster target 
drivers <- raster(here("input_data", "curtis", "Goode_FinalClassification_19_05pcnt_prj", "Goode_FinalClassification_19_05pcnt_prj.tif"))

# crop it to the tropical aoi
drivers <- crop(drivers, SEA_aoi)


#### TURN ALL RASTERS TO INT1U #### 
# Do it now already such that no FLT4S is created/manipulated 

# lower bound layers (not all years)
rasterlist <- list.files(path = here("input_data", "aop_lower", "lower"), 
                         pattern = "op_lower", 
                         full.names = TRUE) %>% as.list()
# DO NOT BRICK, IT DOES NOT WORK, rather stack (don't know why)
aoplb <- stack(rasterlist)

for(l in 1:nlayers(aoplb)){
  calc(aoplb[[l]], 
       fun = function(x){x*100}, 
       filename = here("temp_data", "processed_aop", "SEAsia_aoi", paste0(names(aoplb[[l]]), ".tif")), 
       datatype = "INT1U",
       overwrite = TRUE)
}

# upper bound AND "no midpoint" layers 
rasterlist <- list.files(path = here("input_data", "aop_upper", "upper"), 
                         pattern = "op", # note the difference here with the call above 
                         full.names = TRUE) %>% as.list()
aopub <- stack(rasterlist)

for(l in 1:nlayers(aopub)){
  calc(aopub[[l]], 
       fun = function(x){x*100}, 
       filename = here("temp_data", "processed_aop", "SEAsia_aoi", paste0(names(aopub[[l]]), ".tif")), 
       datatype = "INT1U",
       overwrite = TRUE)
}



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### MIDPOINTS #### 
# make middle point for those years when a lower and an upper bounds are provided
# make it as the mean, such that the "midpoint" will be 0.5 for those pixels where detecting oil palm is not 100% sure. 
# Since we aggregate afterwards, there is not pb of interpretation. 

# IT IS NOT ALL YEARS
rasterlist <- list.files(path = here("temp_data", "processed_aop", "SEAsia_aoi"), 
                         pattern = "op_lower", 
                         full.names = TRUE) %>% as.list()
# DO NOT BRICK, IT DOES NOT WORK, rather stack (don't know why)
aoplb <- stack(rasterlist)

rasterlist <- list.files(path = here("temp_data", "processed_aop", "SEAsia_aoi"), 
                         pattern = "op_upper", 
                         full.names = TRUE) %>% as.list()
aopub <- stack(rasterlist)
dataType(aopub)
names(aopub)
#saop <- sampleRandom(aoplb[[1]], 10000)
#unique(saop[,1])

# both aoplb and aopub are sorted in alphabetical order, as per list.files

beginCluster() # uses by default detectedCores() - 1

for(year in midpoint_years){
  
  # the order in which ub and lb are stacked does not matter
  rs <- stack(aoplb[[paste0("X",year,"_op_lower")]], aopub[[paste0("X",year,"_op_upper")]])
  
  clusterR(rs,
           fun = calc, 
           args = list(mean),
           filename =  here("temp_data", "processed_aop", "SEAsia_aoi", paste0("X", year,"_op_midpoint.tif")), 
           datatype = "INT1U",
           overwrite = TRUE )

}
endCluster()

# MOVE pre-processed annual maps to separate folder
for(year in midpoint_years){
  file.copy(from = here("temp_data", "processed_aop", "SEAsia_aoi", paste0("X", year, "_op_midpoint.tif")), 
            to = here("temp_data", "processed_aop", "SEAsia_aoi", "convenience_folder", paste0("op_", year, ".tif")), 
            overwrite = TRUE)
}

other_years <- all_years[!(all_years %in% midpoint_years)]
for(year in other_years){
  file.copy(from = here("temp_data", "processed_aop", "SEAsia_aoi", paste0("X", year, "_op.tif")), 
            to = here("temp_data", "processed_aop", "SEAsia_aoi", "convenience_folder", paste0("op_", year, ".tif")), 
            overwrite = TRUE)
}


#### IMPOSE UNI-DIRECTIONAL CHANGE #### 
# Such that variation from oil palm plantations that stop being observed is not included in regressions

rasterlist <- list.files(path = here("temp_data", "processed_aop", "SEAsia_aoi", "convenience_folder"), 
                         pattern = "op", 
                         full.names = TRUE) %>% as.list()
aop <- stack(rasterlist)

make_wasNever <- function(x){return(100 - max(x))}
# max(x) is a single layer raster, with values either: 
# 100 if OP detected in previous years, 
# 50 if OP only detected in upper bound in previous midpoint years, 
# or 0 if OP never detected before

# thus, wasNever can be 50, meaning that there is 50% chance that the pixel was actually never cover with OP

beginCluster() # uses by default detectedCores() - 1

for(i in 2:nlayers(aop)){
  previous <- aop[[1:(i-1)]]
  
  clusterR(previous,
           fun = calc, 
           args = list(make_wasNever),
           filename =  here("temp_data", "processed_aop", "SEAsia_aoi", "unidir_making", paste0("wasNever_", all_years[i],".tif")), 
           datatype = "INT1U",
           overwrite = TRUE)
  
  removeTmpFiles(h = 0)
}
endCluster()

# multiplying such 50 pixels with OP extent in year of interest: 
# if there is no OP, then there's just no OP expansion this year, no matter the past history in the pixel
# if there is OP (100) then, the output from the overlay is 100*50*0.01 = 50, which can be interpreted as there is new expansion with 50% certainty, because we know with 50% chance only that there was no OP before.
# if there is OP with 50% certainty (50);, then the output is 25, which can be interpreted as there is expansion in the pixel only if it's actually OP this year (50% chance), and there was actually no oil palm before (50% chance). 
# but that yields a 25 pixel value for all subsequent years until there is certain OP, which is incorrect, as we don't want more than 100% of the pixel to be converted over the full period. 
# thus, reclassify 25 to 0, to count 50% expansion the first time 50 is observed, and the remaining 50 when full pixel (or 100% certainty) is observed. 
make_annualBinary <- function(x, y){return(x*y*0.01)}

# # fake rasters to test 
# wasNever <- r1*100
# aop <- s*100
# i <- 2

beginCluster()
for(i in 2:nlayers(aop)){
  
   wasNever <- raster(here("temp_data", "processed_aop", "SEAsia_aoi", "unidir_making", paste0("wasNever_", all_years[i],".tif")))
   
   rs <- stack(aop[[i]], wasNever)
   
   clusterR(rs, 
            fun = overlay, 
            args = list(fun = make_annualBinary),
            filename = here("temp_data", "processed_aop", "SEAsia_aoi", "unidir_making", paste0("annualBinary_", all_years[i], ".tif")), 
            datatype = "INT1U", 
            overwrite = TRUE)
  
  removeTmpFiles(h=0)
}
endCluster()

for(i in 2:nlayers(aop)){
    
  annualBinary <- raster(here("temp_data", "processed_aop", "SEAsia_aoi", "unidir_making", paste0("annualBinary_", all_years[i], ".tif")))
  
  reclassify(annualBinary, cbind(25, 0), 
             filename = here("temp_data", "processed_aop", "SEAsia_aoi", "unidir_making", paste0("unidir_", all_years[i], ".tif")), 
             datatype = "INT1U", 
             overwrite = TRUE)
  
  removeTmpFiles(h=0)
}



### AGGREGATE AND ALIGNE TO DRIVERS #### 

# Function description
# The function has for inputs annual layers of lucfp events at the pixel level.
# It aggregates these pixels to a parcel size defined by parcel_size (in meters).
# The aggregation operation is the sum of the pixel lucfp events.
# Each annual aggregation is tasked in parallel.

## sequence over which to execute the task.
# We attribute the tasks to CPU "workers" at the annual level and not at the pf_type level.
# Hence, if a worker is done with its annual task before the others it can move on to the next one and workers' labor is maximized wrt.
# attributing tasks at the pf_type level.

unidir_years <- seq(from = 2002, to = 2016, by = 1)

aggregate_pixels <- function(x, na.rm = na.rm){sum(x, na.rm = na.rm)*0.01}

## read the input to the task
# is done within each task because it is each time different here.

## define the task
aggregate_annual <- function(time){
  # Define which process we are in:
  processname <-  here("temp_data", "processed_aop", "SEAsia_aoi", "unidir_making", paste0("unidir_", unidir_years[time], ".tif"))
  
  #set temp directory
  dir.create(paste0(processname,"_Tmp"), showWarnings = FALSE)
  rasterOptions(tmpdir=here(paste0(processname,"_Tmp")))
  
  # read in the input.
  unidir_annual <- raster(processname)
  
  # define output file name
  output_filename <- here("temp_data", "processed_aop", "SEAsia_aoi", paste0("unidir_9km_", unidir_years[time], ".tif"))
  
  raster::aggregate(unidir_annual, fact = 100, # multiplies in both directions by 100. Since current resolution is 100m, target resolution is 10000m, i.e. 10km
                    expand = FALSE,
                    fun = aggregate_pixels,
                    na.rm = FALSE, # NA cells are in margins, see the NOTES part. If FALSE, aggregations at margins that use NA 
                    # are discarded because the sum would be spurious as it would count all NA as 0s while it is not necessary the case.
                    filename = output_filename,
                    datatype = "INT2U", # because the sum may go up to 8100 with aggregation to 9km,
                    overwrite = TRUE)  
  #removes entire temp directory without affecting other running processes (but there should be no temp file now)
  unlink(paste0(processname,"_Tmp"), recursive = TRUE)
  # return the path to this parcels file
  #return(output_filename)
}

## register cluster
registerDoParallel(cores = detectCores() - 1)

##  define foreach object.
foreach(t = 1:length(unidir_years),
        # .combine combine the outputs as a mere character list (by default)
        .inorder = FALSE, # we don't care that the results be combine in the same order they were submitted
        .multicombine = TRUE,
        .export = c("aggregate_pixels", "unidir_years"),
        .packages = c("raster", "rgdal", "here")
) %dopar% aggregate_annual(time = t)



# align to DRIVERS exactly 
rasterlist <- list.files(path = here("temp_data", "processed_aop", "SEAsia_aoi"), 
                         pattern = "unidir_9km_", 
                         full.names = TRUE) %>% as.list()
aggregated <- brick(rasterlist)

op_resampled_output_name <- here("temp_data", "processed_aop", "SEAsia_aoi", "resampled_unidir_0216.tif")

resample(x = aggregated, 
         y = drivers, 
         method = "bilinear", # bilinear or ngb changes nothing 
         filename = pastures_resampled_output_name, 
         overwrite = TRUE)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# OLD STUFF 

# # the rationale is: operate at pixel level, across years. Replace any value (equivalently any 0) by 1 after the first time pixel value is 1. 
# # how to deal with 0.5 midpoints? 
# # if there is a midpoint after a 1, it gets replaced by 1, that is no pb. 
# # if there are midpoints before 1, don't change them: it was not sure that it was oil palm, and new expansion made it surer. Doing so counts the new expansion as the remaining 50% of grid cell. 
# 
# # note that "When using RasterLayer objects, the number of arguments of the function should match the number of Raster objects, 
# # or it should take any number of arguments."
# # so calc is more appropriate
# 
# # compute this quantity in advance so it's not computed for every cell
# aop_length <- nlayers(aop)


# test again how it handles missing and large values other than 1 

## TEST ZONE ON UNI-DIRECTIONAL FUNCTION
# r <- raster(ncol=5, nrow=5)
# r1 <- init(r, fun=sample(x=c(0,1, 0.5), replace = TRUE, size = 25))
# r2 <- init(r, fun=sample(x=c(0,1, 0.5), replace = TRUE, size = 25))
# r3 <- init(r, fun=sample(x=c(0,1, 0.5), replace = TRUE, size = 25))
# r4 <- init(r, fun=sample(x=c(0,1, 0.5), replace = TRUE, size = 25))
# r5 <- init(r, fun=sample(x=c(0,1, 0.5), replace = TRUE, size = 25))
# r6 <- init(r, fun=sample(x=c(0,1, 0.5), replace = TRUE, size = 25))
# r7 <- init(r, fun=sample(x=c(0,1, 0.5), replace = TRUE, size = 25))
# r8 <- init(r, fun=sample(x=c(0,1, 0.5), replace = TRUE, size = 25))
# 
# s <- stack(r1, r2, r3, r4, r5, r6, r7, r8)
# 
# 
# aopunidir <- list()
# 
# ### OR other method, closer to what is done in GEE for mapbiomass and soy unidirectional 
# for(i in 2:nlayers(s)){
#   previous <- s[[1:(i-1)]]
#   previousMax <- max(previous) # so this is a single layer raster, with values either: 
#   # 1 if OP detected in previous years, 
#   # 0.5 if OP only detected in upper bound in previous midpoint years, 
#   # or 0 if OP never detected before
#   
#   # thus, wasNever can be 0.5, meaning that there is 50% chance that the pixel was actually never cover with OP (50% of the pixel area)
#   wasNever <- 1 - previousMax
#   # multiplying such 0.5 pixels with OP extent in year of interest: 
#   # if there is no OP, then there's just no OP expansion this year, no matter the past history in the pixel
#   # if there is OP (1) then, the output is 0.5, which can be interpreted as there is new expansion with 50% certainty, because we know with 50% chance only that there was no OP before.
#   # if there is OP with 50% certainty (0.5);, then the output is 0.25, which can be interpreted as there is expansion in the pixel only if it's actually OP this year (50% chance), and there was actually no oil palm before (50% chance). 
#   # but that yields a 0.25 pixel value for all subsequent years until there is certain OP, which is incorrect, as we don't want more than 100% of the pixel to be converted over the full period. 
#   # thus, reclassify 0.25 to 0, to count 50% expansion the first time 0.5 is observed, and the remaining 0.5 when full pixel (or 100% certainty) is observed. 
#   annualBinary <- overlay(s[[i]], wasNever, fun = function(x, y){return(x*y)})
#   annualBinary <- reclassify(annualBinary, cbind(0.25, 0))
#   aopunidir[[paste0("unidir_", all_years[i])]] <- annualBinary
# }
# 
# 
# unidir <- stack(aopunidir)
# recl <- values(stack(s, unidir))
# 
# 
# 
# impose_unidir <- function(y){ index_1 <- which(y==1)
#                               index_05 <- which(y==0.5)
#                               lth_index_05 <- length(index_05)
#                               if(length(index_1)>0){
# 
#                                 y[min(index_1):length(y)] <- 1
# 
#                                 # all these conditions are necessary to handle cases where index_* objects are empty.
#                                 if(lth_index_05>0){
#                                   if(min(index_05) < min(index_1)){
#                                     # in this case, there is uncertain oil palm, then no oil palm, then oil palm. 
#                                     # we don't want to include variation from 0.5 to 0, so make 0.5 remain until it's 1
#                                     y[min(index_05):(min(index_1)-1)] <- 0.5
#                                     
#                                   }
#                                 }
#                               } else if (lth_index_05>0){ # this is the case where there is only uncertain, midpoint values (0.5s)
#                                 y[min(index_05):length(y)] <- 0.5
#                               }
# 
#                               return(y)
# }
# 
# test <- calc(s, fun = impose_unidir)
# values(stack(s, test))
# 
# 
# 
# beginCluster()
# test <- clusterR(s,
#                  fun = calc,
#                  args = list(impose_unidir)
# )
# 
# endCluster()
# values(stack(s, test))
# 
# b <- rep(c(0, 0.5, 0,0, 0.5, 0,0,0), 2)
# impose_unidir(b)
# r1 <- matrix(b, 4, 4)
# r1 <- raster(r1)
# 
# b0 <- rep(0, 16)
# impose_unidir(b0)



# use overlay with unstack = TRUE










#### MAKE FOREST LOSS IN SOUTH EAST ASIA AOI #### 
# import annual layers of forest loss (in hectares) as computed in GEE (and downloaded from Google Drive to input_data/)
