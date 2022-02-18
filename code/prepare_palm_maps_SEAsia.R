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
dir.create(here("temp_data", "processed_aop", "SEAsia_aoi"), recursive = TRUE)

### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")


### SOUTH EAST ASIA AOI 
aop <- raster(here("input_data", "aop_lower", "lower", "2001_op_lower.tif"))
SEA_ext <- extent(aop) %>% as('SpatialPolygons')
SEA_ext <- st_as_sfc(SEA_ext)
st_crs(SEA_ext) <- crs(aop)
SEA_ext <- st_transform(SEA_ext, crs = 4326)
SEA_ext
r <- raster(here("input_data", "10thLossSumGlass_maxP_SEAsia_2001.tif"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# make middle point for those years when a lower and an upper bounds are provided
# make it as the mean, such that the "midpoint" will be 0.5 for those pixels where detecting oil palm is not 100% sure. 
# Since we aggregate afterwards, there is not pb of interpretation. 

# IT IS NOT ALL YEARS
rasterlist <- list.files(path = here("input_data", "aop_lower", "lower"), 
                         pattern = "op_lower", 
                         full.names = TRUE) %>% as.list()
# DO NOT BRICK, IT DOES NOT WORK, rather stack (don't know why)
aoplb <- stack(rasterlist)

rasterlist <- list.files(path = here("input_data", "aop_upper", "upper"), 
                         pattern = "op_upper", 
                         full.names = TRUE) %>% as.list()
aopub <- stack(rasterlist)

names(aopub)
# both aoplb and aopub are sorted in alphabetical order, as per list.files

midpoint_years <- c("2001", "2002", "2003", "2004", "2005", "2006", "2011", "2012", "2013", "2014")

beginCluster() # uses by default detectedCores() - 1

for(year in midpoint_years){
  
  # the order in which ub and lb are stacked does not matter
  rs <- stack(aoplb[[paste0("X",year,"_op_lower")]], aopub[[paste0("X",year,"_op_upper")]])
  
  clusterR(rs,
           fun = calc, 
           args = list(mean),
           filename =  here("temp_data", "processed_aop", "SEAsia_aoi", paste0(year,"_op.tif")), 
           overwrite = TRUE )

}
endCluster()

r <- raster(here("temp_data", "processed_aop", "SEAsia_aoi", paste0(year,"_op.tif")))

#saop <- sampleRandom(aoplb[[1]], 10000)
#unique(saop[,1])








#### MAKE FOREST LOSS IN SOUTH EAST ASIA AOI #### 
# import annual layers of forest loss (in hectares) as computed in GEE (and downloaded from Google Drive to input_data/)
