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
# r <- raster(here("input_data", "10thLossSumGlass_maxP_SEAsia_2001.tif"))


all_years <- c(2001:2016) %>% as.character()

midpoint_years <- c("2001", "2002", "2003", "2004", "2005", "2006", "2011", "2012", "2013", "2014")


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
# this writes rasters in FLT4S (default) which is needed given that we require values to be 0.5. This makes files much larger. 
# Depending on type of aggregation we can choose another value for midpoints. 

#saop <- sampleRandom(aoplb[[1]], 10000)
#unique(saop[,1])

# this is for automation, but one can also simply copy paste the files from input data to temp data folders
other_years <- all_years[!(all_years %in% midpoint_years)]
for(year in other_years){
  r <- raster(here("input_data", "aop_upper", "upper", paste0(year, "_op.tif")))
  dataType(r) <- "FLT4S"
  writeRaster(r, 
              here("temp_data", "processed_aop", "SEAsia_aoi", paste0(year,"_op.tif")), 
              overwrite = TRUE) # important to match the format from midpoint years
}
# r <- raster(here("temp_data", "processed_aop", "SEAsia_aoi", paste0(year,"_op.tif")))




#### IMPOSE UNI-DIRECTIONAL CHANGE #### 
# Such that variation from oil palm plantations that stop being observed is not included in regressions

rasterlist <- list.files(path = here("temp_data", "processed_aop", "SEAsia_aoi"), 
                         pattern = "_op", 
                         full.names = TRUE) %>% as.list()
aop <- stack(rasterlist)

# the rationale is: operate at pixel level, across years. Replace any value (equivalently any 0) by 1 after the first time pixel value is 1. 
# how to deal with 0.5 midpoints? 
# if there is a midpoint after a 1, it gets replaced by 1, that is no pb. 
# if there are midpoints before 1, don't change them: it was not sure that it was oil palm, and new expansion made it surer. Doing so counts the new expansion as the remaining 50% of grid cell. 

# note that "When using RasterLayer objects, the number of arguments of the function should match the number of Raster objects, 
# or it should take any number of arguments."
# so calc is more appropriate

# compute this quantity in advance so it's not computed for every cell
aop_length <- nlayers(aop)
impose_unidir <- function(y){ index_1 <- which(y==1)
                              index_05 <- which(y==0.5)
                              lth_index_05 <- length(index_05)
                              if(length(index_1)>0){
                                
                                y[min(index_1):16] <- 1
                                
                                # all these conditions are necessary to handle cases where index_* objects are empty.
                                if(lth_index_05>0){
                                  if(min(index_05) < min(index_1)){
                                    # in this case, there is uncertain oil palm, then no oil palm, then oil palm. 
                                    # we don't want to include variation from 0.5 to 0, so make 0.5 remain until it's 1
                                    y[min(index_05):(min(index_1)-1)] <- 0.5
                                    
                                  }
                                }
                              } else if (lth_index_05>0){ # this is the case where there is only uncertain, midpoint values (0.5s)
                                y[min(index_05):16] <- 0.5
                              }
                              
                              return(y)
}

calc(aop, impose_unidir, 
     filename =  here("temp_data", "processed_aop", "SEAsia_aoi", "unidir_brick.tif"), 
     overwrite = TRUE)

dataType(pb)

uni <- brick(here("temp_data", "processed_aop", "SEAsia_aoi", "unidir_op_brick.tif"))

uni@data@max
uni@data@min
pb <- uni$unidir_op_brick.12
pb@data

plot(uni$unidir_op_brick.6)




# test again how it handles missing and large values other than 1 

## TEST ZONE ON UNI-DIRECTIONAL FUNCTION
r <- raster(ncol=5, nrow=5)
r1 <- init(r, fun=sample(x=c(0,1, 0.5), replace = TRUE, size = 25))
r2 <- init(r, fun=sample(x=c(0,1, 0.5), replace = TRUE, size = 25))
r3 <- init(r, fun=sample(x=c(0,1, 0.5), replace = TRUE, size = 25))
r4 <- init(r, fun=sample(x=c(0,1, 0.5), replace = TRUE, size = 25))
r5 <- init(r, fun=sample(x=c(0,1, 0.5), replace = TRUE, size = 25))
r6 <- init(r, fun=sample(x=c(0,1, 0.5), replace = TRUE, size = 25))
r7 <- init(r, fun=sample(x=c(0,1, 0.5), replace = TRUE, size = 25))
r8 <- init(r, fun=sample(x=c(0,1, 0.5), replace = TRUE, size = 25))

s <- stack(r1, r2, r3, r4, r5, r6, r7, r8)
impose_unidir <- function(y){ index_1 <- which(y==1)
                              index_05 <- which(y==0.5)
                              lth_index_05 <- length(index_05)
                              if(length(index_1)>0){

                                y[min(index_1):length(y)] <- 1

                                # all these conditions are necessary to handle cases where index_* objects are empty.
                                if(lth_index_05>0){
                                  if(min(index_05) < min(index_1)){
                                    # in this case, there is uncertain oil palm, then no oil palm, then oil palm. 
                                    # we don't want to include variation from 0.5 to 0, so make 0.5 remain until it's 1
                                    y[min(index_05):(min(index_1)-1)] <- 0.5
                                    
                                  }
                                }
                              } else if (lth_index_05>0){ # this is the case where there is only uncertain, midpoint values (0.5s)
                                y[min(index_05):length(y)] <- 0.5
                              }

                              return(y)
}

test <- calc(s, fun = impose_unidir)
values(stack(s, test))

beginCluster()
test <- clusterR(s,
                 fun = calc,
                 args = list(impose_unidir)
)

endCluster()
values(stack(s, test))

b <- rep(c(0, 0.5, 0,0, 0.5, 0,0,0), 2)
impose_unidir(b)
r1 <- matrix(b, 4, 4)
r1 <- raster(r1)

b0 <- rep(0, 16)
impose_unidir(b0)



# use overlay with unstack = TRUE








### Function description
# The function has for inputs annual layers of lucfp events at the pixel level.
# It aggregates these pixels to a parcel size defined by parcel_size (in meters).
# The aggregation operation is the sum of the pixel lucfp events.
# Each annual aggregation is tasked in parallel.

  
  ## sequence over which to execute the task.
  # We attribute the tasks to CPU "workers" at the annual level and not at the pf_type level.
  # Hence, if a worker is done with its annual task before the others it can move on to the next one and workers' labor is maximized wrt.
  # attributing tasks at the pf_type level.
  years <- seq(from = 2001, to = 2018, by = 1)
  
  ## read the input to the task
  # is done within each task because it is each time different here.
  
  ## define the task
  annual_aggregate <- function(time){
    # Define which process (island, pf_type, and year) we are in:
    processname <- file.path(paste0("temp_data/processed_lu/annual_maps/lucpfip_",island,"_",pf_type,"_", years[time],".tif"))
    
    #set temp directory
    dir.create(paste0(processname,"_Tmp"), showWarnings = FALSE)
    rasterOptions(tmpdir=file.path(paste0(processname,"_Tmp")))
    
    # read in the input.
    lucpfip_annual <- raster(processname)
    
    # define output file name
    output_filename <- file.path(paste0("temp_data/processed_lu/annual_maps/parcel_lucpfip_",island,"_",parcel_size/1000,"km_",pf_type,"_",years[time],".tif"))
    
    # aggregate it from the ~30m cells to parcel_size cells with mean function.
    raster::aggregate(lucpfip_annual, fact = c(parcel_size/res(lucpfip_annual)[1], parcel_size/res(lucpfip_annual)[2]),
                      expand = FALSE,
                      fun = sum,
                      na.rm = FALSE, # NA cells are in margins, see the NOTES part. If FALSE, aggregations at margins that use NA 
                      # are discarded because the sum would be spurious as it would count all NA as 0s while it is not necessary the case.
                      filename = output_filename,
                      datatype = "INT4U", # because the sum may go up to ~ 10 000 with parcel_size = 3000,
                      # but to more than 65k with parcel_size = 10000 so INT4U will be necessary;
                      overwrite = TRUE)
    #removes entire temp directory without affecting other running processes (but there should be no temp file now)
    unlink(file.path(paste0(processname,"_Tmp")), recursive = TRUE)
    #unlink(file.path(tmpDir()), recursive = TRUE)
    # return the path to this parcels file
    #return(output_filename)
  }
  
  ## register cluster
  registerDoParallel(cores = ncores)
  
  ##  define foreach object.
  foreach(t = 1:length(years),
          # .combine combine the outputs as a mere character list (by default)
          .inorder = FALSE, # we don't care that the results be combine in the same order they were submitted
          .multicombine = TRUE,
          .export = c("island", "parcel_size"),
          .packages = c("raster", "rgdal")
  ) %dopar% annual_aggregate(time = t)



### Execute the function to compute the RasterBrick object of 18 annual layers for each primary forest type

pf_typeS <- c("intact", "degraded", "total")
for(pf_type in pf_typeS){
  # run the computation, that writes the layers 
  parallel_aggregate(pf_type = pf_type, ncores = detectCores() - 1)
  
  # brick the layers together and write the brick
  rasterlist <- list.files(path = "temp_data/processed_lu/annual_maps", 
                           pattern = paste0("parcel_lucpfip_",island,"_",parcel_size/1000,"km_",pf_type,"_"), 
                           full.names = TRUE) %>% as.list()
  
  parcels_brick <- brick(rasterlist)
  
  writeRaster(parcels_brick,
              filename = file.path(paste0("temp_data/processed_lu/parcel_lucpfip_",island,"_",parcel_size/1000,"km_",pf_type,".tif")),
              datatype = "INT4U",
              overwrite = TRUE)
  
  rm(rasterlist, parcels_brick)
  removeTmpFiles(h=0)
}

rasterOptions(tmpdir = "temp_data/raster_tmp")    





#### MAKE FOREST LOSS IN SOUTH EAST ASIA AOI #### 
# import annual layers of forest loss (in hectares) as computed in GEE (and downloaded from Google Drive to input_data/)
