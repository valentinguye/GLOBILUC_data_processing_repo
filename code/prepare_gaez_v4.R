
##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages = c("dplyr", "here", "readstata13", 
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
#### PREPARE SUITABILITY INDICES ####
# dir.create(here("temp_data", "GAEZ", "v4", "AES_index_value", "Rain-fed", "High-input"), recursive = TRUE) 
# dir.create(here("temp_data", "GAEZ", "v4", "AES_index_value_current", "Rain-fed", "High-input"), recursive = TRUE) 

datadir <- here("input_data", "GAEZ", "v4", "AES_index_value", "Rain-fed", "High-input")
targetdir <- here("temp_data", "GAEZ", "v4", "AES_index_value", "Rain-fed", "High-input")
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

ext <- extent(c(-180, 179.9167, -30, 30))

## Import most crops
for (j in 1:length(files)) {
  #if (any(crops[j] == cropsToAggregate)) next
  print(files[j])
  # unzip(zipfile = here(datadir, files[j]), exdir = tmpdir)
  dt <- raster(here(datadir, files[j]))
  
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

### CHECK GRASS ### 
# dt <- raster(here("input_data", "GAEZ", "v4", "AEAY_out_density", "Rain-fed", "High-input", "Grass", "data.asc"))
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

writeRaster(gaez_all, here("temp_data", "GAEZ", "v4", "AES_index_value", "Rain-fed", "high_input_all.tif"), 
            overwrite = TRUE)

rm(rasterlist_gaez, gaez_all)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### PREPARE ATTAINABLE YIELDS ####
# dir.create(here("temp_data", "GAEZ", "v4", "AES_index_value", "Rain-fed", "High-input"), recursive = TRUE) 
# dir.create(here("temp_data", "GAEZ", "v4", "AES_index_value_current", "Rain-fed", "High-input"), recursive = TRUE) 

datadir <- here("input_data", "GAEZ", "v4", "AEAY_out_density", "Rain-fed", "High-input")
targetdir <- here("temp_data", "GAEZ", "v4", "AEAY_out_density", "Rain-fed", "High-input")
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

ext <- extent(c(-180, 179.9167, -30, 30))

## Import most crops
for (j in 1:length(files)) {
  #if (any(crops[j] == cropsToAggregate)) next
  print(files[j])
  # unzip(zipfile = here(datadir, files[j]), exdir = tmpdir)
  dt <- raster(here(datadir, files[j]))
  
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

writeRaster(gaez_all, here("temp_data", "GAEZ", "v4", "AEAY_out_density", "Rain-fed", "high_input_all.tif"), 
            overwrite = TRUE)


rm(rasterlist_gaez, gaez_all)
rm(ext, mapmat, mapmat_data)