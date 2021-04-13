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
dir.create(here("temp_data", "processed_phtfloss", "annual_maps"), recursive = TRUE) 

### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")



#### PREPARE PHTF LOSS #### 
  
# import raster brick of phtf loss, annual binary layers
rasterlist <- list.files(path = here("input_data", "phtfLossSumGlass_maxP_crsT_intermediate"), 
                         pattern = "phtfLossSumGlass_maxP_crsT_intermediate_", 
                         full.names = TRUE) %>% as.list()
phtfloss <- brick(rasterlist)


test <- raster(here("input_data", "phtfLossSumGlass_maxP_2002.tif"))
test
plot(test)
summary(values(test))


### Aggregate it to the GAEZ resolution ### 

# Currently, the data are hectares of phtf loss in 5x5km grid cells annually. 
# We aggregate this to 10km grid cells by adding up the hectares of phtf loss

# Read in a GAEZ grid, because we want to align to it. 
gaez <- raster(here("temp_data", "GAEZ", "AES_index_value", "Rain-fed", "High-input", "Alfalfa.tif"))

annual_aggregate_andAllign <- function(year){
  # Read in annual layer input
  input_name <- paste0(here("input_data", "phtfLossSumGlass_maxP", "phtfLossSumGlass_maxP_"), 
                        year,".tif")
  
  phtfloss_annual <- raster(input_name)
  
  
  # define output file name
  aggr_output_name <- paste0(here("temp_data", "processed_phtfloss", "annual_maps", "phtfLossSumGaez_Raggr_"),
                            year,".tif")
  
  # aggregate it from the ~5km cells to ~10km
  raster::aggregate(phtfloss_annual, fact = c(res(gaez)[1]/res(phtfloss_annual)[1], res(gaez)[2]/res(phtfloss_annual)[2]),
                    expand = FALSE,
                    fun = sum,
                    na.rm = TRUE, # NA values are only in the sea. Where there is no phtf loss, like in a city in Brazil, the value is 0 (see with plot())
                    filename = aggr_output_name,
                    # datatype = "INT2U", # let the data be float, as we have decimals in the amount of hectares. 
                    overwrite = TRUE)
  
  
  # Align it to GAEZ exactly
  prepared <- raster(aggr_output_name)
  
  # define output file name
  aligned_output_name <- paste0(here("temp_data", "processed_phtfloss", "annual_maps", "phtfLossSumGaez_Ralign_"),
                                    year,".tif")
  
  resample(x = prepared, 
           y = gaez, 
           method = "ngb", # we use ngb and not bilinear because the output values' summary better fits that of the aggregated layer 
           # and the bilinear interpolation arguably smoothes the reprojection more than necessary given that from and to are already very similar.  
           filename = aligned_output_name, 
           overwrite = TRUE)

  # aligned <- raster(aligned_output_name)
  # 
  # prepared_values <- values(prepared)
  # aligned_values <- values(aligned)
  # 
  # summary(prepared_values[prepared_values != 0 & !is.na(prepared_values)])
  # summary(aligned_values[aligned_values != 0 & !is.na(aligned_values)])
  
}



# Execute the function 
years <- seq(from = 2002, to = 2020, by = 1)

for(t in years){
  annual_aggregate_andAllign(year = t)
}
  

# Brick together the annual layers and with GAEZ crop cross sections 
rasterlist_phtfloss <- list.files(path = here("temp_data", "processed_phtfloss", "annual_maps"), 
                                   pattern = "phtfLossSumGaez_Ralign_", 
                                   full.names = TRUE) %>% as.list()

rasterlist_gaez <- list.files(path = here("temp_data", "GAEZ", "AES_index_value", "Rain-fed", "High-input"), 
                              pattern = "", 
                              full.names = TRUE) %>% as.list()

rasterlist <- append(rasterlist_phtfloss, rasterlist_gaez)

parcels_brick <- brick(rasterlist)






# writeRaster(parcels_brick,
#             filename = file.path(paste0("temp_data/processed_lu/parcel_lucpfip_",island,"_",parcel_size/1000,"km_",pf_type,".tif")),
#             datatype = "INT4U",
#             overwrite = TRUE)








