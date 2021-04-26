
### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages = c("plyr", "dplyr", "foreign", "here",
                   "rgdal", "sf") #"nngeo"
#install.packages("sf", source = TRUE)
# library(sf)
# 
# neededPackages = c("tidyverse","data.table", "readxl","foreign", "data.table", "readstata13", "here",
#                    "rgdal", "raster", "velox","sp", "lwgeom", "rnaturalearth", 
#                    "rlist", "parallel", "foreach", "iterators", "doParallel" )
# 

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

# 3. If the troubling packages could not be loaded ("there is no package called ...) 
#   you should try to install them, preferably in their versions stated in the renv.lock file. 
#   see in particular https://rstudio.github.io/renv/articles/renv.html 


# Define GAEZ variables
gaez_crops <- list.files(path = here("temp_data", "GAEZ", "AES_index_value", "Rain-fed", "High-input"), 
                         pattern = "", 
                         full.names = FALSE)
gaez_crops <- gsub(pattern = ".tif", replacement = "", x = gaez_crops)


origindir <- here("temp_data", "merged_datasets", "tropical_aoi")

dataset_names <- c("glass_gaez_long",
                    "firstloss8320_gaez_long", 
                     "phtfloss_gaez_long")

#### ADD COUNTRY INFORMATION #### 

countries <- st_read(here("input_data", "Global_LSIB_Polygons_Detailed"))
length(unique(countries$COUNTRY_NA)) == nrow(countries)

# # path = dataset_paths[1]
# name = dataset_names[1]
# 
# c1 <- countries[1,]
# c2 <- countries[2,]
# plot(c1[,"COUNTRY_NA"])
# plot(c2[,"COUNTRY_NA"])#, add = TRUE
# point <- st_point(c(28.348235, 9.824566), dim = "XY") %>% st_sfc(crs = 4326) %>% st_sf()
# point$col <- 420
# plot(point, col = "red", add = TRUE)
# test <- st_contains(c1, point, sparse = TRUE)
# class(point)
# 
# st_join(x = c1, y = point, join = st_contains, prepared = TRUE, left = FALSE)


# ~6h for the first one, the two second over night
for(name in dataset_names){
  
  path <- paste0(here(origindir, name), ".Rdata")
  df <- readRDS(path)
  
  # Remove gaez variables
  df <- dplyr::select(df,-gaez_crops)
  
  # Use cross section only
  df_cs <- df[!duplicated(df$grid_id),]
  # Spatial
  df_cs <- st_as_sf(df_cs, coords = c("lon", "lat"), crs = 4326, remove = FALSE)


  # This is much much faster (like 5 minutes vs. 6h).
  # df_cs <- st_join(x = countries[,c("OBJECTID", "COUNTRY_NA")],
  #                  y = df_cs,
  #                  join = st_contains,
  #                  prepared = TRUE,
  #                  left = FALSE)# performs inner join so returns only records that spatially match.
  #

  # However, we use st_nearest_feature so that all points match a country
  df_cs <- st_join(x = df_cs,
                   y = countries,
                   join = st_nearest_feature,
                   left = TRUE)

  names(df_cs)[names(df_cs) == "OBJECTID"] <- "country_id"
  names(df_cs)[names(df_cs) == "COUNTRY_NA"] <- "country_name"

  df_cs <- st_drop_geometry(df_cs)

  df <- left_join(df, df_cs[,c("grid_id", "country_id", "country_name")], by = "grid_id")

  # Create country trends variable
  df$country_year <- paste0(df$country_name, "_", df$year)
  

  saveRDS(df, paste0(here(origindir, name), "_country_nf.Rdata"))  
  rm(df, df_cs)
}



# df_cs_mexico <- df_cs[df_cs$country_name == "Brazil","geometry"]
# plot(df_cs_mexico)



 
# df <- readRDS(here(origindir, "glass_gaez_long_country_nf.Rdata"))
# summary(df$soy)


#### STANDARDIZE SUITABILITY INDICES ####  
 
# name = dataset_names[1]

for(name in dataset_names){
  
  path <- paste0(here(origindir, name), ".Rdata")
  df <- readRDS(path)
  # Use cross section only
  df_cs <- df[!duplicated(df$grid_id),]
  
  # To understand this line, see https://www.tidyverse.org/blog/2020/04/dplyr-1-0-0-colwise/ 
  df_cs <- dplyr::mutate(df_cs, si_sum = rowSums(across(.cols = Alfalfa:Yam)))
  
  # df_cs <- dplyr::mutate(df_cs, across(.cols = Alfalfa:Yam, .fns = ~./si_sum))
  # that adds the 47 new variables, with new names.
  df_cs <- dplyr::mutate(df_cs, across(.cols = Alfalfa:Yam, .fns = ~./si_sum, .names = paste0("{.col}", "_std")))
  # Select only newly constructed variables, and id
  df_cs <- df_cs[,(names(df_cs)=="grid_id" | grepl(pattern = "_std", x = names(df_cs)))]
  
  # merge back to panel
  df <- left_join(df, df_cs, by = "grid_id")
  
  # Remove gaez variables (not useful to save them)
  df <- dplyr::select(df,-gaez_crops)
  
  
  saveRDS(df, paste0(here(origindir, name), "_stdsi.Rdata"))  
  rm(df, df_cs)
}
  