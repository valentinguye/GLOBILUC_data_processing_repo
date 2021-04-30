
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

# dataset_names <- c("glass_gaez_long_country_nf",
#                    "firstloss8320_gaez_long_country_nf",
#                    "phtfloss_gaez_long_country_nf")

# ~6h for the first one, the two second over night
for(name in dataset_names){
  
  path <- paste0(here(origindir, name), ".Rdata")
  df <- readRDS(path)
  
  # Remove gaez variables
  df <- dplyr::select(df,-all_of(gaez_crops))
  
  # Use cross section only
  df_cs <- df[!duplicated(df$grid_id),]
  
  rm(df)
  
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

  # names(df_cs)[names(df_cs) == "OBJECTID"] <- "country_id"
  names(df_cs)[names(df_cs) == "COUNTRY_NA"] <- "country_name"

  df_cs <- st_drop_geometry(df_cs)
  
  # Keep only new variable and id
  df_cs <- df_cs[,c("grid_id", "country_name")]
  
  # # We save the cross section, not the panel, as it is not necessary 
  # # df <- left_join(df, df_cs[,c("grid_id", "country_id", "country_name")], by = "grid_id")

  
  # saveRDS(df_cs, path)
  saveRDS(df_cs, paste0(here(origindir, name), "_country_nf.Rdata"))
  rm(df_cs)
}



# df_cs_mexico <- df_cs[df_cs$country_name == "Brazil","geometry"]
# plot(df_cs_mexico)



 
# df <- readRDS(here(origindir, "glass_gaez_long_country_nf.Rdata"))
# summary(df$soy)


#### STANDARDIZE AND AGGREGATE SUITABILITY INDICES ####  
 
# name = dataset_names[1]

for(name in dataset_names){
  
  path <- paste0(here(origindir, name), ".Rdata")
  df <- readRDS(path)
  # Use cross section only
  df_cs <- df[!duplicated(df$grid_id),]
  
  rm(df)
  
  ### Standardize suitability indexes
  # To understand this line, see https://dplyr.tidyverse.org/articles/rowwise.html#row-wise-summary-functions
  df_cs <- dplyr::mutate(df_cs, si_sum = rowSums(across(.cols = Alfalfa:Yam)))
  
  # df_cs <- dplyr::mutate(df_cs, across(.cols = Alfalfa:Yam, .fns = ~./si_sum))
  # that adds the 47 new variables, with new names.
  df_cs <- dplyr::mutate(df_cs, across(.cols = Alfalfa:Yam, .fns = ~./si_sum, .names = paste0("{.col}", "_std")))
  
  ### Aggregate suitability indexes of similar crops
  # The grouping here corresponds to the categories by GAEZ-IIASA 
  # sugar crops and oil crops could alternatively be categorized as bioenergy feedstock, and Miscanthus etc. as fodder crops (according to Wikipedia).
  # Moreover, we take the maxima of non-standardized SIs as well, for robustness checks that would imply these. 
  df_cs <- df_cs %>% rowwise() %>% mutate(cereal_crops = max(c(Barley, Buckweat, Dryland_rice, Foxtailmillet, Maize, Oat, Pearlmillet, Rye, Sorghum, Wetland_rice, Wheat)), 
                                          oil_crops = max(c(Groundnut, Jatropha, Oilpalm, Olive, Rapeseed, Soybean, Sunflower)),
                                          sugar_crops = max(c(Sugarbeet, Sugarcane)),# Especially necessary to match the price of sugar
                                           fruit_crops = max(c(Banana, Citrus, Cocoa, Coconut)), 
                                           fibre_crops = max(c(Cotton, Flax)),
                                           stimulant_crops = max(c(Coffee, Tea, Tobacco)),
                                           fodder_crops = max(c(Alfalfa)), # To adjust once we have Grass as a GAEZ crop  
                                           bioenergy_crops = max(c(Miscanthus, Reedcanarygrass, Switchgrass)), 
                                          rice_crops = max(c(Dryland_rice, Wetland_rice)) # (this one is not a GAEZ group)
  ) %>% as.data.frame()
  
  # standardize crop groups
  df_cs <- dplyr::mutate(df_cs, across(.cols = contains("_crops"), .fns = ~./si_sum, .names = paste0("{.col}", "_std")))

  # Select only newly constructed variables, and id
  var_names <- names(df_cs)[grepl(pattern = "_std", x = names(df_cs)) | grepl(pattern = "_crops", x = names(df_cs))]
  df_cs <- df_cs[,c("grid_id", var_names)]
  
  # # merge back to panel - NOPE, not anymore
  # df <- left_join(df, df_cs, by = "grid_id")
  
  # Keep only new variables and id
  # df_cs <- df_cs[,(names(df_cs)=="grid_id" | grepl(pattern = "_std", x = names(df_cs)))]
  
  saveRDS(df_cs, paste0(here(origindir, name), "_stdsi.Rdata"))  
  rm(df_cs, path)
}




#### MERGE DATASETS WITH ADDED VARIABLES #### 
for(name in dataset_names){
  
  # Base dataset (including outcome variable(s))
  base_path <- paste0(here(origindir, name), ".Rdata")
  df_base <- readRDS(base_path)
  
  # Remove non-standardized suitability indexes
  df_base <- dplyr::select(df_base,-all_of(gaez_crops))
  
  # Country variable
  country_path <- paste0(here(origindir, name), "_country_nf.Rdata")
  df_country <- readRDS(country_path)
  
  # Merge them and remove to save memory 
  final <- left_join(df_base, df_country, by = "grid_id")
  rm(df_base, df_country)

  # Standardized and aggregated suitability indexes
  stdsi_path <- paste0(here(origindir, name), "_stdsi.Rdata")
  df_stdsi <- readRDS(stdsi_path)  
  
  final <- left_join(final, df_stdsi, by = "grid_id")
  rm(df_stdsi)

  # Create country trends variable
  final <- mutate(final, country_year = paste0(country_name, "_", year))

  saveRDS(final, paste0(here(origindir, name), "_final.Rdata"))
  
  rm(final)
}



### Repeat it for the heavier phtf loss data set, with a slight difference: we do not 
# append to the non standardized suitability indexes (which are only useful for robustness checks)









  