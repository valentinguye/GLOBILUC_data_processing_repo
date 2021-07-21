
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


#### Define GAEZ AESI variables ####
gaez_crops <- list.files(path = here("temp_data", "GAEZ", "v4", "AES_index_value", "Rain-fed", "High-input"), 
                         pattern = "", 
                         full.names = FALSE)
gaez_crops <- gsub(pattern = ".tif", replacement = "", x = gaez_crops)


origindir <- here("temp_data", "merged_datasets", "tropical_aoi")

dataset_names <- c("glass_aesi_long",
                   "firstloss8320_aesi_long", 
                   "phtfloss_aesi_long", 
                   "driverloss_aesi_long")

#### ADD COUNTRY INFORMATION #### 

countries <- st_read(here("input_data", "Global_LSIB_Polygons_Detailed"))
length(unique(countries$COUNTRY_NA)) == nrow(countries)

# #load prices to make country names match
# prices <- readRDS(here("temp_data", "prepared_producer_prices.Rdata"))
# fao_names <- unique(prices$country_name)

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

# dataset_names <- c("glass_aesi_long_country_nf",
#                    "firstloss8320_aesi_long_country_nf",
#                    "phtfloss_aesi_long_country_nf")

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

  
  ### Match country names to those from FAOSTAT 
  # Necessary to match national producer prices 
  # c_names <- unique(df_cs$country_name)
  # c_names[!c_names%in%fao_names]
  
  df_cs$country_name[df_cs$country_name=="United States"] <- "United States of America"
  df_cs$country_name[df_cs$country_name=="Iran"] <- "Iran (Islamic Republic of)"
  df_cs$country_name[df_cs$country_name=="Spain [Canary Is]"] <- "Spain"
  df_cs$country_name[df_cs$country_name=="Burma"] <- "Myanmar"
  df_cs$country_name[df_cs$country_name=="Bahamas, The"] <- "Bahamas"
  # df_cs$country_name[df_cs$country_name=="Taiwan"] 
  df_cs$country_name[df_cs$country_name=="Vietnam"] <- "Viet Nam"
  df_cs$country_name[df_cs$country_name=="Hong Kong (Ch)"] <- "China, Hong Kong SAR"
  df_cs$country_name[df_cs$country_name=="Macau (Ch)"] <- "China, Macao SAR"
  df_cs$country_name[df_cs$country_name=="Turks & Caicos Is (UK)"] <- "United Kingdom of Great Britain and Northern Ireland"
  df_cs$country_name[df_cs$country_name=="Laos"] <- "Lao People's Democratic Republic"
  df_cs$country_name[df_cs$country_name=="Cayman Is (UK)"] <- "United Kingdom of Great Britain and Northern Ireland"
  df_cs$country_name[df_cs$country_name=="Northern Mariana Is (US)"] <- "United States of America"
  df_cs$country_name[df_cs$country_name=="Puerto Rico (US)"] <- "United States of America"
  df_cs$country_name[df_cs$country_name=="Br Virgin Is (UK)"] <- "United Kingdom of Great Britain and Northern Ireland"
  df_cs$country_name[df_cs$country_name=="US Virgin Is (US)"] <- "United States of America"
  df_cs$country_name[df_cs$country_name=="Antigua & Barbuda"] <- "Antigua and Barbuda"
  df_cs$country_name[df_cs$country_name=="St Kitts & Nevis"] <- "Saint Kitts and Nevis"
  df_cs$country_name[df_cs$country_name=="Montserrat (UK)"] <- "United Kingdom of Great Britain and Northern Ireland"
  df_cs$country_name[df_cs$country_name=="Guadeloupe (Fr)"] <- "France"
  df_cs$country_name[df_cs$country_name=="Martinique (Fr)"] <- "France"
  df_cs$country_name[df_cs$country_name=="St Lucia"] <- "Saint Lucia"
  df_cs$country_name[df_cs$country_name=="Guam (US)"] <- "United States of America"
  df_cs$country_name[df_cs$country_name=="St Vincent & the Grenadines"] <- "Saint Vincent and the Grenadines"
  df_cs$country_name[df_cs$country_name=="Venezuela"] <- "Venezuela (Bolivarian Republic of)"
  df_cs$country_name[df_cs$country_name=="Trinidad & Tobago"] <- "Trinidad and Tobago"
  df_cs$country_name[df_cs$country_name=="Cote d'Ivoire"] <- "Côte d'Ivoire"
  df_cs$country_name[df_cs$country_name=="Central African Rep"] <- "Central African Republic"
  # df_cs$country_name[df_cs$country_name=="Micronesia, Fed States of"]
  df_cs$country_name[df_cs$country_name=="French Guiana (Fr)"] <- "France"
  df_cs$country_name[df_cs$country_name=="Congo, Dem Rep of the"] <- "Democratic Republic of the Congo"
  df_cs$country_name[df_cs$country_name=="Brunei"] <- "Brunei Darussalam"
  df_cs$country_name[df_cs$country_name=="Congo, Rep of the"] <- "Congo"
  df_cs$country_name[df_cs$country_name=="Sao Tome & Principe"] <- "Sao Tome and Principe"
  df_cs$country_name[df_cs$country_name=="Tanzania"] <- "United Republic of Tanzania"
  df_cs$country_name[df_cs$country_name=="Solomon Is"] <- "Solomon Islands"
  df_cs$country_name[df_cs$country_name=="French Polynesia (Fr)"] <- "France"
  df_cs$country_name[df_cs$country_name=="Bolivia"] <- "Bolivia (Plurinational State of)"
  df_cs$country_name[df_cs$country_name=="Christmas I (Aus)"] <- "Australia"
  df_cs$country_name[df_cs$country_name=="Mayotte (Fr)"] <- "France"
  df_cs$country_name[df_cs$country_name=="Wallis & Futuna (Fr)"] <- "France"
  df_cs$country_name[df_cs$country_name=="American Samoa (US)"] <- "United States of America"
  df_cs$country_name[df_cs$country_name=="Niue (NZ)"] <- "New Zealand"
  df_cs$country_name[df_cs$country_name=="New Caledonia (Fr)"] <- "France"
  df_cs$country_name[df_cs$country_name=="Reunion (Fr)"] <- "France"
  df_cs$country_name[df_cs$country_name=="Pitcairn Is (UK)"] <- "United Kingdom of Great Britain and Northern Ireland"
  df_cs$country_name[df_cs$country_name=="Swaziland"] <- "Eswatini"
  
  
  # saveRDS(df_cs, path)
  saveRDS(df_cs, paste0(here(origindir, name), "_country_nf.Rdata"))
  rm(df_cs)
}



# df_cs_mexico <- df_cs[df_cs$country_name == "Brazil","geometry"]
# plot(df_cs_mexico)



 
# df <- readRDS(here(origindir, "glass_aesi_long_country_nf.Rdata"))
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
  
  ## Aggregate suitability indexes to match price data
  df_cs <- df_cs %>% rowwise() %>% mutate(Sugar = max(c(Sugarbeet, Sugarcane)),# Especially necessary to match the international price of sugar
                                          Fodder = max(c(Alfalfa, Napiergrass)), # To adjust once we have Grass as a GAEZ crop  
                                          Rice = max(c(Drylandrice, Wetlandrice)),
                                          Sorghum2 = max(c(Sorghum, Sorghumbiomass)),
                                          bioenergy_crops = max(c(Jatropha, Miscanthus, Reedcanarygrass, Switchgrass)),
                                          cereal_crops = max(c(Buckwheat, Foxtailmillet, Pearlmillet, Rye)),
                                          pulses_crops = max(c(Chickpea, Cowpea, Drypea, Gram, Phaseolousbean, Pigeonpea)),
                                          roots_crops = max(c(Cassava, Sweetpotato, Whitepotato, Yam)),
                                          # oil_crops = max(c(Groundnut, Jatropha, Olive, Rapeseed, Sunflower)), # we have prices for all of them
                                          vegetables_crops = max(c(Cabbage, Carrot, Onion, Tomato)),
                                          # fruits_crops = max(c(Banana, Citrus, Coconut)), # For coconut we have price only for the oil but it does not matter because it's the AESI workstream anyway 
                                          industrial_crops = max(c(Flax)) # Rubber, directly responsible for deforestation, are not mixed with these.
                                          # narcotics_crops = max(c(Cocoa, Coffee, Tea, Tobacco)), # We have prices for all of them
                                          
  ) %>% as.data.frame()
  
  # sugar crops and oil crops could alternatively be categorized as bioenergy feedstock, and Miscanthus etc. as fodder crops (according to Wikipedia).

  ## Standardize crops that can be matched with a price and crop groups
  indivcrops_to_std <- c("Banana", "Barley", "Citrus", "Cocoa", "Coconut", "Coffee", "Cotton", "Fodder", 
                         "Groundnut", "Maizegrain", "Oat", "Oilpalm", "Olive", "Rapeseed", "Rice", "Rubber", 
                         "Sorghum2", "Soybean", "Sugar", "Sunflower", "Tea", "Tobacco", "Wheat")
  
  # To understand these lines, see https://dplyr.tidyverse.org/articles/rowwise.html#row-wise-summary-functions
  df_cs <- dplyr::mutate(df_cs, si_sum = rowSums(across(.cols = (contains("_crops") | any_of(indivcrops_to_std)))))

  df_cs <- dplyr::mutate(df_cs, across(.cols = (contains("_crops") | any_of(indivcrops_to_std)),
                                       .fns = ~./si_sum, 
                                       .names = paste0("{.col}", "_std")))

  # Select only newly constructed variables, and id
  var_names <- names(df_cs)[grepl(pattern = "_std", x = names(df_cs))] # this is to add if we want non stded grouped crops too: | grepl(pattern = "_crops", x = names(df_cs))
  # and this is if we want the new variables in non std format too | names(df_cs) %in% c("Fodder", "Rice", "Sugar")
  df_cs <- df_cs[,c("grid_id", var_names)]
  
  saveRDS(df_cs, paste0(here(origindir, name), "_stdsi.Rdata"))  
  rm(df_cs, path)
}




#### MERGE AESI DATASETS WITH ADDED VARIABLES #### 
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





#### Define GAEZ AEAY variables ####
gaez_crops <- list.files(path = here("temp_data", "GAEZ", "v4", "AEAY_out_density", "Rain-fed", "High-input"), 
                         pattern = "", 
                         full.names = FALSE)
gaez_crops <- gsub(pattern = ".tif", replacement = "", x = gaez_crops)


origindir <- here("temp_data", "merged_datasets", "tropical_aoi")

dataset_names <- c("glass_aeay_long",
                   "firstloss8320_aeay_long", 
                   "phtfloss_aeay_long", 
                   "driverloss_aeay_long")

#### GROUP AEAY CROPS #### 
# name = dataset_names[1]

for(name in dataset_names){
  
  path <- paste0(here(origindir, name), ".Rdata")
  df <- readRDS(path)
  # Use cross section only
  df_cs <- df[!duplicated(df$grid_id),]
  
  rm(df)
  
  ## Aggregate suitability indexes to match price data
  df_cs <- df_cs %>% rowwise() %>% mutate(Sugar = max(c(Sugarbeet, Sugarcane)),# Especially necessary to match the international price of sugar
                                          Fodder = max(c(Alfalfa, Napiergrass)),   
                                          Rice = max(c(Drylandrice, Wetlandrice)),
                                          Sorghum2 = max(c(Sorghum, Sorghumbiomass))
                                          # We surely wont need these in the AEAY workflow, as we need to interact these variables with a price and there is no price for those. 
                                          # bioenergy_crops = max(c(Jatropha, Miscanthus, Reedcanarygrass, Sorghumbiomass, Switchgrass)),
                                          # cereal_crops = max(c(Buckwheat, Foxtailmillet, Pearlmillet, Rye)),
                                          # pulses_crops = max(c(Chickpea, Cowpea, Drypea, Gram, Phaseolousbean, Pigeonpea)),
                                          # roots_crops = max(c(Cassava, Sweetpotato, Whitepotato, Yam)),
                                          # # oil_crops = max(c(Groundnut, Jatropha, Olive, Rapeseed, Sunflower)), # we have prices for all of them
                                          # vegetables_crops = max(c(Cabbage, Carrot, Onion, Tomato)),
                                          # # fruits_crops = max(c(Banana, Citrus, Coconut)), # For coconut we have price only for the oil but it does not matter because it's the AESI workstream anyway 
                                          # industrial_crops = max(c(Flax)) # Rubber, directly responsible for deforestation, are not mixed with these.
                                          # # narcotics_crops = max(c(Cocoa, Coffee, Tea, Tobacco)), # We have prices for all of them
  ) %>% as.data.frame()
  
  # Select only newly constructed variables, and id
  df_cs <- df_cs[,c("grid_id", "Fodder", "Rice", "Sugar", "Sorghum2")]
  
  # # merge back to panel - NOPE, not anymore
  # df <- left_join(df, df_cs, by = "grid_id")
  
  # Keep only new variables and id
  # df_cs <- df_cs[,(names(df_cs)=="grid_id" | grepl(pattern = "_std", x = names(df_cs)))]
  
  saveRDS(df_cs, paste0(here(origindir, name), "_groupedcrops.Rdata"))  
  rm(df_cs, path)
}




#### MERGE AEAY DATASETS WITH ADDED VARIABLES #### 
for(name in dataset_names){
  
  # Base dataset (including outcome variable(s))
  base_path <- paste0(here(origindir, name), ".Rdata")
  df_base <- readRDS(base_path)
  
  # Country variable
  # do that because we did not run the spatial join with countries for aeay dataset
  if(name=="glass_aeay_long"){country_path <- paste0(here(origindir, "glass_aesi_long"), "_country_nf.Rdata")}
  if(name=="firstloss8320_aeay_long"){country_path <- paste0(here(origindir, "firstloss8320_aesi_long"), "_country_nf.Rdata")}
  if(name=="phtfloss_aeay_long"){country_path <- paste0(here(origindir, "phtfloss_aesi_long"), "_country_nf.Rdata")}
  
  df_country <- readRDS(country_path)
  
  # Merge them and remove to save memory 
  final <- left_join(df_base, df_country, by = "grid_id")
  rm(df_base, df_country)
  
  # Grouped crops
  grouped_path <- paste0(here(origindir, name), "_groupedcrops.Rdata")
  df_grouped <- readRDS(grouped_path)  
  
  final <- left_join(final, df_grouped, by = "grid_id")
  rm(df_grouped)
  
  # Create country trends variable
  final <- mutate(final, country_year = paste0(country_name, "_", year))
  
  saveRDS(final, paste0(here(origindir, name), "_final.Rdata"))
  
  rm(final)
}




  