
### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages = c("plyr", "dplyr", "foreign", "here",
                   "rgdal", "sf", "tmap") #"nngeo"
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

### GLOBAL CRS USED ### 
mercator_world_crs <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "


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

name <- dataset_names[4]

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

#### ADD "CONTINENT" VARIABLE #### 
# We need to set up 3 rectangular shapes. The extreme latitudes are 30 and -30, so we need 2 lines to divide this subtropical band
# We set the line between America and Africa at -27° lon (includes Cape Verde in Africa, and Fernando in Brazil)
# We set the line between Asia and America at -110° lon (includes Eastern Island in South America)
# We set the line between Africa and Asia at 58° lon (includes Maurice in Africa, but also the Arabic peninsula). 

# writing code this way is necessary for asia to actually go from 58° to -110° and not the other way round, if you go with extent %>% bbox etc. 

asia_coords <- matrix(c(58, 30, 180, 30,
                        180, -30, 58, -30, 
                        58, 30), ncol = 2, byrow = TRUE)

asia_ext <- st_polygon(list(asia_coords)) %>% st_sfc(crs = 4326)

america_coords <- matrix(c(-180, 30, -27, 30,
                        -27, -30, -180, -30, 
                        -180, 30), ncol = 2, byrow = TRUE)

america_ext <- st_polygon(list(america_coords)) %>% st_sfc(crs = 4326)

africa_coords <- matrix(c(-27, 30, 58, 30,
                           58, -30, -27, -30, 
                           -27, 30), ncol = 2, byrow = TRUE)

africa_ext <- st_polygon(list(africa_coords)) %>% st_sfc(crs = 4326)

sfc <- c(asia_ext, america_ext, africa_ext)

# sfc <- st_transform(sfc, crs = mercator_world_crs)

continents <- st_sf(data.frame(continent_name = c("Asia", "America", "Africa"), geom = sfc))

# tm_shape(continents)+tm_borders() +tm_fill(col = "continent_name") + tm_graticules() 


for(name in dataset_names){
  
  path <- paste0(here(origindir, name), ".Rdata")
  df <- readRDS(path)
  
  # Remove gaez variables
  df <- dplyr::select(df,-all_of(gaez_crops))
  
  # Use cross section only
  df_cs <- df[!duplicated(df$grid_id),]
  
  # rm(df)
  
  # Spatial
  df_cs <- st_as_sf(df_cs, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  
  
  # This is much much faster (like 5 minutes vs. 6h).
  df_cs <- st_join(x = continents,
                   y = df_cs,
                   join = st_contains,
                   prepared = TRUE,
                   left = FALSE)# performs inner join so returns only records that spatially match.
  
  
  df_cs <- st_drop_geometry(df_cs)
  
  # Keep only new variable and id
  df_cs <- df_cs[,c("grid_id", "continent_name")]
  
  saveRDS(df_cs, paste0(here(origindir, name), "_continent.Rdata"))
  rm(df_cs)

}

# df_cs_mexico <- df_cs[df_cs$country_name == "Brazil","geometry"]
# plot(df_cs_mexico)



 
# df <- readRDS(here(origindir, "glass_aesi_long_country_nf.Rdata"))
# summary(df$soy)


#### REMAINING FOREST ####

df <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "driverloss_aesi_long.Rdata"))

# Remove gaez variables
df <- dplyr::select(df,-all_of(gaez_crops))

year_list <- list()

# in the first year (2001), the past year accumulated deforestation is null. 
year_list[["2001"]] <- df[df$year == 2001, c("grid_id", "year")] 
year_list[["2001"]][,"accu_defo_since2k"] <- 0

# then, each year's deforestation accumulated in the past is the sum of *past years'* deforestation
years <- 2002:max(df$year)
for(y in years){
  sub_ <- df[df$year < y,]
  year_list[[as.character(y)]] <- ddply(sub_, "grid_id", summarise,
                                        accu_defo_since2k = sum(driven_loss, na.rm = TRUE))
  year_list[[as.character(y)]][,"year"] <- y
}

accu_defo_df <- bind_rows(year_list)

df <- inner_join(df, accu_defo_df, by = c("grid_id", "year"))


# summary(df$accu_lucpfp_since2k)
df <- dplyr::mutate(df, 
                     remaining_fc = fc_2000 - accu_defo_since2k)

# df[df$grid_id == 1267,c("grid_id", "year", "lucpfap_pixelcount", "accu_lucpfp_since2k", "remain_pf_pixelcount")] 

# put keep only new variables in remaining
remaining <- df[,c("grid_id", "year", "remaining_fc", "accu_defo_since2k")]

saveRDS(remaining, here("temp_data", "merged_datasets", "tropical_aoi", "driverloss_aesi_long_remaining.Rdata"))


#### STANDARDIZE AND AGGREGATE SUITABILITY INDICES ####  
 
# name = dataset_names[1]
# Nota: produces warnings ("In max(x[x < 0]) : aucun argument pour max ; -Inf est renvoyé")
# this is not an issue. It comes from obs. that have suitability indexes equal across all crops.
# but in such cases, the -Inf value is not given, it is the max_si value that is given (see below) 

# Note also that the standardizing procedures produce some NaN values for stdandardized variables, due to dividing by 0 
# (for grid cells) that have only 0 SI for all crops. This is not an issue, and it will be handled within the data cleaning process in 
# make_main_reg function in analyses.R

for(name in dataset_names){
  
  path <- paste0(here(origindir, name), ".Rdata")
  df <- readRDS(path)
  # Use cross section only
  df_cs <- df[!duplicated(df$grid_id),]
  
  rm(df)
  
  ### Standardize suitability indexes
  
  ## Aggregate suitability indexes to match price data (for Sugar, Fodder, Rice and Sorghum) 
  # and also aggregate crops with common use, including those that would match a price (eg. oat and cotton) bc we do not match SI with prices
  df_cs <- df_cs %>% rowwise(grid_id) %>% mutate(Sugar = max(c(Sugarbeet, Sugarcane)),# Especially necessary to match the international price of sugar
                                                  Fodder = max(c(Alfalfa, Napiergrass)), # To adjust once we have Grass as a GAEZ crop  
                                                  Rice = max(c(Drylandrice, Wetlandrice)),
                                                  Sorghum2 = max(c(Sorghum, Sorghumbiomass)),
                                                  bioenergy_crops = max(c(Jatropha, Miscanthus, Reedcanarygrass, Switchgrass)),
                                                  cereal_crops = max(c(Buckwheat, Foxtailmillet, Pearlmillet, Rye, Oat)),
                                                  pulses_crops = max(c(Chickpea, Cowpea, Drypea, Gram, Phaseolousbean, Pigeonpea, Groundnut)),
                                                  roots_crops = max(c(Cassava, Sweetpotato, Whitepotato, Yam)),
                                                  # oil_crops = max(c(Groundnut, Jatropha, Olive, Rapeseed, Sunflower)), # we have prices for all of them
                                                  vegetables_crops = max(c(Cabbage, Carrot, Onion, Tomato)),
                                                  # fruits_crops = max(c(Banana, Citrus, Coconut)), # For coconut we have price only for the oil but it does not matter because it's the AESI workstream anyway 
                                                  industrial_crops = max(c(Flax, Cotton)) # Rubber, directly responsible for deforestation, are not mixed with these.
                                                  # narcotics_crops = max(c(Cocoa, Coffee, Tea, Tobacco)), # We have prices for all of them
                                                  
  ) %>% as.data.frame()
  
  # sugar crops and oil crops could alternatively be categorized as bioenergy feedstock, and Miscanthus etc. as fodder crops (according to Wikipedia).
  
  ## No. Rather: standardize crops to compute agro-ecological exposure (/risk) to deforestation for a specific LU
  # indivcrops_to_std <- c("Cocoa", "Coffee", "Fodder", "Oilpalm", "Rubber", "Soybean")
  
  # we need a slightly modified version of gaez_crops, that does not include elements of aggregated crops, but includes aggregated crops
  # (necessary for Fodder, Sugar, Rice, Sorghum2, but not for the _crops variables actually, in the current version that does not 
  # account for the number of crops having equal suitable indexes. 
  all_crops <- gaez_crops[!(gaez_crops %in% c("Alfalfa", "Napiergrass", "Sugarbeet", "Sugarcane", 
                                              "Drylandrice", "Wetlandrice", "Sorghum", "Sorghumbiomass", 
                                              "Jatropha", "Miscanthus", "Reedcanarygrass", "Switchgrass", 
                                              "Buckwheat", "Foxtailmillet", "Pearlmillet", "Rye", 
                                              "Chickpea", "Cowpea", "Drypea", "Gram", "Phaseolousbean", "Pigeonpea", 
                                              "Cassava", "Sweetpotato", "Whitepotato", "Yam", 
                                              "Cabbage", "Carrot", "Onion", "Tomato", 
                                              "Flax"))]
  
  all_crops <- c(all_crops, "Fodder", "Sugar", "Rice", "Sorghum2")#"bioenergy_crops", "cereal_crops", "pulses_crops", "roots_crops", "vegetables_crops", "industrial_crops"
  
  ## First way to standardize, sum (the denominator) over all crops
  # To understand these lines, see https://dplyr.tidyverse.org/articles/rowwise.html#row-wise-summary-functions
  df_cs <- dplyr::mutate(df_cs, si_sum = rowSums(across(.cols = (any_of(all_crops)))))#contains("_crops") |

  df_cs <- dplyr::mutate(df_cs, across(.cols = (any_of(all_crops)),#contains("_crops") | 
                                       .fns = ~./(si_sum), 
                                       .names = paste0("{.col}", "_std")))


  ## Second way to standardize: for each crop that can be matched with a price, 
  # standardize by dividing by the sum of the suitability indexes of the N (N = 1,2) crops with the highest suitability (among all, not only among the six drivers), 
  # and give a 0 value to the crops that are not in the top N suitability index. 
  # if N = 1, this procedure is equivalent to sj = 1[Sj = max(Si)]
  
  # identify the highest suitability index values (in every grid cell)
  df_cs <- df_cs %>% rowwise(grid_id) %>% dplyr::mutate(max_si = max(c_across(cols = any_of(all_crops)))) %>% as.data.frame()

  ## N = 1
  # if N = 1, this procedure is equivalent to sj = 1[Sj = max(Si)]
  df_cs <- dplyr::mutate(df_cs, across(.cols = (any_of(all_crops)),
                                       .fns = ~if_else(.==max_si, true = 1, false = 0), 
                                       .names = paste0("{.col}", "_ismax")))

  # and then standardize by the number of different crops being the highest
  all_crops_ismax <- paste0(all_crops,"_ismax")
  df_cs <- dplyr::mutate(df_cs, n_max = rowSums(across(.cols = (any_of(all_crops_ismax)))))

  df_cs <- dplyr::mutate(df_cs, across(.cols = (any_of(all_crops_ismax)),
                                       .fns = ~./n_max, 
                                       .names = paste0("{.col}", "_std1")))
  
  # remove those columns
  df_cs <- dplyr::select(df_cs, !ends_with("_ismax"))  
  
  # rename new ones 
  names(df_cs)[grepl("_std1", names(df_cs))] <- paste0(all_crops, "_std1")
  
  # _std do not sum up to 1, because the _crops vars weight in the denominator but are not counted in the sum. 
  # _std1 do sum up to 1. 
  # df_cs[876,paste0(indivcrops_to_std, "_std")]%>%sum()
  
  ## N = 2: 2nd highest:
  # work on a separate dataset, because the following modifies the base SI data (needed to make SI of non-top2 crops equal to 0)
  working_df_cs <- df_cs
  
  # loop is not efficient but don't know how to code 2nd highest with dplyr
  working_df_cs$max_si_2nd <- NA
  for(i in 1:nrow(working_df_cs)){
    # vector of interest
    x <- working_df_cs[i,all_crops]
    # max value of the row
    row_max_si <- working_df_cs[i,"max_si"]
    
    # second max value
    # we want to include 2 highest values even if there are more than one crop that have the max value and the 2nd max value, hence the lines below are commented out
    working_df_cs[i,"max_si_2nd"] <- max(x[x<row_max_si])
    # if we wanted to handle cases with more than one max values
    #   if(length(which(x == row_max_si))==1){# i.e. if there is only one crop with the highest value
    #     df_cs[i,"max_si_2nd"] <-  row_max_si_2nd
    #   }else{
    #     df_cs[i,"max_si_2nd"] <- row_max_si
    #   }
    
    # /!\ modify SI values (to 0) for all the crops that are not in the top 2 - this is equivalent to weighting by 1 if crop is in top 2, and by 0 if not. 
    working_df_cs[i,names(x[which(x < working_df_cs[i,"max_si_2nd"])])] <- 0
  }
  rm(row_max_si, x)
  
  # then we are able to sum over all crops with weights = 1 if crop is in top 2, 0 if not. 
  working_df_cs <- dplyr::mutate(working_df_cs, si_sum_top2wgted = rowSums(across(.cols = (any_of(all_crops)))))
  # and standardize 
  working_df_cs <- dplyr::mutate(working_df_cs, across(.cols = (any_of(all_crops)),
                                       .fns = ~./(si_sum_top2wgted), 
                                       .names = paste0("{.col}", "_std2")))
  
  # w <- working_df_cs
  # df_cs[875,all_crops]
  # df_cs[875, c(indivcrops_to_std,"si_sum")]
  # df_cs[875, c(paste0(indivcrops_to_std,"_std"),"si_sum")]
  # w[875, c(all_crops,"max_si", "max_si_2nd", "si_sum_top2wgted")]
  # w[875, paste0(all_crops, "_std1")] %>% sum()
  # working_df_cs <- dplyr::mutate(working_df_cs, sum_2_max_si = rowSums(across(.cols = c(max_si, max_si_2nd))))#contains("_crops") |
  


  # working_df_cs <- dplyr::mutate(working_df_cs, across(.cols = (any_of(all_crops)),#contains("_crops") |
  #                                      .fns = ~./(sum_2_max_si),
  #                                      .names = paste0("{.col}", "_std2")))

  # Add standardized variables to the initial data
  std2_var_names <- names(working_df_cs)[grepl(pattern = "_std2", x = names(working_df_cs))]
  df_cs <- inner_join(df_cs, working_df_cs[,c("grid_id",std2_var_names)], by = "grid_id")
  
  ## Crops we need a standardized version of SI: those that can be matched with a price
  # indivcrops_to_std <- c("Banana", "Barley", "Citrus", "Cocoa", "Coconut", "Coffee", "Cotton", "Fodder",
  #                        "Groundnut", "Maizegrain", "Oat", "Oilpalm", "Olive", "Rapeseed", "Rice", "Rubber",
  #                        "Sorghum2", "Soybean", "Sugar", "Sunflower", "Tea", "Tobacco", "Wheat")
  # Adding up the *_std2 over indivcrops_to_std won't equal 1, but this is not an issue. 
  
  # Select variables to save: only newly constructed variables (not grouped crops), and id
  var_names <- c(paste0(all_crops, "_std"), paste0(all_crops, "_std1"), paste0(all_crops, "_std2"))
  # this is to add if we want non stded grouped crops too: | grepl(pattern = "_crops", x = names(df_cs_saved))
  # and we want the new variables in non std format too 
  df_cs <- df_cs[,c("grid_id", "Fodder", "Rice", "Sugar", "Sorghum2", var_names)]
  
  saveRDS(df_cs, paste0(here(origindir, name), "_stdsi.Rdata"))  
  rm(df_cs, working_df_cs, path)
}


# df_cs[,c("si_sum", "max_si", "max_si_2nd", var_names)]


#### MERGE AESI DATASETS WITH ADDED VARIABLES #### 
for(name in dataset_names){
  
  # Base dataset (including outcome variable(s))
  base_path <- paste0(here(origindir, name), ".Rdata")
  df_base <- readRDS(base_path)
  # length(unique(df_base$lon))
  # length(unique(df_base$lat))
  # 
  # df_base$lon <- round(df_base$lon, 6)
  # df_base$lat <- round(df_base$lat, 6)
  # length(unique(df_base$lon))
  # length(unique(df_base$lat))
  # df_base$lat <- as.character(df_base$lat)
  # df_base$lon <- as.character(df_base$lon)
  # 
  # commodity_driven_path <- paste0(here(origindir, "driverloss_commodity_aesi_long.Rdata"))
  # commodity_driven <- readRDS(commodity_driven_path)
  # length(unique(commodity_driven$lon))
  # length(unique(commodity_driven$lat))
  # commodity_driven$lon <- round(commodity_driven$lon, 6)
  # commodity_driven$lat <- round(commodity_driven$lat, 6)
  # length(unique(commodity_driven$lon))
  # length(unique(commodity_driven$lat))
  # commodity_driven$lat <- as.character(commodity_driven$lat)
  # commodity_driven$lon <- as.character(commodity_driven$lon)
  
  # final <- inner_join(df_base, commodity_driven, by = c("lon", "lat", "year"))
 
  # # Remove non-standardized suitability indexes
  # df_base <- dplyr::select(df_base,-all_of(gaez_crops))
  
  # Country variable
  country_path <- paste0(here(origindir, name), "_country_nf.Rdata")
  df_country <- readRDS(country_path)
  # Merge them and remove to save memory 
  final <- left_join(df_base, df_country, by = "grid_id")
  rm(df_base, df_country)

  # Continent variable
  continent_path <- paste0(here(origindir, name), "_continent.Rdata")
  df_continent <- readRDS(continent_path)

  final <- left_join(final, df_continent, by = "grid_id")
  rm(df_continent)
  
  # Standardized and aggregated suitability indexes
  stdsi_path <- paste0(here(origindir, name), "_stdsi.Rdata")
  df_stdsi <- readRDS(stdsi_path)  
  
  final <- left_join(final, df_stdsi, by = "grid_id")
  rm(df_stdsi)

  # Pasture 2000 share of area and remaining variables
  if(name == "driverloss_aesi_long"){
    # pastures
    df_pasture <- readRDS(here("temp_data", "processed_pasture2000", "tropical_aoi", "pasture_4_driverloss_df.Rdata")) 
    names(df_pasture)[names(df_pasture)=="driverloss_masked_pasture"] <- "pasture_share"
    
    final <- left_join(final, df_pasture, by = c("lon", "lat")) # /!\ THIS ACTUALLY DOES NOT MATCH PERFECTLY. HANDLE IF WE REALLY WANT TO USE PASTURE SHARES
    rm(df_pasture)
    
    # remaining
    df_remain <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "driverloss_aesi_long_remaining.Rdata"))
    final <- inner_join(final, df_remain, by = c("grid_id", "year")) # no issue with using grid_id as a key here, bc df_remain was computed just above from the df_base data
    rm(df_remain)
  }
  
  # Create country year fixed effect
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
                                          # # fruits_crops = max(c(Banana, Citrus, Coconut)), # For coconut we have price only for the oil while it is not expressed in oil in GAEZ. 
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
  if(name=="driverloss_aeay_long"){country_path <- paste0(here(origindir, "driverloss_aesi_long"), "_country_nf.Rdata")}
  
  df_country <- readRDS(country_path)
  
  # Merge them and remove to save memory 
  final <- left_join(df_base, df_country, by = "grid_id")
  rm(df_base, df_country)
  
  # Continent variable
  if(name=="glass_aeay_long"){continent_path <- paste0(here(origindir, "glass_aesi_long"), "_continent.Rdata")}
  if(name=="firstloss8320_aeay_long"){continent_path <- paste0(here(origindir, "firstloss8320_aesi_long"), "_continent.Rdata")}
  if(name=="phtfloss_aeay_long"){continent_path <- paste0(here(origindir, "phtfloss_aesi_long"), "_continent.Rdata")}
  if(name=="driverloss_aeay_long"){continent_path <- paste0(here(origindir, "driverloss_aesi_long"), "_continent.Rdata")}
  
  df_continent <- readRDS(continent_path)
  
  final <- left_join(final, df_continent, by = "grid_id")
  rm(df_continent)
  
  # Grouped crops
  grouped_path <- paste0(here(origindir, name), "_groupedcrops.Rdata")
  df_grouped <- readRDS(grouped_path)  
  
  final <- left_join(final, df_grouped, by = "grid_id")
  rm(df_grouped)
  
  # Pasture 2000 share of area variable
  if(name == "driverloss_aeay_long"){
    df_pasture <- readRDS(here("temp_data", "processed_pasture2000", "tropical_aoi", "pasture_4_driverloss_df.Rdata")) 
    names(df_pasture)[names(df_pasture)=="driverloss_masked_pasture"] <- "pasture_share"
    
    final <- left_join(final, df_pasture, by = c("lon", "lat"))
    rm(df_pasture)
    
    # remaining
    # was run only on aesi, hence aesi in the name, but works for aeay too. 
    df_remain <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "driverloss_aesi_long_remaining.Rdata"))
    final <- left_join(final, df_remain, by = c("grid_id", "year"))
    rm(df_remain)
  }
  
  # Create country trends variable
  final <- mutate(final, country_year = paste0(country_name, "_", year))
  
  saveRDS(final, paste0(here(origindir, name), "_final.Rdata"))
  
  rm(final)
}




  