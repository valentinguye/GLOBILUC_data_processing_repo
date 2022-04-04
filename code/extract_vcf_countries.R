
##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 

neededPackages <- c("data.table", "plyr", "tidyr", "dplyr",  "Hmisc", "sjmisc", "stringr",
                    "here", "readstata13", "foreign", "readxl", "writexl",
                    "raster", "rgdal", "sp", "spdep", "sf","stars", "gfcanalysis",  "nngeo", # "osrm", "osrmr",
                    "lubridate","exactextractr",
                    "doParallel", "foreach", "snow", 
                    "knitr", "kableExtra",
                    "DataCombine", 
                    "fixest", 
                    "boot", "fwildclusterboot", "sandwich",
                    "ggplot2", "leaflet", "tmap", "dotwhisker")


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
dir.create(here("temp_data", "processed_vcf", "no_aoi"), recursive = TRUE)
dir.create(here("temp_data", "processed_vcf", "southam_aoi"), recursive = TRUE)
dir.create(here("temp_data", "processed_vcf", "no_aoi"), recursive = TRUE)


### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


# GLASS GLC codebook 
# 0 = No data, 10 = cropland, 20 = forest, 30 = Grassland, 40 = Shrubland, 70 = Tundra, 90 = Barren land, 100 = Snow/ice

#### WRITE TO BRICK #### 
# Read data in
rasterlist <- list.files(path = "input_data/VCF5KYR", 
                         pattern = paste0("VCF5KYR"), 
                         full.names = TRUE) %>% as.list()

# each of these rasters has 3 bands
for(ras in rasterlist){

  r <- stack(ras)
  # select band of interest 
  r <- r[[1]]

  vr <- values(r)    
}
parcels_brick <- stack(rasterlist)

names(parcels_brick)
# Write
writeRaster(parcels_brick, here("temp_data", "processed_vcf", "no_aoi", "brick_no_aoi.tif"), 
            overwrite = TRUE)


no_aoi <- brick(here("temp_data", "processed_vcf", "no_aoi", "brick_no_aoi.tif"))

#### 83-15 TRACK: FIRST FOREST LOSS ####
# Create annual layers of forest loss defined as: class is not 20 in a given year while it was 20 in *all the previous year*.
# This restricts annual loss to that occurring for the first time

# construct previous: a collection of annual layers, each giving the mean of the class value in the previous years. 
previous_years <- seq(1982, 2014, 1) 
for(t in 1:length(previous_years)){ # goes only up to 2014, as we don't need the average up to 2015.
  calc(no_aoi[[1:t]], fun = mean, 
       filename = here("temp_data", "processed_vcf", "no_aoi", paste0("past_mean_lu_",previous_years[t], ".tif")), 
       datatype = "FLT4S", # necessary so that a 19.9 mean is not counted as a 20 (i.e. so far undisturbed forest pixel)
       overwrite = TRUE)
}

# this is a raster of 33 layers, giving the mean of GLC class value in 1982, 1982-83, 1982-84, ..., 1982-2014. 
rasterlist <- list.files(path = here("temp_data", "processed_vcf", "no_aoi"), 
                         pattern = "^past_mean_lu_", 
                         full.names = TRUE) %>% as.list()
previous <- stack(rasterlist)

#unique(values(previous))

make_first_loss <- function(previous, current){if_else(condition = (previous == 20 & current != 20), 
                                                       true = 1, false = 0)}

years <- seq(1982, 2015, 1) 
for(t in 2:length(years)){ # starts from 1983 as we need t-1 and thus t starts from 2
  overlay(previous[[t-1]], no_aoi[[t]], fun = make_first_loss, 
          filename = here("temp_data", "processed_vcf", "no_aoi", paste0("first_loss_",years[t], ".tif")), 
          datatype = "INT1U", 
          overwrite = TRUE) 
}

# first loss no disturbance: here we count forest loss as deforestation if it is the first time forest loss AND it is not 
# a disturbance in forest cover observation (as told by the next year being not forest cover either) 
make_nd_first_loss <- function(previous, current, sbqt){if_else(condition = (previous == 20 & current != 20 & sbqt != 20), 
                                                                true = 1, false = 0)}

years <- seq(1982, 2014, 1) 
for(t in 2:length(years)){ # starts from 1983 as we need t-1 and thus t starts from 2. And ends in 2014 as we need 2015 for sbqt. 
  overlay(previous[[t-1]], no_aoi[[t]], no_aoi[[t+1]], fun = make_nd_first_loss, 
          filename = here("temp_data", "processed_vcf", "no_aoi", paste0("nd_first_loss_",years[t], ".tif")), # nd for no disturbance
          datatype = "INT1U", 
          overwrite = TRUE) 
}

#### 83-15 TRACK: CONVERT TO AREA #### 
# so now, after make_nd_first_loss, the values are binary
make_area <- function(values, areas){
  # *100 to convert km2 (returned by raster::area) to hectares, to match phtfl 
  return(values*areas*100) 
}

## SIMPLE FIRST LOSS 
rasterlist <- list.files(path = here("temp_data", "processed_vcf", "no_aoi"),
                         pattern = paste0("^first_loss_"), 
                         full.names = TRUE) %>% as.list()
first_loss <- stack(rasterlist)

cell_area <- area(first_loss)

overlay(first_loss, cell_area, 
        fun = make_area, 
        filename = here("temp_data", "processed_vcf", "no_aoi", "ha_first_loss.tif"), 
        overwrite = TRUE)


## NO DISTURBANCE 
rasterlist <- list.files(path = here("temp_data", "processed_vcf", "no_aoi"),
                         pattern = paste0("^nd_first_loss_"), 
                         full.names = TRUE) %>% as.list()
nd_first_loss <- stack(rasterlist)

cell_area <- area(nd_first_loss)

overlay(nd_first_loss, cell_area, 
        fun = make_area, 
        filename = here("temp_data", "processed_vcf", "no_aoi", "ha_nd_first_loss.tif"), 
        overwrite = TRUE)

removeTmpFiles(h = 0)


#### EXTRACT IN COUNTRIES #### 
# at this point, values in pixels are in hectares of the pixels, if the pixel has first loss. 
first_loss_ha <- brick(here("temp_data", "processed_vcf", "no_aoi", "ha_first_loss.tif"))

## Prepare countries
countries <- st_read(here("input_data", "Global_LSIB_Polygons_Detailed"))
length(unique(countries$COUNTRY_NA)) == nrow(countries)
# countries <- st_transform(countries, crs = mercator_world_crs)
countries_sp <- as(countries, "Spatial")

# all.equal(crs(countries), crs(first_loss_ha))

# two ways to go from here, depending on what's fastest 
# either raster::extract, or as.data.frame first_loss_ha, removing pixels in the sea, and st_joint-st_contains. 

df_cs <- raster::extract(first_loss_ha, 
                         y = countries_sp, 
                         fun = sum, 
                         na.rm = TRUE,
                         exact = FALSE,
                         nl = nlayers(first_loss_ha), 
                         df = TRUE, 
                         sp = TRUE)

saveRDS(df_cs, here("temp_data", "processed_glass_glc", "country_extract.Rdata"))


df_cs_sf <- st_as_sf(df_cs)
df_cs_df <- st_drop_geometry(df_cs_sf)

write.csv(df_cs, here("temp_data", "processed_glass_glc", "country_extract.csv"), row.names=FALSE, quote=FALSE)


nd_first_loss_ha <- brick(here("temp_data", "processed_vcf", "no_aoi", "ha_nd_first_loss.tif"))

# two ways to go from here, depending on what's fastest 
# either raster::extract, or as.data.frame first_loss_ha, removing pixels in the sea, and st_joint-st_contains. 

df_cs_nd <- raster::extract(nd_first_loss_ha, 
                            y = countries_sp, 
                            fun = sum, 
                            na.rm = TRUE,
                            exact = FALSE,
                            nl = nlayers(first_loss_ha), 
                            df = TRUE, 
                            sp = TRUE) 
# df_cs <- st_join(x = countries[,c("OBJECTID", "COUNTRY_NA")],
#                  y = df_cs,
#                  join = st_contains,
#                  prepared = TRUE,
#                  left = FALSE)# performs inner join so returns only records that spatially match.





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

