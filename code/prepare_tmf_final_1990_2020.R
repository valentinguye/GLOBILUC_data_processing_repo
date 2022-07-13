##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages <- c("data.table", "plyr", "tidyr", "dplyr",  "Hmisc", "sjmisc", "stringr",
                    "here", "readstata13", "foreign", "readxl", "writexl",
                    "raster", "rgdal", "sp", "spdep", "sf","gfcanalysis",  "nngeo", "stars", # "osrm", "osrmr",
                    "lubridate","exactextractr",
                    "doParallel", "foreach", "snow", "parallel",
                    "knitr", "kableExtra",
                    "DataCombine", 
                    "fixest", 
                    "boot", "fwildclusterboot", "sandwich", "MASS",
                    "ggplot2", "leaflet", "tmap",  "dotwhisker", "viridis", "hrbrthemes")
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
dir.create(here("temp_data", "processed_tmf", "tmf_aoi"), recursive = TRUE)
dir.create(here("temp_data", "processed_pasture2000", "tmf_aoi"), recursive = TRUE)
dir.create(here("temp_data", "merged_datasets", "tmf_aoi"), recursive = TRUE)

### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")

### GLOBAL CRS USED ### 
mercator_world_crs <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "

### COUNTRIES ###
countries <- st_read(here("input_data", "Global_LSIB_Polygons_Detailed"))
length(unique(countries$COUNTRY_NA)) == nrow(countries)

### GAEZ OBJECTS
# in this script, GAEZ is the target raster of all aggregations / resamplings
gaez_dir <- here("temp_data", "GAEZ", "v4", "AEAY_out_density",  "Rain-fed")
gaez_crops <- list.files(path = here(gaez_dir, "High-input"), 
                         pattern = "", 
                         full.names = FALSE)
gaez_crops <- gsub(pattern = ".tif", replacement = "", x = gaez_crops)

#### GAEZ #### 
# Athough GAEZ data is kind of in its native resolution and extent, it's data are from a model and we prefer to resample THEM, 
# rather than resampling TMF to GAEZ. This is motivated by resampled tmf rasters yielding systematically lower aggregated statistics than 
# just-aggregate ones. The latter fit correctly Vancutsem et al. statistics. 

# BUT, from an econometric point of view, it's better to have the noise in the LHS than in the RHS, because the former yields only noise, and not attenuation bias. 

## THIS IS GAEZ IN FULL TROPICAL AOI 
# a priori no issue if it's gaez in global aoi, i.e. not croped in prepare_gaez.R, it only needs to be a larger aoi than continental ones given above. 
# besides, note that we brick the file that was already saved as a single brick of raster layers. 
# Otherwise, calling brick on multiple layers takes some time, and calling stack on multiple layers implies that the object is kept in R memory, and it's ~.06Gb
gaez <- brick(here(gaez_dir, "high_input_all.tif"))

# restore gaez crop names 
names(gaez) <- gaez_crops


### TMF DATA ###
# So out of GEE, TMF data is in 9 stacks, 3 continents * 3 drivers (plantation, flood, agri), and each stack has 1990-2020 layers. 
# Resolution is 3km; 
# Unit is hectares of Tropical Moist Forest (TMF) extent, or deforested and replaced with either agriculture/urbanism (agri), tree plantation (plantation), or water (flood)
# there are no NAs. 

# vt <- values(tmf_agri_am[[4]])
# summary(vt)

### TMF AOI ### 
# The procedure adopted in this script currently is to make all preparation by continental area of interest 

ext_am <- stack(here("input_data", "TMF", "TMF_defo_agri_3km_America.tif"))%>% extent()
ext_af <- stack(here("input_data", "TMF", "TMF_defo_agri_3km_Africa.tif"))%>% extent()
ext_as <- stack(here("input_data", "TMF", "TMF_defo_agri_3km_Asia.tif"))%>% extent()

ext_list <- list(America = ext_am, 
                 Africa = ext_af, 
                 Asia = ext_as)
ext_list

# the smallest xmin (most western point) is that of America
# the highest xmax (most eastern point) is that of Asia
# the smallest ymin (most southern point) is the same for all (~ -30°)
# the highest ymax (most northern point) is the same for all (~ 30°)
# tmf_aoi <- extent(ext_am[1], ext_as[2], ext_af[3], ext_af[4])
# this is not the same AOI as tropical_aoi used in other scripts. 

### YEARS ### 
# ---------------------------------- For this project, we need ALL years (1990-2020) --------------------------------------------------------------------
#----------------------------------- THIS IS THE MAIN DIFFERENCE WITH prepare_tmf_final_2010_2020.R script ---------------------------------------------- 
t0 <- 1990
tT <- 2020
wanted_TMF_layers <- paste0("Dec",c(t0:tT))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


# Store continental dataframes in long format in list 
long_df_list <- list()

# type <- "agri"
# CNT <- "America"
transition_types <- c("agricommo", "plantationcommo", "flood")
continents <- c("America", "Africa", "Asia")#  

for(CNT in continents){ # the order of the loops matter for the mask making operation, see below. 
  
  ### AGGREGATE AND RESAMPLE TO GAEZ RESOLUTION ####
  # Currently, the data are hectares of deforestation in 3x3km grid cells annually. 
  # We aggregate this to ~9km grid cells by adding up the hectares of deforestation
  
  # First, crop gaez to the TMF-continental aoi 
  croped_gaez <- crop(gaez, extent(ext_list[[CNT]]))
  
  for(type in transition_types){
    # input stack 
    tmf <- brick(here("input_data", "TMF", paste0("TMF_defo_",type,"_3km_",CNT,".tif")))
    # plot(tmf[[15]])
    # vtmf <- values(tmf[[4]])
    # summary(vtmf)
    
    # Select only years used in this project
    tmf <- tmf[[wanted_TMF_layers]]
    
    # define output file name
    aggr_output_name <- here("temp_data", "processed_tmf", "tmf_aoi", paste0("tmf_",type,"_9km_",CNT,"_",t0,"_",tT,".tif"))
    
    # aggregate it from the 3km cells to 9km
    raster::aggregate(tmf, fact = 3, # c(res(croped_gaez)[1]/res(tmf)[1], res(croped_gaez)[2]/res(tmf)[2]),
                      expand = TRUE,
                      fun = sum,
                      na.rm = TRUE, # NA values are only in the sea. Where there is no forest loss, like in a city in Brazil, the value is 0 (see with plot())
                      filename = aggr_output_name,
                      # datatype = "INT2U", # let the data be float, as we have decimals in the amount of hectares. 
                      overwrite = TRUE)
    
    tmf <- brick(aggr_output_name)
    
    resampled_output_name <- here("temp_data", "processed_tmf", "tmf_aoi", paste0("resampledgaez_tmf_",type,"_9km_",CNT,"_",t0,"_",tT,".tif"))
    
    resample(x = tmf, 
             y = croped_gaez, 
             method = "ngb", # we use ngb and not bilinear because the output values' summary better fits that of the aggregated layer 
             filename = resampled_output_name, 
             overwrite = TRUE)
    
    rm(aggr_output_name, resampled_output_name)
  } # closes the loop over deforestation types
  
  # output nameS of the above loop (this is a length-3 vector)
  resampled_ouput_nameS <- here("temp_data", "processed_tmf", "tmf_aoi",  paste0("resampledgaez_tmf_",
                                                                                  transition_types,"_9km_",CNT,"_",t0,"_",tT,".tif"))
  names(resampled_ouput_nameS) <- transition_types
  
  #### PREPARE MASK BASED ON ALWAYS ZERO FOR EVERY TYPE ####
  # this is done outside the loop over transition types, as it needs to work on all of them at once
  # by checking whether there is a deforestation event related to ANY of the three types.

  # the order of the stacked layers does not matter here
  # and it is necessary to use stack here, since there are several sources.
  any_type <- stack(resampled_ouput_nameS)
  # names(any_type)

  always_zero <- function(y){if_else(condition = (sum(y)==0), true = 0, false = 1)}

  mask_output_name <- here("temp_data", "processed_tmf", "tmf_aoi", paste0("always_zero_mask_any_tmf_type_",CNT,"_",t0,"_",tT,".tif"))

  overlay(x = any_type,
          fun = always_zero,
          filename = mask_output_name,
          na.rm = TRUE, # but there is no NA anyway
          overwrite = TRUE)

  rm(any_type)

  
  #### TMF ANNUAL EXTENTS #### 
  # --> NEED TO MAKE IT IN GEE - COMPUTE TMF EXTENT FOR ALL YEARS SO DON'T HAVE TO MAKE REMAINING HERE 
  tmfext <- brick(here("input_data", "TMF", paste0("TMF_extent_3km_",CNT,".tif")))
  
  # Select only years used in this project
  tmfext <- tmfext[[wanted_TMF_layers]]
  
  tmfext_aggr_output_name <- here("temp_data", "processed_tmf", "tmf_aoi", paste0("tmf_extent_9km_",CNT,"_",t0,"_",tT,".tif"))
  # aggregate it from the 3km cells to ~9km
  raster::aggregate(tmfext, fact = 3, # c(res(croped_gaez)[1]/res(tmfext)[1], res(croped_gaez)[2]/res(tmfext)[2]),
                    expand = TRUE,
                    fun = sum,
                    na.rm = TRUE, # NA values are only in the sea. Where there is no forest loss, like in a city in Brazil, the value is 0 (see with plot())
                    filename = tmfext_aggr_output_name,
                    # datatype = "INT2U", # let the data be float, as we have decimals in the amount of hectares. 
                    overwrite = TRUE)
  
  # # Keep one layer as a destination raster to resample other data to. 
  # tmf_meta <- raster(tmfext_aggr_output_name)
  tmfext <- brick(tmfext_aggr_output_name)
  
  tmfext_resampled_output_name <- here("temp_data", "processed_tmf", "tmf_aoi", paste0("resampledgaez_tmf_extent_9km_",CNT,"_",t0,"_",tT,".tif"))
  
  resample(x = tmfext, 
           y = croped_gaez, 
           method = "ngb", # we use ngb and not bilinear because the output values' summary better fits that of the aggregated layer 
           filename = tmfext_resampled_output_name, 
           overwrite = TRUE)
  
  rm(tmfext_aggr_output_name)
  
  
  #### 2000 PASTURE SHARE ####
  pst2k <- raster(here("input_data", "CroplandPastureArea2000_Geotiff", "Pasture2000_5m.tif"))
  croped_pst2k <- crop(pst2k, extent(ext_list[[CNT]]))
  
  # resample directly (without aggregating first), because pasture data already at the same resolution as croped_gaez. 
  # (aggregate does not work bc resolutions are to close)
  
  # transform NAs to 0, such that the more numerous NAs in pasture data do not force losing information when the stack is turned to a data frame. 
  croped_pst2k <- reclassify(croped_pst2k, cbind(NA, 0))
  
  pst2k_resampled_output_name <- here("temp_data", "processed_pasture2000", "tmf_aoi", paste0("resampledgaez_pasture_2000_",CNT,".tif"))
  
  resample(x = croped_pst2k, 
           y = croped_gaez, 
           method = "ngb", # we use ngb and not bilinear because the output values' summary better fits that of the aggregated layer 
           # and the bilinear interpolation arguably smoothes the reprojection more than necessary given that from and to are already very similar.  
           filename = pst2k_resampled_output_name, 
           overwrite = TRUE)
  
  rm(croped_pst2k, pst2k)
  
  
  #### STACK RASTERS TO MERGE ####
  # Read layers to be stacked
  agri <- brick(resampled_ouput_nameS["agricommo"])
  plantation <- brick(resampled_ouput_nameS["plantationcommo"])
  flood <- brick(resampled_ouput_nameS["flood"])
  tmfext <- brick(tmfext_resampled_output_name)
  #gaez <- brick(gaez_resampled_output_name) # we use croped_gaez here
  pst2k <- raster(pst2k_resampled_output_name)
  
  # It is important to explicitly rename layers that are going to be stacked and then called to reshape the data frame 
  # for time varying variables, the dot is important. 
  names(agri) <- paste0("tmf_agri.",seq(t0, tT, 1)) 
  names(plantation) <- paste0("tmf_plantation.",seq(t0, tT, 1)) 
  names(flood) <- paste0("tmf_flood.",seq(t0, tT, 1)) 
  names(tmfext) <- paste0("tmf_ext.",seq(t0, tT, 1)) 
  
  #names(gaez) <- gaez_crops
  names(pst2k) <- "pasture_share_2000"
  
  # Stack together the annual layers TMF data and GAEZ crop and pasture2000 cross sections 
  continental_stack <- stack(agri,
                             plantation, 
                             flood, 
                             tmfext, 
                             croped_gaez, # NOTICE THIS  
                             pst2k)
  # save those names 
  continental_stack_names <- names(continental_stack)
  
  ### MASK THE STACK TO REMOVE ALWAYS ZERO PIXELS AND LIGHTEN THE DATA FRAMES ###
  # Nope, we don't do that in 1990-2020, in order to have the most general data frame as possible. 
  mask <- raster(mask_output_name)

  masked_stack_output_name <- here("temp_data", "merged_datasets", "tmf_aoi", paste0("anytype_masked_stack_",CNT,"_",t0,"_",tT,".tif"))

  mask(x = continental_stack,
       mask = mask,
       maskvalue = 0, # necessary here, because the there is no NA in the mask, only 0 and 1
       updatevalue = NA,
       filename = masked_stack_output_name,
       overwrite = TRUE)

  rm(mask)
  # (note that masking changes the summary/aggregated values)
  
  
  ### RASTER TO DATAFRAME ### 
  
  continental_stack <- stack(masked_stack_output_name)
  
  # all names are lost through the writing/reading operation, so rename layers 
  names(continental_stack) <- continental_stack_names
  
  # na.rm = TRUE is key here, as it removes previously masked pixels (NA) and ensures the output is not too large (memory intensive)
  # We also set long to false because we reshape with a proper function for more control
  wide_df <- raster::as.data.frame(continental_stack, na.rm = TRUE, xy = TRUE, centroids = TRUE, long = FALSE) 
  
  rm(continental_stack)
  
  # Rename coordinate variables
  names(wide_df)
  head(wide_df[,c("x", "y")])
  wide_df <- dplyr::rename(wide_df, lon = x, lat = y)
  
  
  ### WIDE TO LONG ### 
  
  # grid ID cannot simply be a sequence, because they wouldn't be unique once continental dataframes row-binded 
  # thus add a continental prefix
  wide_df$grid_id <- paste0(CNT, "_", seq(1, nrow(wide_df), 1) )
  
  # create also the continent variable now
  wide_df$continent_name <- CNT
  
  # the dot is, by construction of all variable names, only in the names of time varying variables. 
  # fixed = TRUE is necessary (otherwise the dot is read as a regexp I guess)
  # Note also that it is important that it is structured in a LIST when there are several varying variables in the *long* format
  # Because: "Notice that the order of variables in varying is like x.1,y.1,x.2,y.2."
  varying_vars <- list(names(agri),
                       names(plantation),
                       names(flood),
                       names(tmfext))
  
  # reshape to long.
  long_df <- stats::reshape(wide_df,
                            varying = varying_vars,
                            v.names = c("tmf_agri", "tmf_plantation", "tmf_flood", "tmf_ext"),
                            sep = ".",
                            timevar = "year",
                            idvar = "grid_id", # don't put "lon" and "lat" in there, otherwise memory issue (see https://r.789695.n4.nabble.com/reshape-makes-R-run-out-of-memory-PR-14121-td955889.html)
                            ids = "grid_id", # lonlat is our cross-sectional identifier.
                            direction = "long",
                            new.row.names = NULL)#seq(from = 1, to = nrow(ibs_msk_df)*length(years), by = 1)
  
  names(long_df)
  # replace the indices from the raster::as.data.frame with actual years.
  
  years <- seq(t0, tT, 1)
  long_df <- mutate(long_df, year = years[year])
  
  long_df <- dplyr::arrange(long_df, grid_id, year)
  
  # d <- long_df[long_df$driven_loss != long_df$driven_loss_commodity,]
  # d <- mutate(d, diff = driven_loss - driven_loss_commodity)
  # summary(d$diff)
  # d[d$diff>0 , c("driven_loss", "driven_loss_commodity")]
  
  long_df_list[[CNT]] <- long_df
  
  rm(agri, plantation, flood, tmfext, croped_gaez, pst2k)
  
  rm(resampled_ouput_nameS, tmfext_aggr_output_name, tmfext_resampled_output_name, pst2k_resampled_output_name, 
     mask_output_name, masked_stack_output_name, wide_df, long_df) 
} # closes the loop over continents


tropical_long_df <- bind_rows(long_df_list)

saveRDS(tropical_long_df, here("temp_data", "merged_datasets", "tmf_aoi", paste0("tmf_gaezresampled_aeay_pantrop_long_",t0,"_",tT,".Rdata")))

rm(tropical_long_df, long_df_list)

#### COUNTRY VARIABLE #### 
countries <- st_read(here("input_data", "Global_LSIB_Polygons_Detailed"))
length(unique(countries$COUNTRY_NA)) == nrow(countries)

path <- here("temp_data", "merged_datasets", "tmf_aoi",  paste0("tmf_aeay_pantrop_long_",t0,"_",tT,".Rdata"))
df <- readRDS(path)

# Remove gaez variables
df <- dplyr::select(df,-all_of(gaez_crops))

# Use cross section only
df_cs <- df[!duplicated(df$grid_id),]
# Spatial
df_cs <- st_as_sf(df_cs, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# nrow(df_cs)*length(unique(df$year)) == nrow(df)
rm(df)

# transform shapes such that st_nearest_feature works with projected data
df_cs <- st_transform(df_cs, crs = mercator_world_crs)
countries <- st_transform(countries, crs = mercator_world_crs)

# This is much much faster (like 5 minutes vs. 6h).
df_cs <- st_join(x = countries[,c("OBJECTID", "COUNTRY_NA")],
                 y = df_cs,
                 join = st_contains,
                 prepared = TRUE,
                 left = FALSE)# performs inner join so returns only records that spatially match.


# However, we use st_nearest_feature so that all points match a country
# df_cs <- st_join(x = df_cs,
#                  y = countries,
#                  join = st_nearest_feature,
#                  left = TRUE)

rm(countries)

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
saveRDS(df_cs, here("temp_data", "merged_datasets", "tmf_aoi", "tmf_pantrop_cs_country_nf.Rdata"))
rm(df_cs)


#### BIGGER CELL VARIABLES #### 
## Prepare base data
path <- here("temp_data", "merged_datasets", "tmf_aoi",  paste0("tmf_aeay_pantrop_long_",t0,"_",tT,".Rdata"))
df <- readRDS(path)

# Remove gaez variables
df <- dplyr::select(df,-all_of(gaez_crops))
# Use cross section only
df_cs <- df[!duplicated(df$grid_id),]
# Spatial
df_cs <- st_as_sf(df_cs, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

rm(df)


## Prepare bigger square grids
# They are based on the resolution of the raster that has been converted to data frame.  
# Available only by continent here 
for(CNT in c("America", "Africa", "Asia")){
  grid_base <- raster(here("temp_data", "merged_datasets", "tmf_aoi", paste0("anytype_masked_stack_",CNT,"_",t0,"_",tT,".tif")))
  
  # for ~45km grid cells (5 times larger grid cells in both dimensions, hence 25 times larger)
  bigger_5 <- aggregate(grid_base, fact = 5, expand = TRUE, fun = sum)
  bigger_5_stars <- st_as_stars(bigger_5)
  bigger_5_sf <- st_as_sf(bigger_5_stars, as_points = FALSE, merge = FALSE, na.rm = FALSE, long = FALSE)
  
  # Make the grid cell index. It is continent-specific, otherwise it won't be unique across continents. 
  bigger_5_sf$grid_id_5 <- paste0(CNT, "_", seq(from = 1, to = nrow(bigger_5_sf)))
  bigger_5_sf <- bigger_5_sf[,c("grid_id_5", "geometry")]
  
  # repeat for ~90km grid cells (10 times larger grid cells in both dimensions, hence 10 times larger)
  bigger_10 <- aggregate(grid_base, fact = 10, expand = TRUE, fun = sum)
  bigger_10_stars <- st_as_stars(bigger_10)
  bigger_10_sf <- st_as_sf(bigger_10_stars, as_points = FALSE, merge = FALSE, na.rm = FALSE, long = FALSE)
  
  # Make the grid cell index. It is continent-specific, otherwise it won't be unique across continents. 
  bigger_10_sf$grid_id_10 <- paste0(CNT, "_", seq(from = 1, to = nrow(bigger_10_sf)))
  bigger_10_sf <- bigger_10_sf[,c("grid_id_10", "geometry")]
  
  # repeat for ~180km grid cells (20 times larger grid cells in both dimensions, hence 20 times larger)
  bigger_20 <- aggregate(grid_base, fact = 20, expand = TRUE, fun = sum)
  bigger_20_stars <- st_as_stars(bigger_20)
  bigger_20_sf <- st_as_sf(bigger_20_stars, as_points = FALSE, merge = FALSE, na.rm = FALSE, long = FALSE)
  
  # Make the grid cell index. It is continent-specific, otherwise it won't be unique across continents. 
  bigger_20_sf$grid_id_20 <- paste0(CNT, "_", seq(from = 1, to = nrow(bigger_20_sf)))
  bigger_20_sf <- bigger_20_sf[,c("grid_id_20", "geometry")]
  
  # DO NOT TRANSFORM, because it makes the inner_join associate two bigger squares for each smaller one, 
  # I don't really know why, but it does not do so with geographic coordinates, and there is no mismatch, and no problem of imprecision due to 
  # geographic coordinates being used as planer because aoi is not near the pole. 
  # df_cs <- st_transform(df_cs, crs = mercator_world_crs)
  # bigger_50km_sf <- st_transform(bigger_50km_sf, crs = mercator_world_crs)
  # bigger_100km_sf <- st_transform(bigger_100km_sf, crs = mercator_world_crs)
  
  # join the variables
  df_cs <- st_join(x = df_cs,
                    y = bigger_5_sf,
                    join = st_within,
                    prepared = FALSE, # tests shows that this changes nothing, whether shapes are transformed or not
                    left = TRUE)
  
  df_cs <- st_join(x = df_cs,
                    y = bigger_10_sf,
                    join = st_within,
                    prepared = FALSE, # tests shows that this changes nothing, whether shapes are transformed or not
                    left = TRUE)
  
  df_cs <- st_join(x = df_cs,
                    y = bigger_20_sf,
                    join = st_within,
                    prepared = FALSE, # tests shows that this changes nothing, whether shapes are transformed or not
                    left = TRUE)
}
# left = FALSE  performs inner join so returns only records that spatially match.
# df_cs1 <- df_cs
df_cs <- st_drop_geometry(df_cs)

df_cs[,grepl("grid_id_", names(df_cs))] <- sapply(df_cs[,grepl("grid_id_", names(df_cs))], function(x) {if_else(is.na(x), "", x)})

# use max, as paste0 does not work...
df_cs <- df_cs %>% rowwise(grid_id) %>% dplyr::mutate(newname_5 = max(c_across(cols = starts_with("grid_id_5"))),
                                                      newname_10 = max(c_across(cols = starts_with("grid_id_10"))),
                                                      newname_20 = max(c_across(cols = starts_with("grid_id_20")))) %>% as.data.frame()

# head(df_cs) 
# Keep only new variable and id
df_cs <- df_cs[,c("grid_id", "newname_5", "newname_10", "newname_20")]
names(df_cs) <- c("grid_id", "grid_id_5", "grid_id_10", "grid_id_20")
# length(unique(df_cs$grid_id_5))
length(unique(df_cs$grid_id_10))
length(unique(df_cs$grid_id_20))


saveRDS(df_cs, here("temp_data", "merged_datasets", "tmf_aoi", "tmf_pantrop_cs_biggercells.Rdata"))

rm(df_cs)



#### GROUP AND STANDARDIZE AEAY CROPS #### 
# all groupings in this section are motivated on the GAEZ v4 model documentation, and in particular Table A4-1.3

df <- readRDS(here("temp_data", "merged_datasets", "tmf_aoi",  paste0("tmf_aeay_pantrop_long_",t0,"_",tT,".Rdata")))
# Use cross section only
df_cs <- df[!duplicated(df$grid_id),]

rm(df)

# NOW, prices are needed only for SOME crops, those that are grouped but GAEZ yield in dry weights are not comparable. 
prices <- readRDS(here("temp_data", "prepared_international_prices.Rdata"))

### This matrix is used for maping crops from GAEZ with commodities from price data sets
mapmat_data <- c(
  "Banana","Banana",
  "Barley", "Barley",
  "Crude_oil", "Biomass", # these crop categories are gonna be created in the present script
  "Orange", "Citrus", # Citrus sinensis in both GAEZ and FAO
  "Cocoa", "Cocoa",
  "Coconut_oil", "Coconut", # Coconut not excluded as we don't use prices anymore. See below, in conversion part, why we would exclude it if we needed price scaling 
  "Coffee", "Coffee",
  "Cotton", "Cotton",
  "Beef", "Fodder", # these crop categories are gonna be created in the present script 
  "Groundnuts", "Groundnut",
  "Maize", "Maizegrain",
  "Oat", "Oat",
  "Olive_oil", "Olive",  
  "Palm_oil", "Oilpalm",
  "Rapeseed_oil", "Rapeseed",
  "Rice", "Rice",
  "Rubber", "Rubber",
  "Sorghum", "Sorghum", # this will be matched with barley and wheat, i.e. grains we have price data on.
  "Soybean", "Soybean",
  "Soybean_meal", "Soybean_meal", # these crop categories are gonna be created in the present script
  "Soybean_oil", "Soybean_oil", # these crop categories are gonna be created in the present script
  "Sugar", "Sugar", # these crop categories are gonna be created in the present script
  "Sugar", "Sugarbeet",
  "Sugar", "Sugarcane",
  "Sunflower_oil", "Sunflower",
  "Tea", "Tea",
  "Tobacco", "Tobacco", 
  "Wheat", "Wheat")

mapmat <- matrix(data = mapmat_data, 
                 nrow = length(mapmat_data)/2,
                 ncol = 2, 
                 byrow = TRUE)

colnames(mapmat) <- c("Prices", "Crops")

# crops to group based on potential REVENUE
crops2grp <- c("Barley", "Sorghum", "Wheat", "Cocoa", "Coffee", "Groundnut", "Rapeseed", "Sunflower")

# crops to standardize. There is not fodder, rubber, citrus, banana, cocoa, coffee, olive and tea
eaear2std <- paste0("eaear_", c("Cereals", "Oilfeed_crops", "Cotton", "Maizegrain", "Oat", "Oilpalm", "Rice",
                                "Soy_compo", "Sugar", "Tobacco")) 
# add cocoa, coffee and tea for std2 
eaear2std_bis <- paste0("eaear_", c("Banana", "Biomass", "Cereals", "Oilfeed_crops", "Cocoa_Coffee", "Cotton", 
                                    "Maizegrain", "Oat", "Olive", "Oilpalm", "Rice",
                                    "Soy_compo", "Sugar", "Tea", "Tobacco")) 

# this vector is used to convert potential production in dry weight from GAEZ into weight of traded commodities. 
# Olive, Oil palm, and Sugar crops are already provided in GAEZ in units of traded goods (oil or sugar)
conv_fac <- c(Wheat = 0.87, 
              Rice = 0.87, 
              Maizegrain = 0.86, 
              Sorghum = 0.87, 
              Barley = 0.87, 
              Oat = 0.87, 
              Soybean = 0.90, 
              Rapeseed = 0.90, 
              Sunflower = 0.92, 
              Groundnut = 0.65, # this is applied to go from shelled groundnuts to GAEZ dry weight 
              # Cotton = 0.33, we won't use this, see below paragraph on cotton
              Banana = 0.25, # from banana to dry weight banana
              Tobacco = 0.75) # from traded tobacco dry leaves to GAEZ dry weight (i.e. traded leaves have water content of 25%)  


# no vlaue missing for any of these commodities in the broad period 1991-2019
# working_prices <- prices[,c("year", mapmat[,"Prices"])]%>%filter(year>1990 & year < 2020)

# get prices for 2000
price_avg <- prices %>% 
  filter(year>=1995 & year <= 2004) %>% 
  #filter(year==2000) %>% 
  summarise(across(.cols = any_of(mapmat[,"Prices"]), 
                   .fns = mean, na.rm = TRUE))


## Aggregate suitability indexes to match price data
# makes sense for SUgar fodder and rice subcrops because yields expressed in comaprable units  
df_cs <- df_cs %>% rowwise() %>% mutate(Biomass = max(c(Miscanthus, Reedcanarygrass, Sorghumbiomass*10, Switchgrass)), # sorghum biomass is expressed in kg/ha, and not 10kg/ha as it is the case for the three other crops
                                        #  don't include Jatropha because it is expressed in seeds and not above ground biomass in GAEZ. 
                                        Fodder = max(c(Alfalfa, Napiergrass)),   
                                        Rice = max(c(Drylandrice, Wetlandrice)),
                                        Sugar = max(c(Sugarbeet, Sugarcane)) # Especially necessary to match the international price of sugar
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

# keep only crops of interest 
df_cs <- df_cs[, names(df_cs) %in% c("grid_id", "lon", "lat", mapmat[,"Crops"])]
## TONS OF WHAT? For some crops, the AEAY quantity does not necessarily match the market price unit  
# "For most crops the agro-climatic potential yield is given as kg dry weight per hectare. 
# For alfalfa, miscanthus, switchgrass, reed canary grass, napier grass, pasture legumes and grasses the yield is given in 10 kg dry weight per hectare. 
# For sugar beet and sugarcane (and hence Sugar, the max of them) yields are in kg sugar per hectare, 
# and for oil palm and olives in kg oil per hectare. Cotton yield is given as kg lint per hectare." 
# https://gaez.fao.org/pages/theme-details-theme-3

# BANANA. 
# Just convert from dry matter weight
df_cs <- dplyr::mutate(df_cs, Banana = Banana / conv_fac["Banana"])

# BARLEY. 
# Just convert from dry matter weight
df_cs <- dplyr::mutate(df_cs, Barley = Barley / conv_fac["Barley"])

# COCONUT. 
# Coconut cannot be converted to coconut oil, as it is made of only a by product of coconut; 
# We would not correctly estimate the value of coconut if we counted the whole coconut yield as the copra byproduct;   

# COTTON.
# From FAO documentation and v4 data visualization, the unit of cotton AEAY is kg lint/ha. 
# cotton lint is "raw ginned cotton which is ready for baling" (baling is 'packing' cotton in standard volumes), see there for definitions http://agropedia.iitk.ac.in/content/glossary-useful-terms-related-cotton
# the price in pink sheet is in $/kg and comes from Cotton Outlook A Index, which is expressed for 'raw cotton', see there https://www.cotlook.com/information-2/the-cotlook-indices-an-explanation/
# which I understand as being raw ginned cotton, as it is given for a particular grade, and ginned (baled) cotton is graded, not pure raw cotton (before ginning)
# THUS, the price is expressed in the same unit as the AEAY
# NOTE that the conversion factor of 0.33 used in GAEZ should not be applied, because it converts from seed cotton to cotton lint to match with FAOSTAT data, which is expressed in seed cotton weight. 

# FODDER. 
# Convert fodder crop yield into beef by a feed conversion ratio. 
# Following Galloway et al. 2007 who get a ratio of 20 (feed to meat conversion rate of 0.05) for ruminants (beef and sheep and goats) on non arable land (i.e. for fodder, not feed from crops which is more efficient)
# This is quite similar to Alexander et al. 2016, 25 tons of feed transform into 1 ton of beef meat (in edible weight). 
# Galloway being more specific about FCR of non-arable land feed, we retain this. 
# Lower (more efficient) FCR found in the literature, typically below 10, represent feed, not fodder (often necessary to compare with non-ruminants)
# Thus, every ton of fodder (coming from agro-climatically achievable yields in ton/ha) is scaled to 1/20 ton of beef meat  
# BUT needs to account for a dry matter content conversion to actually grazed feed. 
# We apply an dry matter content of 0.271 as reported for summer green chop alfalfa in https://www.ccof.org/sites/default/files/Feed%20Type%20DMI%20Table%20Final.pdf
# (which has a similar dry matter content as napier grass)
df_cs <- dplyr::mutate(df_cs, Fodder = Fodder * 0.05 / 0.271)

# GROUNDNUT. 
# Oil content is 45-56% in https://link.springer.com/chapter/10.1007%2F978-94-011-0733-4_6 as reporte by https://link.springer.com/article/10.1007/s11746-017-2981-3
# it is 31% in table 32 in https://www.ers.usda.gov/webdocs/publications/41880/33132_ah697_002.pdf
# Because it's a bit unclear what the value is in the literature, and price is available in PS for groudnuts (not oil), we don't convert to oil.
# Nevertheless, we convert from DM weight kernel (in GAEZ) to shelled weight (Pink Sheet), using GAEZ conversion factor
df_cs <- dplyr::mutate(df_cs, Groundnut = Groundnut / conv_fac["Groundnut"])

# MAIZE. 
# Just convert from dry matter weight
df_cs <- dplyr::mutate(df_cs, Maizegrain = Maizegrain / conv_fac["Maizegrain"])

# OAT
df_cs <- dplyr::mutate(df_cs, Oat = Oat / conv_fac["Oat"])

# OLIVE AND OIL PALM. 
# They are expressed in oil already, as in the price data. Hence nothing to do. 

# RAPESEED.
# Seeds contain around 41% oil see Table 4 in Yasar 2018, and https://www.agmrc.org/commodities-products/grains-oilseeds/rapeseed 
df_cs <- dplyr::mutate(df_cs, Rapeseed = Rapeseed * 0.41 / conv_fac["Rapeseed"])

# RUBBER. 
# "In general, latex contains about 30-40 % of rubber particles and 55-65 % of water. However, fresh latex shows 15-45 % of rubber hydrocarbon and about
# 2-4 % of non-rubber ingredients [2]. Latex is usually sold either in the form of dry rubber sheet or concentrated rubber solution."
# https://www.measurement.sk/2014/Kerdtongmee.pdf
# So the RSS3 product the PS price is for, is Ribbed Smoked Sheet (dry rubber sheet), i.e. a raw form of latex, that is said to have a dry rubber content of 30-40%. 
# Dry rubber content from this page (60%) is for concentrated rubber/latex, not RSS http://www.unistarglobal.com/natural_rubber.php 
# So here 1 ton of dry rubber is diluted to produce a higher weight of RSS. 
df_cs <- dplyr::mutate(df_cs, Rubber = Rubber / 0.35)

# SORGHUM. 
# Just convert from dry matter weight
df_cs <- dplyr::mutate(df_cs, Sorghum = Sorghum / conv_fac["Sorghum"])

# SUNFLOWER
# The oil content is set at 42%, https://www.sciencedirect.com/science/article/pii/B9780123849472006747
# it is 40.5% in table 32 in https://www.ers.usda.gov/webdocs/publications/41880/33132_ah697_002.pdf and https://www.agmrc.org/commodities-products/grains-oilseeds/sunflower-profile 
# in Table 4 in Yasar 2018 it's 40-50%, 
# 22–55% oil content (Flagella et al. 2002; Gonzalez-Martin et al. 2013) as reported in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5976617/
df_cs <- dplyr::mutate(df_cs, Sunflower = Sunflower * 0.42 / conv_fac["Sunflower"])

# SOY
# 18% of the seed is extracted as oil and 82% as meal, so we make two distinct yields 
df_cs <- dplyr::mutate(df_cs, Soybean = Soybean / conv_fac["Soybean"])
df_cs <- dplyr::mutate(df_cs, Soybean_oil = Soybean * 0.18 / conv_fac["Soybean"])
df_cs <- dplyr::mutate(df_cs, Soybean_meal = Soybean * 0.82 / conv_fac["Soybean"])
# Compute also the value for mere Soy beans

# RICE 
# Just convert from dry matter weight
df_cs <- dplyr::mutate(df_cs, Rice = Rice / conv_fac["Rice"])

# TOBACCO. 
# prices are for unmanufactured tobacco 
# Just convert from dry matter weight
df_cs <- dplyr::mutate(df_cs, Tobacco = Tobacco / conv_fac["Tobacco"])

# WHEAT. 
# Just convert from dry matter weight
df_cs <- dplyr::mutate(df_cs, Wheat = Wheat / conv_fac["Wheat"])


## CONVERT FROM kg/ha to ton/ha
# All prices have been converted to $/ton in prepare_prices.R but all yields are expressed in kg/ha
df_cs <- dplyr::mutate(df_cs, across(.cols = all_of(mapmat[,"Crops"]),
                                     .fns = ~./1000)) 
# Convert Nappier grass and alfalfa (components of Fodder), and Biomass from now 10 tons to tons. 
df_cs <- dplyr::mutate(df_cs, across(.cols = all_of(c("Fodder", "Biomass")),
                                     .fns = ~.*10)) 

## Interact with average prices to get Expected Agro-Ecological Attainable Revenue (EAEAR)

# Prices have been converted to $/t in prepare_prices.R
# we actually need to convert yields to revenues only for crops that we want to group
# yet, for sthe sake of generality, we convert all crops to their revenues
# for non-grouped crops, it's not a big deal if the conversion to $/ha is not accurate, as comparisons will only be within crops, between grid cells 
for(aeay_i in mapmat[,"Crops"]){
  price_i <- price_avg[mapmat[mapmat[,"Crops"]==aeay_i,"Prices"]]%>%as.numeric()
  eaear_i <- paste0("eaear_", aeay_i)
  df_cs <- dplyr::mutate(df_cs, 
                         !!as.symbol(eaear_i) := !!as.symbol(aeay_i) * price_i)
}

df_cs <- df_cs %>% rowwise() %>% mutate(eaear_Cereals = max(c(eaear_Barley, eaear_Sorghum, eaear_Wheat)), 
                                        eaear_Oilfeed_crops = max(c(eaear_Groundnut, eaear_Rapeseed, eaear_Sunflower)), 
                                        eaear_Cocoa_Coffee = max(c(eaear_Cocoa, eaear_Coffee))) %>% as.data.frame()

# and for Soy commodities: 
df_cs <- dplyr::mutate(df_cs, eaear_Soy_compo =  eaear_Soybean_meal + eaear_Soybean_oil)

# the highest values for Soy_compo represent the value added from processing into oil/meals
summary(df_cs$eaear_Soybean)
summary(df_cs$eaear_Soy_compo)
# Retain processed soy commodities, because this is more traded, and more comparable to the other oil seeds that are expressed in oil too. 
df_cs <- dplyr::select(df_cs, -eaear_Soybean, -eaear_Soybean_meal, -eaear_Soybean_oil)

### STANDARDIZE ### 

# # divide each revenue variable by the sum of them, to standardize.
# # To understand this line, see https://dplyr.tidyverse.org/articles/rowwise.html#row-wise-summary-functions
# df_cs <- dplyr::mutate(df_cs, eaear_sum = rowSums(across(.cols = any_of(eaear2std))))
# 
# df_cs <- dplyr::mutate(df_cs, across(.cols = any_of(eaear2std),
#                                      .fns = ~./eaear_sum, 
#                                      .names = paste0("{.col}", "_std"))) 
# 
# df_cs <- dplyr::select(df_cs, -eaear_sum)
# 
# ## Second way to standardize: for each crop that can be matched with a price, 
# # standardize by dividing by the sum of the suitability indexes of the N (N = 1,2) crops with the highest suitability (among all, not only among the six drivers), 
# # and give a 0 value to the crops that are not in the top N suitability index. 
# # if N = 1, this procedure is equivalent to sj = 1[Sj = max(Si)]
# 
# # this code is intricate but it handles potential but unlikely cases where two crops have equal EAEAR
# 
# # identify the highest suitability index values (in every grid cell)
# df_cs <- df_cs %>% rowwise(grid_id) %>% dplyr::mutate(max_eaear = max(c_across(cols = any_of(eaear2std)))) %>% as.data.frame()
# 
# # and for the alternative set of crops
# df_cs <- df_cs %>% rowwise(grid_id) %>% dplyr::mutate(max_eaearbis = max(c_across(cols = any_of(eaear2std_bis)))) %>% as.data.frame()
# 
# ## N = 1
# # if N = 1, this procedure is equivalent to sj = 1[Sj = max(Si)]
# df_cs <- dplyr::mutate(df_cs, across(.cols = any_of(eaear2std),
#                                      .fns = ~if_else(.==max_eaear, true = 1, false = 0), 
#                                      .names = paste0("{.col}", "_ismax")))
# 
# # and then standardize by the number of different crops being the highest
# all_crops_ismax <- paste0(eaear2std,"_ismax")
# df_cs <- dplyr::mutate(df_cs, n_max = rowSums(across(.cols = (any_of(all_crops_ismax)))))
# 
# unique(df_cs$n_max) # 16 is when GAEZ is NA
# # df_cs[df_cs$n_max==16,]
# 
# df_cs <- dplyr::mutate(df_cs, across(.cols = (any_of(all_crops_ismax)),
#                                      .fns = ~./n_max, 
#                                      .names = paste0("{.col}", "_std1")))
# 
# # remove those columns
# df_cs <- dplyr::select(df_cs, !ends_with("_ismax"))  
# 
# # rename new ones 
# names(df_cs)[grepl("_std1", names(df_cs))] <- paste0(eaear2std, "_std1")
# 
# # _std do sum up to 1. 
# # df_cs[87687,paste0(eaear2std, "_std2")]%>%sum()
# 
# ## N = 2: 2nd highest:
# # helper function to apply within dplyr rowwise framework
# max2nd <- function(x){max(x[x!=max(x)])}
# 
# df_cs <- df_cs %>% rowwise() %>% mutate(max_eaear_2nd = max2nd(c_across(cols = (any_of(eaear2std))))) %>% as.data.frame()
# # this returns a warning about aucun argument pour max ; -Inf est renvoyé" when there are only zero values for every crop in the grid cell
# # not a problem
# 
# # new column for each crop, telling whether it's in the top 2 or not
# df_cs <- dplyr::mutate(df_cs, across(.cols = any_of(eaear2std),
#                                      .fns = ~if_else(.>=max_eaear_2nd, true = 1, false = 0), 
#                                      .names = paste0("{.col}", "_istop2")))
# 
# all_crops_istop2 <- paste0(eaear2std,"_istop2")
# 
# df_cs <- dplyr::mutate(df_cs, n_top2 = rowSums(across(.cols = (any_of(all_crops_istop2)))))
# unique(df_cs$n_top2) # 16 is when GAEZ is NA
# # df_cs[df_cs$n_max==16,]
# # in other words there is always only 2 crops in the top 2     
# 
# # multiply istop2 columns with their corresponding EAEAR columns
# for(crop in eaear2std){
#   df_cs <- mutate(df_cs, !!as.symbol(paste0(crop, "_istop2")) := !!as.symbol(paste0(crop, "_istop2")) * !!as.symbol(crop) )
# }
# 
# # sum them, and standardize with this sum
# df_cs <- dplyr::mutate(df_cs, sum_top2 = rowSums(across(.cols = (any_of(all_crops_istop2)))))
# 
# df_cs <- dplyr::mutate(df_cs, across(.cols = (any_of(all_crops_istop2)),
#                                      .fns = ~./sum_top2, 
#                                      .names = paste0("{.col}", "_std2")))
# 
# # remove those columns
# df_cs <- dplyr::select(df_cs, !ends_with("_istop2"))  
# # rename columns
# names(df_cs)[grepl("_istop2_std2", names(df_cs))] <- paste0(eaear2std, "_std2")  
# 
# 
# ### ### ### 
# ## Repeat std2 for the broader set of crops 
# df_cs <- df_cs %>% rowwise() %>% mutate(max_eaear_2nd_bis = max2nd(c_across(cols = (any_of(eaear2std_bis))))) %>% as.data.frame()
# 
# # new column for each crop, telling whether it's in the top 2 or not
# df_cs <- dplyr::mutate(df_cs, across(.cols = any_of(eaear2std_bis),
#                                      .fns = ~if_else(.>=max_eaear_2nd_bis, true = 1, false = 0), 
#                                      .names = paste0("{.col}", "_istop2")))
# 
# all_crops_istop2_bis <- paste0(eaear2std_bis,"_istop2")
# 
# # multiply istop2 columns with their corresponding EAEAR columns
# for(crop in eaear2std_bis){
#   df_cs <- mutate(df_cs, !!as.symbol(paste0(crop, "_istop2")) := !!as.symbol(paste0(crop, "_istop2")) * !!as.symbol(crop) )
# }
# 
# # sum them, and standardize with this sum
# df_cs <- dplyr::mutate(df_cs, sum_top2_bis = rowSums(across(.cols = (any_of(all_crops_istop2_bis)))))
# 
# df_cs <- dplyr::mutate(df_cs, across(.cols = (any_of(all_crops_istop2_bis)),
#                                      .fns = ~./sum_top2_bis, 
#                                      .names = paste0("{.col}", "_std2bis")))
# 
# # remove those columns
# df_cs <- dplyr::select(df_cs, !ends_with("_istop2"))  
# # rename columns
# names(df_cs)[grepl("_istop2_std2bis", names(df_cs))] <- paste0(eaear2std_bis, "_std2bis")  
# 
# ### ### ### ### ### 
# 
# # Select variables to save: all the eaear variables, standardized or not
# df_cs <- dplyr::select(df_cs, -max_eaear_2nd, - max_eaear_2nd_bis)
var_names <- grep(pattern = "eaear_", names(df_cs), value = TRUE) 
df_cs <- df_cs[,c("grid_id", var_names)]

saveRDS(df_cs, here("temp_data", "merged_datasets", "tmf_aoi", "tmf_pantrop_cs_stdeaear.Rdata"))  
rm(df_cs, prices, price_avg, mapmat)



#### MERGE ADDITIONAL VARIABLES ####  
# Base dataset (including outcome variable(s))
df_base <- readRDS(here("temp_data", "merged_datasets", "tmf_aoi",  paste0("tmf_aeay_pantrop_long_",t0,"_",tT,".Rdata")))

# Remove unprocesed gaez variables, for memory purpose
df_base <- dplyr::select(df_base,-all_of(gaez_crops))

## EAEAR
df_stdeaear <- readRDS(here("temp_data", "merged_datasets", "tmf_aoi",  "tmf_pantrop_cs_stdeaear.Rdata"))  

final <- left_join(df_base, df_stdeaear, by = "grid_id") # no issue with using grid_id as a key here, bc df_remain was computed just above from the df_base data
rm(df_stdeaear)

## COUNTRY
df_country <- readRDS(here("temp_data", "merged_datasets", "tmf_aoi", "tmf_pantrop_cs_country_nf.Rdata"))

# Merge them and remove to save memory 
final <- left_join(final, df_country, by = "grid_id")
rm(df_country)

# Create country-year identifier
final <- mutate(final, country_year = paste0(country_name, "_", year))

## BIGGER CELLS
df_biggercells <- readRDS(here("temp_data", "merged_datasets", "tmf_aoi", "tmf_pantrop_cs_biggercells.Rdata"))

final <- left_join(final, df_biggercells, by = "grid_id")
rm(df_biggercells)

# Create bigger cell-year identifier
final <- mutate(final, grid_id_5_year = paste0(grid_id_5, "_", year))
final <- mutate(final, grid_id_10_year = paste0(grid_id_10, "_", year))
final <- mutate(final, grid_id_20_year = paste0(grid_id_20, "_", year))
# length(unique(final$grid_id_50km_year))==length(unique(final$grid_id_50km))*length(unique(final$year))


saveRDS(final, here("temp_data", "merged_datasets", "tmf_aoi", paste0("tmf_aeay_pantrop_long_final_",t0,"_",tT,".Rdata")))

rm(final)



# #### COUNTRY ANNUAL AVERAGES #### 
# final <- readRDS(here("temp_data", "merged_datasets", "tmf_aoi", paste0("tmf_aeay_pantrop_long_final_",t0,"_",tT,".Rdata")))
# 
# country_final <- dplyr::select(final, year, continent_name, country_name, country_year, tmf_agri, tmf_plantation, tmf_flood, tmf_ext)
# rm(final)
# country_avges <- ddply(country_final, "country_year", summarise, 
#                        agriurba_defo = sum(tmf_agri, na.rm = TRUE),
#                        plantation_defo = sum(tmf_plantation, na.rm = TRUE),
#                        flood_defo = sum(tmf_flood, na.rm = TRUE),
#                        tmf_extent = sum(tmf_ext, na.rm = TRUE), 
#                        .progress = "text")
# 
# country_year_final <- country_final[!duplicated(country_final$country_year),] 
# 
# country_avges <- inner_join(country_avges, country_year_final[,c("year", "continent_name", "country_name", "country_year")], by = "country_year")
# 
# country_avges <- dplyr::mutate(country_avges, deforestation = agriurba_defo + plantation_defo + flood_defo)
# country_avges <- dplyr::mutate(country_avges, deforestation = deforestation/1e6)
# country_avges <- dplyr::mutate(country_avges, tmf_extent = tmf_extent/1e6)
# 
# 
# 
# # some checks against statistics in Vancutsem Science article
# # In Tab. 10, they report 39.6, 88.7, and 60.9 Mha deforestation over 1990-2020, in Africa, America and Asia resp.
# # (excluding regrowth, but including prior degradation)
# country_avges %>% filter(continent_name=="Africa") %>% dplyr::select(deforestation) %>% sum()
# country_avges %>% filter(continent_name=="America") %>% dplyr::select(deforestation) %>% sum()
# country_avges %>% filter(continent_name=="Asia") %>% dplyr::select(deforestation) %>% sum()
# # here, we obtain 33.8, 83.4, and 56 Mha respectively
# 
# # Using aggregated and resampled rasters instead 
# # Aggregated by fact = 3: 
# # America: 88.4 Mha 
# agri <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("tmf_agri_9km_America_1990_2020.tif")))
# plantation <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("tmf_plantation_9km_America_1990_2020.tif")))
# flood <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("tmf_flood_9km_America_1990_2020.tif")))
# any <- overlay(stack(agri, plantation, flood), 
#                fun = sum)
# vany <- values(any)
# sum(vany)/1e6
# 
# # Asia : 59.8
# agri <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("tmf_agri_9km_Asia_1990_2020.tif")))
# plantation <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("tmf_plantation_9km_Asia_1990_2020.tif")))
# flood <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("tmf_flood_9km_Asia_1990_2020.tif")))
# any <- overlay(stack(agri, plantation, flood), 
#                fun = sum)
# vany <- values(any)
# sum(vany)/1e6
# 
# # Resampled raster: 83.4 Mha too 
# agri <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("resampledgaez_tmf_agri_America_1990_2020.tif")))
# plantation <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("resampledgaez_tmf_plantation_America_1990_2020.tif")))
# flood <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("resampledgaez_tmf_flood_America_1990_2020.tif")))
# any <- overlay(stack(agri, plantation, flood), 
#                fun = sum)
# vany <- values(any)
# sum(vany)/1e6
# 
# # Aggregated raster: 88.47 Mha, so very very similar to Vancutsem et al. 
# agri <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("aggrgaez_tmf_agri_America_1990_2020.tif")))
# plantation <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("aggrgaez_tmf_plantation_America_1990_2020.tif")))
# flood <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("aggrgaez_tmf_flood_America_1990_2020.tif")))
# any <- overlay(stack(agri, plantation, flood), 
#                fun = sum)
# vany <- values(any)
# sum(vany)/1e6
# 
# # So the difference comes from the resampling. 
# 
# # resampling with bilinear: 83.3
# croped_gaez <- crop(gaez, ext_list[["America"]])
# bilany <- resample(any, croped_gaez, method = "bilinear")
# sum(values(bilany))/1e6
# 
# # and ngb (what we used, but layer by layer) : 83.4
# ngbany <- resample(any, croped_gaez, method = "ngb")
# sum(values(ngbany))/1e6
# 
# # so it's really just a matter of resampling, not the method used, or what is done before; 
# # try project raster instead of resample: same results. 
# ngbany <- projectRaster(any, croped_gaez, method = "bilinear")
# sum(values(ngbany))/1e6
# 
# # trying different resampling from terra: does not help. 
# any <- terra::rast(any)
# croped_gaez <- terra::rast(croped_gaez)
# terra::resample(any, croped_gaez, method = "near") %>% values() %>% sum()/1e6 # same as raster::resample
# terra::resample(any, croped_gaez, method = "bilinear") %>% values() %>% sum()/1e6 # same as raster::resample
# terra::resample(any, croped_gaez, method = "cubic") %>% values() %>% sum()/1e6 # 83.279
# terra::resample(any, croped_gaez, method = "cubicspline") %>% values() %>% sum()/1e6 # 83.278
# terra::resample(any, croped_gaez, method = "lanczos") %>% values() %>% sum()/1e6 # 83.281
# 
# 
# 
# 
# # TMF EXTENT 
# # In Tab. 2, they report 273.4, 705.1, 311,1 Mha TMF (undisturbed and degraded) in 1990 in Africa, America and ASia resp. 
# # They do not count regrowth, while I incude it. However, in 1990 this is not observed in their data anyways. 
# # Since we masked places where no deforestation ever occured, we should not use the final data frame to make such comparisons. 
# # Resampled raster: 677.07 Mha
# am_tmfext <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("resampledgaez_tmfext_America_1990_2020.tif")))
# am_tmf90 <- values(am_tmfext[[1]])
# sum(am_tmf90)/1e6
# # Aggregated raster: 719.6 Mha
# am_tmfext <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("aggrgaez_tmfext_America_1990_2020.tif")))
# am_tmf90 <- values(am_tmfext[[1]])
# sum(am_tmf90)/1e6
# 
# # with bilinear resampling, very similar sum: 677.29
# croped_gaez <- crop(gaez, ext_list[["America"]])
# aggr <- am_tmfext[[1]]
# resamp_bil <- resample(aggr, croped_gaez, method = "bilinear") 
# sum(values(resamp_bil))/1e6
# 
# # Resampled Africa: 263.3 Mha
# af_tmfext <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("resampledgaez_tmfext_Africa_1990_2020.tif")))
# af_tmf90 <- values(af_tmfext[[1]])
# sum(af_tmf90)/1e6
# 
# # Resampled Asia: 303.4 Mha
# as_tmfext <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("resampledgaez_tmfext_Asia_1990_2020.tif")))
# as_tmf90 <- values(as_tmfext[[1]])
# sum(as_tmf90)/1e6
# 
# 
# 
# # Old tests --------------
# 
# # they report 238.7, 616.4 and 244.1 Mha of TMF (undisturbed and degraded) in 2015 in Africa, Latin America, and Asia-Oceania resp. 
# country_avges %>% filter(continent_name=="Africa" & year == 1989) %>% dplyr::select(tmf_extent) %>% sum()
# country_avges %>% filter(continent_name=="America" & year == 1989) %>% dplyr::select(tmf_extent) %>% sum()
# country_avges %>% filter(continent_name=="Asia" & year == 1989) %>% dplyr::select(tmf_extent) %>% sum()
#     
# head(country_avges)
# brazil <- dplyr::filter(country_avges, country_name == "Brazil" & year > 2000)
# indonesia <- dplyr::filter(country_avges, country_name == "Indonesia" & year > 2000)
# 
# 
# plot(x = brazil$year, y = brazil$deforestation)
# plot(x = indonesia$year, y = indonesia$deforestation)
# plot(x = brazil$year, y = brazil$tmf_extent)
# plot(x = indonesia$year, y = indonesia$tmf_extent)
# 
# # Compare with FIg. 3 in Vancutsem et al. 2021 Science. Indonesian panel
# # The present measure of deforestation is to be compare with the sum of 
# # - direct defo not followed by regrowth
# # - conversion to plantations
# # - conversion to water
# # - defor after degradatation not followed by regrowth (i.e. degradation -> no regrowth -> defor)
# # defor after degradation followed by regrowth (i.e. i.e. degradation -> regrowth -> defor)
# # degradation before defor
# 
# 
# America <- dplyr::filter(country_avges, continent_name == "America" & year > 1994)
# Asia <- dplyr::filter(country_avges, continent_name == "Asia" & year > 1994)
# America <- ddply(America, "year", summarise, 
#                  deforestation = sum(deforestation))
# Asia <- ddply(Asia, "year", summarise, 
#                  deforestation = sum(deforestation))
# plot(x = America$year, y = America$deforestation)
# plot(x = Asia$year, y = Asia$deforestation)
# 
# 
# # compare with Table 11
# brazil[brazil$year > 2000 & brazil$year<2013,"deforestation"] %>% mean()
# America[America$year > 2000 & America$year<2013,"deforestation"] %>% mean()
# America[America$year > 2000 & America$year<2020,"deforestation"] %>% mean()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  ### ### ### ### ### ### ### ### ### ### ### 
