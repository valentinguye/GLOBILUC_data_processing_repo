##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY STARTS WHERE THIS SCRIPT IS LOCATED (see package "here")
# input data is searched in a folder called "./input_data/" 
# output data is stored in a folder called "./temp_data/" 

### PACKAGES ###

# These are the packages needed in this particular script. 
neededPackages <- c("data.table", "plyr", "tidyr", "dplyr", "writexl",
                    "here",
                    "raster", "rgdal", "sp", "sf")

# Load them 
lapply(neededPackages, library, character.only = TRUE)


### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")

### COUNTRIES
countries <- st_read(here("input_data", "Global_LSIB_Polygons_Detailed"))
length(unique(countries$COUNTRY_NA)) == nrow(countries)

### TMF DATA
# TMF deforestation data is in 9 stacks: 3 continents * 3 drivers (plantation, flood, agri), and each stack has 1990-2020 layers. 
# Preparation in GEE, at: https://code.earthengine.google.com/

# TMF extent is in 3 stacks (for each continent), each with 1990-2020 layers.
# Preparation in GEE, at: https://code.earthengine.google.com/

# For all rasters: 
# Resolution is 3km; 
# Unit is hectares of Tropical Moist Forest (TMF) extent, or deforested and replaced with either agriculture/urbanism (agri), tree plantation (plantation), or water (flood)
# there are no NAs. 


### YEARS
# ---------------------------------- Here, need ALL years (1990-2020) --------------------------------------------------------------------
t0 <- 1990
tT <- 2020
wanted_TMF_layers <- paste0("Dec",c(t0:tT))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


# Store continental dataframes in long format in list 
long_df_list <- list()

# type <- "agri"
# CNT <- "America"
transition_types <- c("agri", "plantation", "flood")
continents <- c("America", "Africa", "Asia")#  

for(CNT in continents){ # the order of the loops matter for the mask making operation, see below. 
  
  ### AGGREGATE AND alignE TO GAEZ RESOLUTION ####
  # Currently, the data are hectares of deforestation in 3x3km grid cells annually. 
  # We aggregate this to ~9km grid cells by adding up the hectares of deforestation
  for(type in transition_types){
    # input 
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

    rm(aggr_output_name)

  } # closes the loop over deforestation types
  
  # output nameS of the above loop (this is a length-3 vector)
  aggregated_ouput_nameS <- here("temp_data", "processed_tmf", "tmf_aoi",  paste0("tmf_",transition_types,"_9km_",CNT,"_",t0,"_",tT,".tif"))
  names(aggregated_ouput_nameS) <- transition_types
  
  
  #### TMF ANNUAL EXTENTS #### 
  # --> NEED TO MAKE IT IN GEE - COMPUTE TMF EXTENT FOR ALL YEARS SO DON'T HAVE TO MAKE REMAINING HERE 
  tmfext <- brick(here("input_data", "TMF", paste0("TMF_extent_classes1-2_3km_",CNT,".tif")))
  
  # Select only years used in this project
  tmfext <- tmfext[[wanted_TMF_layers]]
  
  tmfext_aggr_output_name <- here("temp_data", "processed_tmf", "tmf_aoi", paste0("TMF_extent_classes1-2_9km_",CNT,"_",t0,"_",tT,".tif"))
  
  # aggregate it from the 3km cells to ~9km
  raster::aggregate(tmfext, fact = 3, # c(res(croped_gaez)[1]/res(tmfext)[1], res(croped_gaez)[2]/res(tmfext)[2]),
                    expand = TRUE,
                    fun = sum,
                    na.rm = TRUE, # NA values are only in the sea. Where there is no forest loss, like in a city in Brazil, the value is 0 (see with plot())
                    filename = tmfext_aggr_output_name,
                    # datatype = "INT2U", # let the data be float, as we have decimals in the amount of hectares. 
                    overwrite = TRUE)

  #### EXTRACT VALUES IN COUNTRIES #### 

  # Read layers to be stacked
  agri <- brick(aggregated_ouput_nameS["agri"])
  plantation <- brick(aggregated_ouput_nameS["plantation"])
  flood <- brick(aggregated_ouput_nameS["flood"])
  tmfext <- brick(tmfext_aggr_output_name)
  
  # It is important to explicitly rename layers that are going to be stacked and then called to reshape the data frame 
  # for time varying variables, the dot is important. 
  names(agri) <- paste0("tmf_agri.",seq(t0, tT, 1)) 
  names(plantation) <- paste0("tmf_plantation.",seq(t0, tT, 1)) 
  names(flood) <- paste0("tmf_flood.",seq(t0, tT, 1)) 
  names(tmfext) <- paste0("tmf_extent.",seq(t0, tT, 1)) 
  
  # Stack together the annual layers TMF data and GAEZ crop and pasture2000 cross sections 
  continental_stack <- stack(agri,
                             plantation, 
                             flood, 
                             tmfext)
  # save those names 
  continental_stack_names <- names(continental_stack)
  
  # Extract sum in countries 
  country_tmf <- raster::extract(x = continental_stack, 
                                  y = countries, 
                                  fun = sum, 
                                  na.rm = TRUE, 
                                  layer = 1, 
                                  nl = nlayers(continental_stack))
  removeTmpFiles(h=0)
  
  country_tmf %>% class()
  ncol(country_tmf) == nlayers(continental_stack)
  nrow(country_tmf) == nrow(countries)
  
  colnames(country_tmf) <- names(continental_stack)
  # summary(country_tmf[,"tmf_ext.2018"])
  
  country_tmf_df <- as.data.frame(country_tmf) 
  country_tmf_df$country_name <- countries$COUNTRY_NA
  
  country_tmf_df$continent_name <- CNT

  # the dot is, by construction of all variable names, only in the names of time varying variables. 
  # fixed = TRUE is necessary (otherwise the dot is read as a regexp I guess)
  # Note also that it is important that it is structured in a LIST when there are several varying variables in the *long* format
  # Because: "Notice that the order of variables in varying is like x.1,y.1,x.2,y.2."
  varying_vars <- list(names(agri),
                       names(plantation),
                       names(flood),
                       names(tmfext))
  
  # reshape to long.
  long_country_tmf_df <- stats::reshape(country_tmf_df,
                                        varying = varying_vars,
                                        v.names = c("tmf_agri", "tmf_plantation", "tmf_flood", "tmf_extent"),
                                        sep = ".",
                                        timevar = "year",
                                        idvar = "country_name", # don't put "lon" and "lat" in there, otherwise memory issue (see https://r.789695.n4.nabble.com/reshape-makes-R-run-out-of-memory-PR-14121-td955889.html)
                                        ids = "country_name", # lonlat is our cross-sectional identifier.
                                        direction = "long",
                                        new.row.names = NULL)#seq(from = 1, to = nrow(ibs_msk_df)*length(years), by = 1)
  
  names(long_country_tmf_df)

  # Some polishing
  years <- seq(t0, tT, 1)
  long_country_tmf_df <- mutate(long_country_tmf_df, year = years[year])
  
  long_country_tmf_df <- dplyr::arrange(long_country_tmf_df, country_name, year)

  long_df_list[[CNT]] <- long_country_tmf_df
  
  rm(agri, plantation, flood, tmfext)
  rm(aggregated_ouput_nameS, tmfext_aggr_output_name, wide_df, long_country_tmf_df, country_tmf_df, country_tmf) 
} # closes the loop over continents

# stack data frames from each continent
tropical_long_df <- bind_rows(long_df_list)

# total deforestation (for any driver)
tropical_long_df <- dplyr::mutate(tropical_long_df, tmf_deforestation = tmf_agri + tmf_plantation + tmf_flood)

# country-year obs. are duplicated across continents, due to the loop over continents.
# (the spatial join was over the full set of countries each time) 
# but among the 3 duplicates, it's not always the case that 2 are empty and one is full, because a country can extract values from two continent-tmf-rasters
# Therefore, I collapse the values across continents, for every country. 
tropical_long_df <- readRDS(here("temp_data", paste0("tmf_countries_pantrop_nomask_",t0,"_",tT,".Rdata"))) 

country_tmf_final <- tropical_long_df %>% ddply(.variables = c("country_name", "year"), 
                                                .fun = summarise, 
                                                tmf_extent = sum(tmf_extent, na.rm = TRUE),
                                                tmf_deforestation = sum(tmf_deforestation, na.rm = TRUE),
                                                tmf_agri = sum(tmf_agri, na.rm = TRUE),
                                                tmf_plantation = sum(tmf_plantation, na.rm = TRUE),
                                                tmf_flood = sum(tmf_flood, na.rm = TRUE))
# check it's a balanced panel 
length(unique(country_tmf_final$country_name))*length(unique(country_tmf_final$year)) == nrow(country_tmf_final)

# remove also countries that are completely outside of the tropical belt, and thus completely empty of information 
anyNA(country_tmf_final)

empty <- country_tmf_final %>% 
                            ddply(.variables = c("country_name"), 
                                  .fun = summarise, 
                                  empty = if_else(sum(tmf_extent)==0, 1, 0))

country_tmf_final <- left_join(country_tmf_final, empty, "country_name")
country_tmf_final <- dplyr::filter(country_tmf_final, empty == 0)
country_tmf_final <- dplyr::select(country_tmf_final, -empty)

# check it's still a balanced panel 
length(unique(country_tmf_final$country_name))*length(unique(country_tmf_final$year)) == nrow(country_tmf_final)

saveRDS(country_tmf_final, here("temp_data", paste0("tmf_countries_pantrop_nomask_",t0,"_",tT,".Rdata")))

write_xlsx(country_tmf_final, here("temp_data", paste0("tmf_countries_pantrop_nomask_",t0,"_",tT,".xlsx")))


country_tmf_final <- readRDS(here("temp_data", paste0("tmf_countries_pantrop_nomask_",t0,"_",tT,".Rdata")))



### Some checks against statistics in Vancutsem Science article, at continental level #### 

# Deforestation, including after degradation, but excluding followed by regrowth

# America
agri <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("tmf_agri_9km_America_1990_2020.tif")))
plantation <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("tmf_plantation_9km_America_1990_2020.tif")))
flood <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("tmf_flood_9km_America_1990_2020.tif")))
any <- overlay(stack(agri, plantation, flood), 
               fun = sum)
vany <- values(any)
sum(vany)/1e6 # = 88.5
# pretty well matches the 88.7 Mha deforested in Latin America over 1990-2020 in Vancutsem et al. 2021 (Science), Tab. 10 first column 

# Asia
agri <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("tmf_agri_9km_Asia_1990_2020.tif")))
plantation <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("tmf_plantation_9km_Asia_1990_2020.tif")))
flood <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("tmf_flood_9km_Asia_1990_2020.tif")))
any <- overlay(stack(agri, plantation, flood), 
               fun = sum)
vany <- values(any)
sum(vany)/1e6 # = 59.8
# pretty well matches the 60.9 Mha deforested in Asia-Oceania over 1990-2020 in Vancutsem et al. 2021 (Science), Tab. 10 first column 

# TMF, undisturbed and degraded, America in 1990
tmfext <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("TMF_extent_classes1-2_9km_America_1990_2020.tif")))
tmfext90 <- tmfext[[1]]
vtmfext90 <- values(tmfext90)
sum(vtmfext90)/1e6
# 719 Mha, so a bit more than the 705.1 Mha in Vancutsem et al. Tab. 2. 
# Maybe they removed some regions from what they call "Latin America", like French Guyane for instance

# TMF, undisturbed and degraded, Africa in 1990
tmfext <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("TMF_extent_classes1-2_9km_Africa_1990_2020.tif")))
tmfext90 <- tmfext[[1]]
vtmfext90 <- values(tmfext90)
sum(vtmfext90)/1e6
# 280 Mha, so a bit more than the 273 Mha in Vancutsem et al. Tab. 2. 

# TMF, undisturbed and degraded, Asia in 1990
tmfext <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("TMF_extent_classes1-2_9km_Asia_1990_2020.tif")))
tmfext90 <- tmfext[[1]]
vtmfext90 <- values(tmfext90)
sum(vtmfext90)/1e6
# 322 Mha, so a bit more than the 311 Mha in Vancutsem et al. Tab. 2. 


