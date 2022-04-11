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
dir.create(here("temp_data", "processed_tmf", "tmf_aoi"), recursive = TRUE)
dir.create(here("temp_data", "processed_pasture2000", "tmf_aoi"), recursive = TRUE)
dir.create(here("temp_data", "merged_datasets", "tmf_aoi"), recursive = TRUE)

### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")

### GLOBAL CRS USED ### 
mercator_world_crs <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "


### GAEZ OBJECTS
# in this script, GAEZ is the target raster of all aggregations / resamplings
gaez_dir <- here("temp_data", "GAEZ", "v4", "AEAY_out_density",  "Rain-fed")
gaez_crops <- list.files(path = here(gaez_dir, "High-input"), 
                         pattern = "", 
                         full.names = FALSE)
gaez_crops <- gsub(pattern = ".tif", replacement = "", x = gaez_crops)



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
transition_types <- c("agri", "plantation", "flood")
continents <- c("America", "Africa", "Asia")#  

for(CNT in continents){ # the order of the loops matter for the mask making operation, see below. 
  
  ### AGGREGATE AND alignE TO GAEZ RESOLUTION ####
  # Currently, the data are hectares of deforestation in 3x3km grid cells annually. 
  # We aggregate this to ~9km grid cells by adding up the hectares of deforestation
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
    
    rm(aggr_output_name)
    
  } # closes the loop over deforestation types
  
  # output nameS of the above loop (this is a length-3 vector)
  aggregated_ouput_nameS <- here("temp_data", "processed_tmf", "tmf_aoi",  paste0("tmf_",transition_types,"_9km_",CNT,"_",t0,"_",tT,".tif"))
  names(aggregated_ouput_nameS) <- transition_types
  
  
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
  names(tmfext) <- paste0("tmf_ext.",seq(t0, tT, 1)) 
  
  # Stack together the annual layers TMF data and GAEZ crop and pasture2000 cross sections 
  continental_stack <- stack(agri,
                             plantation, 
                             flood, 
                             tmfext)
  # save those names 
  continental_stack_names <- names(continental_stack)
  
  # Extract sum in countries 
  countries <- st_read(here("input_data", "Global_LSIB_Polygons_Detailed"))
  length(unique(countries$COUNTRY_NA)) == nrow(countries)
  
  country_tmf <- raster::extract(x = continental_stack, 
                                  y = countries, 
                                  fun = sum, 
                                  na.rm = TRUE, 
                                  layer = 1, 
                                  nl = nlayers(continental_stack))
                  
  
  ### RASTER TO DATAFRAME ### 
  
  continental_stack <- stack(masked_stack_output_name)
  
  # all names are lost through the writing/reading operation, so rename layers 
  names(continental_stack) <- continental_stack_names
  
  # na.rm = TRUE does not matter much here, because there is no NA. 
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
  
  rm(agri, plantation, flood, tmfext)
  
  rm(aggregated_ouput_nameS, tmfext_aggr_output_name, wide_df, long_df) 
} # closes the loop over continents


tropical_long_df <- bind_rows(long_df_list)

saveRDS(tropical_long_df, here("temp_data", "merged_datasets", "tmf_aoi", paste0("tmf_pantrop_nomask_long_",t0,"_",tT,".Rdata")))

rm(tropical_long_df, long_df_list)




#### COUNTRY ANNUAL AVERAGES #### 
final <- readRDS(here("temp_data", "merged_datasets", "tmf_aoi", paste0("tmf_aeay_pantrop_long_final_",t0,"_",tT,".Rdata")))

country_final <- dplyr::select(final, year, continent_name, country_name, country_year, tmf_agri, tmf_plantation, tmf_flood, tmf_ext)
rm(final)
country_avges <- ddply(country_final, "country_year", summarise, 
                       agriurba_defo = sum(tmf_agri, na.rm = TRUE),
                       plantation_defo = sum(tmf_plantation, na.rm = TRUE),
                       flood_defo = sum(tmf_flood, na.rm = TRUE),
                       tmf_extent = sum(tmf_ext, na.rm = TRUE), 
                       .progress = "text")

country_year_final <- country_final[!duplicated(country_final$country_year),] 

country_avges <- inner_join(country_avges, country_year_final[,c("year", "continent_name", "country_name", "country_year")], by = "country_year")

country_avges <- dplyr::mutate(country_avges, deforestation = agriurba_defo + plantation_defo + flood_defo)
country_avges <- dplyr::mutate(country_avges, deforestation = deforestation/1e6)
country_avges <- dplyr::mutate(country_avges, tmf_extent = tmf_extent/1e6)



# some checks against statistics in Vancutsem Science article
# In Tab. 10, they report 39.6, 88.7, and 60.9 Mha deforestation over 1990-2020, in Africa, America and Asia resp.
# (excluding regrowth, but including prior degradation)
country_avges %>% filter(continent_name=="Africa") %>% dplyr::select(deforestation) %>% sum()
country_avges %>% filter(continent_name=="America") %>% dplyr::select(deforestation) %>% sum()
country_avges %>% filter(continent_name=="Asia") %>% dplyr::select(deforestation) %>% sum()
# here, we obtain 33.8, 83.4, and 56 Mha respectively

# Using aggregated and resampled rasters instead 
# Resampled raster: 83.4 Mha too 
agri <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("resampledgaez_tmf_agri_America_1990_2020.tif")))
plantation <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("resampledgaez_tmf_plantation_America_1990_2020.tif")))
flood <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("resampledgaez_tmf_flood_America_1990_2020.tif")))
any <- overlay(stack(agri, plantation, flood), 
               fun = sum)
vany <- values(any)
sum(vany)/1e6

# Aggregated raster: 88.47 Mha, so very very similar to Vancutsem et al. 
agri <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("aggrgaez_tmf_agri_America_1990_2020.tif")))
plantation <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("aggrgaez_tmf_plantation_America_1990_2020.tif")))
flood <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("aggrgaez_tmf_flood_America_1990_2020.tif")))
any <- overlay(stack(agri, plantation, flood), 
               fun = sum)
vany <- values(any)
sum(vany)/1e6

# So the difference comes from the resampling. 

# resampling with bilinear: 83.3
croped_gaez <- crop(gaez, ext_list[["America"]])
bilany <- resample(any, croped_gaez, method = "bilinear")
sum(values(bilany))/1e6

# and ngb (what we used, but layer by layer) : 83.4
ngbany <- resample(any, croped_gaez, method = "ngb")
sum(values(ngbany))/1e6

# so it's really just a matter of resampling, not the method used, or what is done before; 
# try project raster instead of resample: same results. 
ngbany <- projectRaster(any, croped_gaez, method = "bilinear")
sum(values(ngbany))/1e6

# trying different resampling from terra: does not help. 
any <- terra::rast(any)
croped_gaez <- terra::rast(croped_gaez)
terra::resample(any, croped_gaez, method = "near") %>% values() %>% sum()/1e6 # same as raster::resample
terra::resample(any, croped_gaez, method = "bilinear") %>% values() %>% sum()/1e6 # same as raster::resample
terra::resample(any, croped_gaez, method = "cubic") %>% values() %>% sum()/1e6 # 83.279
terra::resample(any, croped_gaez, method = "cubicspline") %>% values() %>% sum()/1e6 # 83.278
terra::resample(any, croped_gaez, method = "lanczos") %>% values() %>% sum()/1e6 # 83.281




# TMF EXTENT 
# In Tab. 2, they report 273.4, 705.1, 311,1 Mha TMF (undisturbed and degraded) in 1990 in Africa, America and ASia resp. 
# They do not count regrowth, while I incude it. However, in 1990 this is not observed in their data anyways. 
# Since we masked places where no deforestation ever occured, we should not use the final data frame to make such comparisons. 
# Resampled raster: 677.07 Mha
am_tmfext <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("resampledgaez_tmfext_America_1990_2020.tif")))
am_tmf90 <- values(am_tmfext[[1]])
sum(am_tmf90)/1e6
# Aggregated raster: 719.6 Mha
am_tmfext <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("aggrgaez_tmfext_America_1990_2020.tif")))
am_tmf90 <- values(am_tmfext[[1]])
sum(am_tmf90)/1e6

# with bilinear resampling, very similar sum: 677.29
croped_gaez <- crop(gaez, ext_list[["America"]])
aggr <- am_tmfext[[1]]
resamp_bil <- resample(aggr, croped_gaez, method = "bilinear") 
sum(values(resamp_bil))/1e6

# Resampled Africa: 263.3 Mha
af_tmfext <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("resampledgaez_tmfext_Africa_1990_2020.tif")))
af_tmf90 <- values(af_tmfext[[1]])
sum(af_tmf90)/1e6

# Resampled Asia: 303.4 Mha
as_tmfext <- brick( here("temp_data", "processed_tmf", "tmf_aoi", paste0("resampledgaez_tmfext_Asia_1990_2020.tif")))
as_tmf90 <- values(as_tmfext[[1]])
sum(as_tmf90)/1e6



# Old tests --------------

# they report 238.7, 616.4 and 244.1 Mha of TMF (undisturbed and degraded) in 2015 in Africa, Latin America, and Asia-Oceania resp. 
country_avges %>% filter(continent_name=="Africa" & year == 1989) %>% dplyr::select(tmf_extent) %>% sum()
country_avges %>% filter(continent_name=="America" & year == 1989) %>% dplyr::select(tmf_extent) %>% sum()
country_avges %>% filter(continent_name=="Asia" & year == 1989) %>% dplyr::select(tmf_extent) %>% sum()

head(country_avges)
brazil <- dplyr::filter(country_avges, country_name == "Brazil" & year > 2000)
indonesia <- dplyr::filter(country_avges, country_name == "Indonesia" & year > 2000)


plot(x = brazil$year, y = brazil$deforestation)
plot(x = indonesia$year, y = indonesia$deforestation)
plot(x = brazil$year, y = brazil$tmf_extent)
plot(x = indonesia$year, y = indonesia$tmf_extent)

# Compare with FIg. 3 in Vancutsem et al. 2021 Science. Indonesian panel
# The present measure of deforestation is to be compare with the sum of 
# - direct defo not followed by regrowth
# - conversion to plantations
# - conversion to water
# - defor after degradatation not followed by regrowth (i.e. degradation -> no regrowth -> defor)
# defor after degradation followed by regrowth (i.e. i.e. degradation -> regrowth -> defor)
# degradation before defor


America <- dplyr::filter(country_avges, continent_name == "America" & year > 1994)
Asia <- dplyr::filter(country_avges, continent_name == "Asia" & year > 1994)
America <- ddply(America, "year", summarise, 
                 deforestation = sum(deforestation))
Asia <- ddply(Asia, "year", summarise, 
              deforestation = sum(deforestation))
plot(x = America$year, y = America$deforestation)
plot(x = Asia$year, y = Asia$deforestation)


# compare with Table 11
brazil[brazil$year > 2000 & brazil$year<2013,"deforestation"] %>% mean()
America[America$year > 2000 & America$year<2013,"deforestation"] %>% mean()
America[America$year > 2000 & America$year<2020,"deforestation"] %>% mean()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  ### ### ### ### ### ### ### ### ### ### ### 
