# Script to prepare CROPGRIDS data 

### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
neededPackages <- c("data.table", "plyr", "tidyr", "dplyr",  "Hmisc", "sjmisc", "stringr",
                    "here", "readstata13", "foreign", "readxl", "writexl",
                    "raster", "terra", "sp", "spdep", "sf","gfcanalysis",  "nngeo", "stars", # "osrm", "osrmr",
                    "lubridate","exactextractr",
                    "doParallel", "foreach", "snow", "parallel",
                    "knitr", "kableExtra",
                    "DataCombine", 
                    "fixest", 
                    "boot", "fwildclusterboot", "sandwich", "MASS",
                    "ggplot2", "ggpubr", "leaflet",  "dotwhisker", "viridis", "hrbrthemes")
# "tmap",
# "pglm", "multiwayvcov", "clusterSEs", "alpaca", "clubSandwich",

# Install them in their project-specific versions
renv::restore(packages = neededPackages)

# Often useful to upgrade renv to debug: 
# renv::upgrade()

# Load them
lapply(neededPackages, library, character.only = TRUE)

dir.create(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "grouped"), recursive = TRUE)
dir.create(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "masks"), recursive = TRUE)
dir.create(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "shares"), recursive = TRUE)

### Raster global options
rasterOptions(timer = TRUE, 
              progress = "text")

### GLOBAL CRS USED ### 
mercator_world_crs <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "

### TROPICAL AOI 
tropical_aoi <- extent(c(-180, 179.9167, -30, 30))

land <- st_read(here("input_data", "ne_50m_land"))
unique(land$scalerank)
land <- land[land$scalerank==0, c("geometry")]

# The target raster to mask 
# croplandcommo <- brick(
#   list.files(path = here("input_data", "10thLoss_croplandcommo_3km"), 
#                          # pattern = paste0("^10thLoss_",type,"_3km"), 
#                          full.names = TRUE) %>% as.list()
# )
lossdriver = brick(here("temp_data", "processed_lossdrivers", "tropical_aoi", paste0("loss_croplandcommo_3km_0119.tif")))

# CROPGRIDS names mapped to group names.  
gridcrops = c(
  # "alfalfa", "mixedgrass", # we will do pastures separately, because not handled in cropgrids apparently.
  # "oilpalm" # no need either
  # "banana", "plantain", # for now, don't bother cat. 6 crops. 

  # Rice
  Rice = "rice",  
  # Cereals there's also "rye", "ryefor",   "millet", cerealnes
  Cereals = "barley", 
  Cereals = "sorghumfor", 
  Cereals = "wheat", 
  Cereals = "buckwheat", 

  # Roots crops - there's also "rootnes"
  Roots = "cassava", 
  Roots = "potato", 
  Roots = "sweetpotato", 
  Roots = "yam", 
  # Soybean
  Soy_compo = "soybean",
  # Oil crops there's also "oilseed", "oilseedfor", not sure what this captures? 
  Oilfeed_crops = "rapeseed", 
  Oilfeed_crops = "sunflower", 
  Oilfeed_crops = "groundnut", 
  # Cotton
  Cotton = "cotton",   
  # Biomass crops
  Biomass = "sorghum", 
  Biomass = "canaryseed", 
  # Sugar
  Sugarcane = "sugarcane", # "sugarbeet", 
  # Citrus - NOT IN croplands from Potatpov, so don't include here
  # Citrus = "citrusnes",
  
  Maizegrain = "maize", 
  Maizegrain = "maizefor",
  
  # Tobacco
  Tobacco = "tobacco"
)

(all_group_names = unique(names(gridcrops)))


# Turn to TIF ----------
for(cropgridcrop in gridcrops){
  cg = rast(here("input_data", "CROPGRIDSv1.08_NC_maps", "CROPGRIDSv1.08_NC_maps", paste0("CROPGRIDSv1.08_",cropgridcrop,".nc")))
  
  # keep only croparea (in HECTARES)
  cg = cg$croparea

  # crop to trop ical area
  cg = crop(cg, tropical_aoi)
  
  # write as tif, to read in raster, while reclassifying negative values
  classify(cg, 
             rcl = cbind(-Inf,0, 0), # Negative values are in the SEA OR LAND, where value is 0... 
           # I haven't found where they explain these values, but I guess it is to differ between where the model ran and estimated 0 and where it did'nt run;  
             filename = here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", paste0(cropgridcrop, ".tif")), 
             overwrite =TRUE) 
}


# Match crop Types ----------
for(cropgroup in all_group_names){
  group_crops = gridcrops[names(gridcrops)==cropgroup]
  # necessary to split the process, because calc(sum) behaves differently if x is a single layer or a stack/brick. 
  if(length(group_crops) == 1){
    group_layer = raster(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", paste0(group_crops, ".tif")))
    names(group_layer) = group_crops
    
    # write
    writeRaster(group_layer, 
                filename = here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "grouped", paste0(cropgroup, ".tif")), 
                overwrite =TRUE)
  }
  else{
    # Read all crops within group
    group_names = as.list(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", paste0(group_crops, ".tif")))
    group_stack = raster::stack(group_names)
    names(group_stack) = group_crops
    
    # sum up their crop area in each cell
    calc(group_stack, 
         fun = sum, na.rm = TRUE, # not necessary anymore, but was for the NAs introduced in reclassification above
         filename = here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "grouped", paste0(cropgroup, ".tif")), 
         overwrite =TRUE)
  }
  rm(group_crops, cropgroup)
}
cereals = raster(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "grouped", paste0("Cereals", ".tif")))
cereals %>% values() %>% summary()
rice = raster(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "grouped", paste0("Rice", ".tif")))
rice %>% values() %>% summary()


## Repeat for all crops --------------
all_group_paths = as.list(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "grouped", paste0(all_group_names, ".tif")))
all_group_stack = raster::stack(all_group_paths)
# sum up their crop area
calc(all_group_stack, 
     fun = sum, 
     filename = here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "grouped", "All_cropland_crops.tif"), 
     overwrite =TRUE
)

all_group_stack$Maizegrain %>% area() %>% values() %>% max() # This in km2, i.e. 3077 hectares.  
# --> confirms that areas of every crop is such that they can be summed up. 
all_group_stack$Tobacco %>% area() %>% values() %>% max() # This in km2, i.e. 3077 hectares.  


# PREPARE SHARE & MASK OF MAIN CROPGRIDS CROP -------------

allcrops = raster(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "grouped", "All_cropland_crops.tif"))
allcrops
plot(allcrops)  
plot(land, add = T)

# There are several ways to do this... 

## As a binary mask ----------------
# Here, we consider that to attribute a loss-to-commodity cropland event to a specific crop, 
# this crop must be the main one, i.e. have the largest area compared to the other crops that can count as commodity cropland. 

# But before using which.max, we need to cast differently all those cells in the sea etc. 
# which have value 0 or NA for all crops will be given the index of the "first layer with maximum value for this cell"
# and this behavior is an issue. 
# max_layer_idx = raster::which.max(all_group_stack)
# 
# # PLOT
# max_layer_idx_df <- as.data.frame(max_layer_idx, xy=TRUE)
# # Create a factor for the raster values
# max_layer_idx_df$layer_fact <- factor(max_layer_idx_df$layer, levels = 1:11, labels = all_group_names)
# max_layer_idx_df$layer_fact %>% table()
# 
# # Plot using ggplot
# ggplot(max_layer_idx_df, aes(x=x, y=y, fill=layer_fact)) +
#   geom_raster() +
#   # scale_fill_manual(values = rainbow(11), name = "Crops") +
#   scale_colour_brewer(na.value = "transparent",
#                       type = "qual",
#                       direction = -1,
#                       palette = "Set3", 
#                       name = "Crops") +
#   labs(title = "Main crop in area by cell, among those shown, in CROPGRIDS") +
#   theme_minimal()
# 
# # make binary masks
# all_group_names = names(all_group_stack)
# for(cropgroup in all_group_names){
#   index_target = match(cropgroup, all_group_names)
#   calc(
#     x = max_layer_idx,
#     fun = function(x){
#       return(if_else(x == index_target, 1, 0))
#     }, 
#     filename = here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "masks", paste0(all_group_names[index_target], "_mask.tif")), 
#     overwrite = TRUE 
#   )
# }

## As a continuous mask ------------
# An alternative is to weight annually lost area to commodity cropland by 
# the share of every crop area in the total cropland area in the cell at the end of the period. 
# For this, we need to make a layer for every crop, where values are this crop's share
tot = raster(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "grouped", "All_cropland_crops.tif"))

# Stack all individual groups and the total 
# tot_groups = stack(
#   tot,
#   all_group_stack
# )

make_shares = function(x, y){
  # Handling 0s in the denominator like this crashes overlay, 
  # and apparently this is handle by raster by setting NA to all those cases
  # if(y == 0){
  #   ret = 0
  # }
  # else{
    ret = x / y
    # so we remove these NAs (I have checked that this yields same sum in the output values as without.  
    ret[!is.finite(ret)] = 0 
  # }
  return(ret)
}

# make binary masks
for(cropgroup in all_group_names){
  group_layer = raster(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "grouped", paste0(cropgroup, ".tif")))
  overlay(
    x = group_layer,
    y = tot,
    fun = make_shares,
    filename = here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "shares", paste0(cropgroup, "_share.tif")),
    overwrite = TRUE
  )
}
# This introduces quite many NAs
# there is none in x or y, but ~7M in the resulting shares!
group_layer %>% values() %>% summary()
tot %>% values() %>% summary()
rice = raster(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "shares", paste0("Rice", "_share.tif")))
rice %>% values() %>% summary()

roots = raster(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "shares", paste0("Roots", "_share.tif")))
roots %>% values() %>% sum(na.rm = TRUE)

plot(rice)

# REMOVE THESE NAs, to align with the other outcome maps



## Downscale to ~3km and align to target loss raster as prepared in GEE ---------

# It makes no difference to disaggregate only now, because the within-cell distribution of crop area is not explicit, 
# i.e., if there's 50/50 between two crops in the 5km cell, it's going to be exactly the same distribution in the cell divided by 2. 
# or, if a crop is the max in the 5km res. raster, it is going to be exactly in the same subcells in the disaggregated one. 

# So, we don't have to disaggregate and resample all groups 
# ... but only the mask 
for(cropgroup in all_group_names){
  group_share = raster(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "shares", paste0(cropgroup, "_share.tif")))
  group_share_disag = disaggregate(group_share, fact = 2)
  resample(x = group_share_disag,
           y = lossdriver,
           filename = here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "shares", paste0(cropgroup, "_share_resampled10thLoss.tif")),
           overwrite =TRUE)
  
  # group_mask = raster(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "masks", paste0(cropgroup, "_mask.tif")))
  # group_mask_disag = disaggregate(group_mask, fact = 2)
  # resample(x = group_mask_disag,
  #          y = lossdriver,
  #          filename = here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "masks", paste0(cropgroup, "_mask_resampled10thLoss.tif")),
  #          overwrite =TRUE)
  
}
# store values for inspection 
# # store values for inspection
# stats_list = list(resampled_stats = list(), 
#                   native_stats = list())
# for(cropgroup in all_group_names){
#   native = raster(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "masks", paste0(cropgroup, "_mask.tif")))
#   resampled = raster( here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "masks", paste0(cropgroup, "_mask_resampled10thLoss.tif")))
#   stats_list[["native_stats"]][[cropgroup]] = native %>% values() %>% summary()
#   stats_list[["resampled_stats"]][[cropgroup]] = resampled %>% values() %>% summary()
# }
# Note: resampling introduces a couple (33k out of 30M) of Nas.  
# but means are very similar. 
# stats_list

# maizeMask_native = raster(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "masks", paste0("Maizegrain_mask.tif")))
# maizeMask_resampled = raster(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "masks", paste0("Maizegrain_mask_3km_resampled10thLoss.tif")))
# maizeMask_native %>% values() %>% summary()
# maizeMask_resampled %>% values() %>% summary()
# 
# native = raster(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "masks", paste0("Citrus_mask.tif")))
# resampled = raster(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "masks", paste0("Citrus_mask_resampled10thLoss.tif")))
# native %>% values() %>% summary()
# resampled %>% values() %>% summary()


# APPLY MASK ---------

# Mask loss-to-commodity cropland by the main CROPGRIDS crop in cell

for(cropgroup in all_group_names){
  cg_share = raster(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "shares", paste0(cropgroup, "_share_resampled10thLoss.tif")))
  
  # write in there because this is where prepare_loss_cropcommoMain_opindusnotrans.R reads to aggregate and resample to GAEZ.
  overlay(x =lossdriver, 
          y = cg_share, # will be recycled to every layer of lossdriver, i.e. every year.
          fun = function(x,y){x*y}, 
          recycle = TRUE,
          filename = here("temp_data", "processed_lossdrivers", "tropical_aoi", paste0("loss_", cropgroup, "_CGweighted", "_3km_0119.tif")),
          overwrite = TRUE
  )
  
  # Binary mask
  # cg_mask = raster(here("temp_data", "CROPGRIDSV1.08", "tropical_aoi", "masks", paste0(cropgroup, "_mask_resampled10thLoss.tif")))
  # # not good practice to write in input_data, but this is to fit the naming convention inherited from the 
  # # workflow for other outcomes measurement that were produced in GEE and written in input_data
  # dir.create(here("input_data", "10thLoss_croplandcommoCGmain_3km"))
  # mask(lossdriver, 
  #      mask = cg_mask, 
  #      maskvalue = 0, # the value in the mask, that identifies cells in the input that will have a different value in the output
  #      updatevalue = 0, # the different value
  #      filename = here("input_data", "10thLoss_croplandcommoCGmain_3km", paste0("10thLoss_", cropgroup, "_CGmain", "_3km.tif")),
  #      overwrite = TRUE
  #      )
}
# Check
# ld_rice_weighted = brick(here("input_data", "10thLoss_croplandcommoCGweighted_3km", paste0("10thLoss_", cropgroup, "_CGweighted", "_3km.tif")))
lossdriver

# lossdriver$loss_croplandcommo_3km_0119_1 %>% values %>% mean(na.rm = TRUE)
# ld_rice_weighted$X10thLoss_Tobacco_CGweighted_3km_1 %>% values() %>% mean(na.rm = TRUE)
# cg_share %>% values() %>% mean(na.rm = TRUE)



