

##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
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
# troublePackages <- c() 
# # Attempt to load packages from user's default libraries.
# lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

# 3. If the troubling packages could not be loaded ("there is no package called ...") 
#   you should try to install them, preferably in their versions stated in the renv.lock file. 
#   see in particular https://rstudio.github.io/renv/articles/renv.html 


# # # /!\ THIS BREAKS THE PROJECT REPRODUCIBILITY GUARANTY /!\
# troublePackages <- c("leaflet", "leaflet.providers", "png")
# # Attempt to load packages from user's default libraries.
# lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing

### NEW FOLDERS USED IN THIS SCRIPT 
dir.create(here("temp_data","processed_trade_exposures"))


# NOTE THERE ARE BIG DIFFERENCES NOW BETWEEN THIS SCRIPT AND prepare_exposures.R in RFSFOOD project. 
# all the dictionaries below give, in particular, the final names that should match eaear_* crop names in GLOBILUC project. 

# adding up nutritional contents of these categories yields the grand total. 
categories <- c("Cereals - Excluding Beer" = "all_cereals", 
                "Starchy Roots" = "roots", 
                "Sugar Crops" = "sugarcane", # NOTE sugarcane here, to match the choice made to study only sugar cane 
                "Sugar & Sweeteners" = "sweeteners",
                "Pulses" = "pulses",
                "Treenuts" = "treenuts",
                "Oilcrops" = "oilcrops",
                "Vegetable Oils" = "vegetable_oils",
                "Vegetables" = "vegetables",
                "Fruits - Excluding Wine" = "fruits",
                "Stimulants" = "stimulants",
                "Spices" = "spices",
                "Alcoholic Beverages" = "alcool",
                "Meat" = "fodder", # NOTE fodder here, to match other parts of the code
                "Offals" = "offals",
                "Animal fats" = "animal_fats",
                "Eggs" = "eggs",
                "Milk - Excluding Butter" = "milk",
                "Fish, Seafood" = "fish",
                "Aquatic Products, Other" = "aquatic_products"
                # "Miscellaneous" exclude Miscellaneous category because blind nutrient conversion not relevant for it, and the item is not available for all elements (missing on Production)
)


cereals <- c("Wheat and products" = "wheat", 
             "Rice (Milled Equivalent)" = "rice", 
             "Barley and products" = "barley",
             "Maize and products" = "maizegrain",
             "Rye and products" = "rye", 
             "Oats" = "oats", 
             "Millet and products" = "millet", 
             "Sorghum and products" = "sorghum", 
             "Cereals, Other" = "cereals_other")

oilcrops <- c("Soyabeans" = "soybeans",
              "Groundnuts (Shelled Eq)" = "groundnuts",
              "Rape and Mustardseed" = "rapeseed", 
              "Sunflower seed" = "sunflowerseed", 
              "Cottonseed" = "cottonseed",
              "Coconuts - Incl Copra" = "coconuts",
              "Sesame seed" = "sesameseed",
              "Palm kernels" = "palm_kernels",
              "Olives (including preserved)" = "olives",
              "Oilcrops, Other" = "oilcrops_other")

vegetable_oils <- c("Soyabean Oil" = "soybean_oil",
                    "Cottonseed Oil" = "cottonseed_oil",
                    "Groundnut Oil" = "groundnut_oil",
                    "Sunflowerseed Oil" = "sunflowerseed_oil", 
                    "Rape and Mustard Oil" = "rapeseed_oil",
                    
                    "Palm Oil" = "oilpalm", 
                    "Coconut Oil" = "coconut",
                    "Sesameseed Oil" = "sesameseed_oil",
                    "Olive Oil" = "olive_oil",
                    "Maize Germ Oil" = "maizegerm_oil",
                    "Oilcrops Oil, Other" = "oilcrops_oil_other")

others <- c("Lemons, Limes and products" = "citrus", 
            "Coffee and products" = "coffee",
            "Cocoa Beans and products" = "cocoa", 
            "Tea (including mate)" = "tea",
            "Bananas" = "banana")

all_items_fb <- c("Grand Total" = "grand_total", 
               categories,  
               cereals, 
               oilcrops,
               vegetable_oils, 
               others)

# there are some repetitions (*seed) with food balances that we'll need to remove from commodity balances later
all_items_cb <- c("Brans" = "brans",
                  "Cotton lint" = "cotton_lint", 
                  "Cottonseed" = "cottonseed", 
                  "Cottonseed Cake" = "cottonseed_cake", 
                  "Jute" = "jute", 
                  "Oilseed Cakes, Other" = "oilseed_cake_other", 
                  "Sesameseed Cake" = "sesameseed_cake", 
                  "Silk" = "silk",
                  "Sunflowerseed Cake" = "sunflowerseed_cake", 
                  "Tobacco" = "tobacco", 
                  "Wool (Clean Eq.)" = "whool", 
                  "Alcohol, Non-Food" = "alcohol",     
                  "Copra Cake" = "copra_cake", 
                  "Groundnut Cake" = "groundnut_cake", 
                  "Hard Fibres, Other" = "hard_fibres_other", 
                  "Palm kernels" = "palm_kernels", 
                  "Palmkernel Cake" = "palm_kernel_cake", 
                  "Rape and Mustard Cake" = "rapeseed_cake", 
                  "Rape and Mustardseed" = "rapeseed",
                  "Rubber" = "rubber", 
                  "Sisal" = "sisal", 
                  "Soft-Fibres, Other" = "soft_fibres_other", 
                  "Soyabean Cake" = "soybean_cake", 
                  "Jute-Like Fibres" = "jute_fibres",  
                  "Abaca" = "abaca")

all_items <- c(all_items_fb, all_items_cb)

# The selected items are all those we need to match the crops. Some will be aggregated first 
selected_items <- c("Meat", 
                    
                    "Maize and products",
                    
                    # for those, we will just average with equal weights
                    "Wheat and products", 
                    "Barley and products",
                    "Sorghum and products",
                    
                    "Rice (Milled Equivalent)", 
                    
                    # for soy, cotton, groundnuts, rapeseed and sunflowerseed, we will convert to seed equivalents, then sum up and take one aggregated  measure. 
                    "Soyabeans",
                    "Soyabean Oil",
                    "Soyabean Cake",
                    
                    "Cottonseed Oil",
                    "Cotton lint", 
                    "Cottonseed",
                    "Cottonseed Cake",
                    
                    "Groundnuts (Shelled Eq)",
                    "Groundnut Oil",
                    "Groundnut Cake",
                    "Rape and Mustardseed", 
                    "Rape and Mustard Oil",
                    "Rape and Mustard Cake",
                    "Sunflower seed", 
                    "Sunflowerseed Oil", 
                    "Sunflowerseed Cake",
                    
                    "Sugar Crops",
                    
                    "Tobacco", 
                    
                    "Lemons, Limes and products", 
                    
                    "Palm Oil", 
                    "Palm kernels", 
                    "Palmkernel Cake",
                    
                    "Coconut Oil", 
                    
                    "Coffee and products",
                    "Cocoa Beans and products",
                    
                    "Tea (including mate)",
                    
                    "Bananas", 
                    
                    "Rubber")

final_items <- c("fodder", 
                 "maizegrain",
                 "cereals",
                 "rice", 
                 "soy_compo",
                 "cotton",
                 "oilfeed_crops",
                 "sugarcane",
                 "tobacco", 
                 "citrus", 
                 "oilpalm", 
                 "coconut", 
                 "cocoa_coffee",
                 "tea",
                 "banana", 
                 "rubber")

# note that we won't need food elements here, but it's just that data is from the same download as the food security project
all_elements_fb <- c("Production (1000 tonnes)" = "production_ktonnes",
                     "Import Quantity (1000 tonnes)" = "import_ktonnes",
                     "Stock Variation (1000 tonnes)" = "stock_var_ktonnes", 
                     "Export Quantity (1000 tonnes)" = "export_ktonnes", 
                     "Domestic supply quantity (1000 tonnes)" = "domestic_supply_ktonnes",
                     "Food (1000 tonnes)" = "fsupply_ktonnes", 
                     "Food supply quantity (kg/capita/yr)" = "fsupply_kg/capita/day",
                     "Food supply (kcal/capita/day)" = "fsupply_kcal/capita/day", 
                     "Protein supply quantity (g/capita/day)" = "fsupply_gprot/capita/day",
                     "Fat supply quantity (g/capita/day)" = "fsupply_gfat/capita/day")

# they are expressed in tonnes, initially, but are converted here in ktonnes, thus we match the initial names (tonnes) to *_ktonnes
all_elements_cb <- c("Production (tonnes)" = "production_ktonnes",
                  "Import Quantity (tonnes)" = "import_ktonnes",
                  "Stock Variation (tonnes)" = "stock_var_ktonnes", 
                  "Export Quantity (tonnes)" = "export_ktonnes", 
                  "Domestic supply quantity (tonnes)" = "domestic_supply_ktonnes",
                  "Food supply quantity (tonnes)" = "fsupply_ktonnes")

#### PREPARE FOOD BALANCES  --------------------------------------------------------------------------------------------

pretreatment_years <- c(2001:2007)
year <- 2001
# list to store annual data sets prepared
wide_fb_list <- list()

for(year in pretreatment_years){
  fb <- read.csv(here("input_data", "trade_exposure", paste0("FAOSTAT-oldfoodbalances_allcountries_aggritems_",year,".csv")))
  
  # unique(fb$Item)[grepl("meal", unique(fb$Item))]
  
  # the first row names gets weird "ï.." prefixe
  # names(fb)[1] <- "Domain.Code"
  
  # the domain is the common to the whole data set "Food Balances (-2013, old methodology and population)" so we can remove it
  unique(fb$Domain) 
  
  fb <- fb[, !grepl("Domain", names(fb))]
  
  ## Area and area code are bijective 
  length(unique(fb$Area)) == length(unique(fb$Area.Code)) 
  fb <- dplyr::select(fb, -Area.Code)
  
  names(fb)[names(fb)=="Area"] <- "country_name"
  
  ## Element
  unique(fb$Element)
  # "Food supply (kcal/capita/day)"          
  # "Protein supply quantity (g/capita/day)" 
  # "Fat supply quantity (g/capita/day)"   
  # "Production"
  # "Import Quantity"      
  # "Stock Variation" 
  # "Export Quantity" 
  # "Domestic supply quantity"               
  # "Food"                                  
  # "Food supply quantity (kg/capita/yr)"   
  
  # Element and Element.Code are bijective 
  length(unique(fb$Element)) == length(unique(fb$Element.Code))
  fb <- dplyr::select(fb, -Element.Code)
  
  ## Item 
  unique(fb$Item)
  
  # there are some Items that have more than one Item.Code
  length(unique(fb$Item)) == length(unique(fb$Item.Code))
  
  spot_doublones <- sapply(unique(fb$Item), function(itm){itm_length <- fb[fb$Item==itm,"Item.Code"] %>% unique() %>% length() 
  itm_spot <- if_else(itm_length > 1, true = itm, false = "")
  return(itm_spot) })
  spot_doublones[spot_doublones != ""]
  
  fb[fb$Item == "Eggs","Item.Code"] %>% unique()
  fb[fb$Item == "Milk - Excluding Butter","Item.Code"] %>% unique()
  fb[fb$Item == "Miscellaneous","Item.Code"] %>% unique()
  
  fb[fb$Item == "Eggs",] 
  fb[fb$Item == "Miscellaneous",] 
  
  # they seem to be the same figures, or almost, but with different Flags (not always) --> let's keep only one instance
  
  # Items with several item codes are duplicates (within the same country_name and Element)
  fb[duplicated(fb[,c("country_name", "Element", "Item")]), ] %>% nrow() # 2951 obs. in 2001
  # we want those that are not duplicates
  fb <- fb[!duplicated(fb[,c("country_name", "Element", "Item")]), ]
  
  if(!(length(unique(fb$Item)) == length(unique(fb$Item.Code)))){
    stop("there are still duplicates in Item variable")
  }
  
  fb <- dplyr::select(fb, -Item.Code)
  
  # keep only Items specified
  fb <- dplyr::filter(fb, Item %in% names(all_items_fb))
  
  ## Year
  fb <- dplyr::select(fb, -Year.Code)
  # Year is not useful either
  fb <- dplyr::select(fb, -Year)
  
  
  ## Unit 
  # For "Import Quantity", "Domestic supply quantity" & "Food", the unit is not given in the Element. 
  
  # some checks that the units are expressed in a sound way
  prod_u <- fb[fb$Element=="Production", "Unit"] %>% unique()
  import_u <- fb[fb$Element=="Import Quantity", "Unit"] %>% unique()
  export_u <- fb[fb$Element=="Export Quantity", "Unit"] %>% unique()
  stock_u <- fb[fb$Element=="Stock Variation", "Unit"] %>% unique()
  dom_supply_u <- fb[fb$Element=="Domestic supply quantity", "Unit"] %>% unique()
  food_u <- fb[fb$Element=="Food", "Unit"] %>% unique()
  
  
  if(length(import_u) > 1 | length(dom_supply_u) > 1 | length(food_u) > 1 | 
     length(prod_u) > 1 | length(export_u) > 1 | length(stock_u) > 1 ){
    stop("different units used within Elements")
  }
  if(!all.equal(import_u, dom_supply_u, food_u, prod_u, export_u, stock_u)){
    stop("different units used across Elements")
  }
  if(import_u != "1000 tonnes"){
    stop("different units used across YEARS")
  }
  
  u_less_slct <- fb$Element %in% c("Production", "Import Quantity", "Export Quantity", "Stock Variation", 
                                   "Domestic supply quantity", "Food")
  fb[u_less_slct,] <- mutate(fb[u_less_slct,], Element = paste0(Element, " (",Unit,")"))
  
  fb <- dplyr::select(fb, -Unit)
  
  ## For the moment, do not bother Flags
  fb <- fb[, !grepl("Flag", names(fb))]
  
  ## Change item and element strings
  fb$Item <- sapply(fb$Item, FUN = function(i){spaceless <- all_items_fb[i]
  names(spaceless) <- NULL
  return(spaceless)})
  
  fb$Element <- sapply(fb$Element, FUN = function(i){spaceless <- all_elements_fb[i]
  names(spaceless) <- NULL
  return(spaceless)})
  
  ## RESHAPE
  # First split data by Element
  unique(fb$Element)
  # "Food (1000 tonnes)" is the annual quantity available
  # "Food supply quantity (kg/capita/yr)" is the annual quantity available, but divided by population
  elmt_wide_ds_list <- list()
  for(elmt in unique(fb$Element)){
    long_ds <- fb[fb$Element==elmt, c("country_name", "Item", "Value")]
    
    wide_ds <- stats::reshape(long_ds,
                              # varying = unique(long_ds$Item),
                              # v.names = c("Value"),
                              sep = ".",
                              timevar = "Item",
                              idvar = "country_name", 
                              direction = "wide",
                              new.row.names = NULL)  
    
    vars_slct <- grepl("Value.", names(wide_ds))
    
    # those variables that have been reshaped, give the Element identifier to their names
    names(wide_ds)[vars_slct] <- paste0(elmt,"_",names(wide_ds)[vars_slct])
    
    # remove "Value." part in names
    names(wide_ds)[vars_slct] <- gsub("Value.", "", 
                                      x = names(wide_ds)[vars_slct])
    
    elmt_wide_ds_list[[elmt]] <- wide_ds
  }
  rm(wide_ds)
  
  # and then join them back based on country_name key
  wide_fb <- elmt_wide_ds_list[[1]]
  for(i in 2:length(elmt_wide_ds_list)){
    wide_fb <- left_join(wide_fb, elmt_wide_ds_list[[i]], by = "country_name")
  }
  # at this point, 175 rows, one for each country_name, and, if no Item has been removed, 808 columns, one for each type Element*Item 
  
  length(all_items_fb)*7 == ncol(wide_fb) - 1
  # not all selected items are available for every 7 elements. 
  # in particular, grand total and vegetable and animal product totals are available only in nutrient, not in weights 
  
  ### Clean some country_name related things
  unique(wide_fb$country_name)
  # Handle China: get Taiwan apart (makes sense in food security context)
  # wide_fb[grepl("China", wide_fb$country_name), c("country_name", "rice_production_ktonnes")]
  wide_fb$country_name[wide_fb$country_name == "China, Taiwan Province of"] <- "Taiwan"
  
  # and remove Hong Kong and Macao 
  wide_fb <- dplyr::filter(wide_fb, country_name != "China, Hong Kong SAR")
  wide_fb <- dplyr::filter(wide_fb, country_name != "China, Macao SAR") 
  
  # Remove China (which aggregates China mainland and Taiwan), and keep only China mainland 
  wide_fb <- dplyr::filter(wide_fb, country_name != "China")
  
  # and call it just China to match main_data
  wide_fb$country_name[wide_fb$country_name == "China, mainland"] <- "China"
  
  # Remove oversea territories
  # wide_fb[grepl("Fr", wide_fb$country_name), c("country_name", "rice_production_ktonnes")]
  wide_fb <- dplyr::filter(wide_fb, country_name != "French Polynesia")
  wide_fb <- dplyr::filter(wide_fb, country_name != "Netherlands Antilles (former)")
  wide_fb <- dplyr::filter(wide_fb, country_name != "Bermuda")
  wide_fb <- dplyr::filter(wide_fb, country_name != "New Caledonia")
  
  # Handle some weird names 
  # "TÃ¼rkiye" and "CÃ´te d'Ivoire" 
  wide_fb$country_name[wide_fb$country_name == "Türkiye"] <- "Turkey"
  wide_fb$country_name[wide_fb$country_name == "TÃ¼rkiye"] <- "Turkey"
  wide_fb$country_name[wide_fb$country_name == "T?rkiye"] <- "Turkey"
  # wide_fb[grepl("Tur", wide_fb$country_name), c("country_name", "rice_production_ktonnes")]
  # wide_fb[grepl("?", wide_fb$country_name), c("country_name", "rice_production_ktonnes")]
  wide_fb$country_name[wide_fb$country_name == "Côte d'Ivoire"] <- "Côte d'Ivoire"
  wide_fb$country_name[wide_fb$country_name == "CÃ´te d'Ivoire"] <- "Côte d'Ivoire"
  wide_fb$country_name[wide_fb$country_name == "C?te d'Ivoire"] <- "Côte d'Ivoire"
  
  # CONGO - this is special: there is only one Congo in the data, named simply "Congo". 
  # Checking by the population size in the data, it is the Republic of the Congo (i.e. Congo Brazzaville)
  # Which is named Congo, as well, in main_data

  # take only Serbia, after it splitted with Montenegro, and call it as it was prior splitting
  # grep(pattern = "erbia", x = unique(wide_fb$country_name), value = TRUE) 
  wide_fb$country_name[wide_fb$country_name=="Serbia"] <- "Serbia and Montenegro" # Montenegro independent since 2006, (almost only) after our data period
  wide_fb <- dplyr::filter(wide_fb, country_name != "Montenegro")

  # just change sudan so that it matches main_data, don't bother more
  wide_fb$country_name[wide_fb$country_name == "Sudan (former)"] <- "South Sudan"
  
  wide_fb$country_name[wide_fb$country_name == "Gambia"] <- "Gambia, The"
  
  # Keep track of the year 
  wide_fb$year <- year
  
  wide_fb_list[[match(year, pretreatment_years)]] <- wide_fb
} # Close loop over years here, because we will need imputations from other years for NAs in the next steps

# stack to make a panel of food balance data  
pfb <- bind_rows(wide_fb_list)
# rm(wide_fb_list, wide_fb)
# unique(pfb$country_name)




#  Data above miss RDC, LIBYA, SYRIA & PAPUA NEW GUINEA (entre autres mineurs) 
# mais tant pis, on ne les prend pas dans 2010- FAOSTAT, car ces pays manquent aussi dans commodity balances, 
# qui est plus important car ça comprend tous les meals. 
# (il vaut mieux avoir une varaiable mieux mesurée - i.e. qui comprend les usages non-food comme meals, mais qui ne couvre pas tous les pays, que l'inverse, surtout étant donnés les pays qui manquent, qui ne sont pas cruciaux en termes de déforestation)

#### PREPARE COMMODITY BALANCES --------------------------------------------------------------------------------------------------

wide_cb_list <- list()

for(year in pretreatment_years){
  cb <- read.csv(here("input_data", "trade_exposure", paste0("FAOSTAT-commoditybalances_allcountries_allitems_",year,".csv")))
  
  # unique(cb$Item)[grepl("meal", unique(cb$Item))]
  
  # the first row names gets weird "ï.." prefixe
  # names(cb)[1] <- "Domain.Code"
  
  # the domain is the common to the whole data set "Food Balances (-2013, old methodology and population)" so we can remove it
  unique(cb$Domain) 
  
  cb <- cb[, !grepl("Domain", names(cb))]
  
  ## Area and area code are bijective 
  length(unique(cb$Area)) == length(unique(cb$Area.Code..FAO.)) 
  cb <- dplyr::select(cb, -Area.Code..FAO.)
  
  names(cb)[names(cb)=="Area"] <- "country_name"
  
  ## Element
  unique(cb$Element)
  # "Food supply (kcal/capita/day)"          
  # "Protein supply quantity (g/capita/day)" 
  # "Fat supply quantity (g/capita/day)"   
  # "Production"
  # "Import Quantity"      
  # "Stock Variation" 
  # "Export Quantity" 
  # "Domestic supply quantity"               
  # "Food"                                  
  # "Food supply quantity (kg/capita/yr)"   
  
  # Element and Element.Code are bijective 
  length(unique(cb$Element)) == length(unique(cb$Element.Code))
  cb <- dplyr::select(cb, -Element.Code)
  
  ## Item 
  unique(cb$Item)
  
  # there are some Items that have more than one Item.Code
  length(unique(cb$Item)) == length(unique(cb$Item.Code))
  
  spot_doublones <- sapply(unique(cb$Item), function(itm){itm_length <- cb[cb$Item==itm,"Item.Code"] %>% unique() %>% length() 
  itm_spot <- if_else(itm_length > 1, true = itm, false = "")
  return(itm_spot) })
  spot_doublones[spot_doublones != ""]
  
  
  # they seem to be the same figures, or almost, but with different Flags (not always) --> let's keep only one instance
  
  # Items with several item codes are duplicates (within the same country_name and Element)
  cb[duplicated(cb[,c("country_name", "Element", "Item")]), ] %>% nrow() # 2951 obs. in 2001
  # we want those that are not duplicates
  cb <- cb[!duplicated(cb[,c("country_name", "Element", "Item")]), ]
  
  if(!(length(unique(cb$Item)) == length(unique(cb$Item.Code)))){
    stop("there are still duplicates in Item variable")
  }
  
  cb <- dplyr::select(cb, -Item.Code)
  
  # keep only Items specified
  cb <- dplyr::filter(cb, Item %in% names(all_items_cb))
  
  ## Year
  cb <- dplyr::select(cb, -Year.Code)
  # Year is not useful either
  cb <- dplyr::select(cb, -Year)
  
  
  ## Unit 
  # For "Import Quantity", "Domestic supply quantity" & "Food", the unit is not given in the Element. 
  
  # some checks that the units are expressed in a sound way
  # NOTE: EVERYTHING IS IN "tonnes" NOT IN 1000 tonnes  AS IN FOOD BALANCES
  prod_u <- cb[cb$Element=="Production", "Unit"] %>% unique()
  import_u <- cb[cb$Element=="Import Quantity", "Unit"] %>% unique()
  export_u <- cb[cb$Element=="Export Quantity", "Unit"] %>% unique()
  stock_u <- cb[cb$Element=="Stock Variation", "Unit"] %>% unique()
  dom_supply_u <- cb[cb$Element=="Domestic supply quantity", "Unit"] %>% unique()
  food_u <- cb[cb$Element=="Food supply quantity (tonnes)", "Unit"] %>% unique()
  
  
  if(length(import_u) > 1 | length(dom_supply_u) > 1 | length(food_u) > 1 | 
     length(prod_u) > 1 | length(export_u) > 1 | length(stock_u) > 1 ){
    stop("different units used within Elements")
  }
  if(!all.equal(import_u, dom_supply_u, food_u, prod_u, export_u, stock_u)){
    stop("different units used across Elements")
  }
  if(import_u != "tonnes"){
    stop("different units used across YEARS")
  }
  
  u_less_slct <- cb$Element %in% c("Production", "Import Quantity", "Export Quantity", "Stock Variation", 
                                   "Domestic supply quantity")# don't feature "Food supply quantity (tonnes)" here, as the unit is already in the original name
  cb[u_less_slct,] <- mutate(cb[u_less_slct,], Element = paste0(Element, " (",Unit,")"))
  
  cb <- dplyr::select(cb, -Unit)
  
  # so elements are ALL in tonnes here in commodity balances, but we want them all in ktonnes (and names accordingly)
  cb <- mutate(cb, Value = Value/1000)
  # naming is just below, with all_elements_cb sapply
  
  ## For the moment, do not bother Flags
  cb <- cb[, !grepl("Flag", names(cb))]
  
  ## Change item and element strings
  cb$Item <- sapply(cb$Item, FUN = function(i){spaceless <- all_items_cb[i]
  names(spaceless) <- NULL
  return(spaceless)})
  
  cb$Element <- sapply(cb$Element, FUN = function(i){spaceless <- all_elements_cb[i]
  names(spaceless) <- NULL
  return(spaceless)})
  
  ## RESHAPE
  # First split data by Element
  unique(cb$Element)
  # "Food (1000 tonnes)" is the annual quantity available
  # "Food supply quantity (kg/capita/yr)" is the annual quantity available, but divided by population
  elmt_wide_ds_list <- list()
  for(elmt in unique(cb$Element)){
    long_ds <- cb[cb$Element==elmt, c("country_name", "Item", "Value")]
    
    wide_ds <- stats::reshape(long_ds,
                              # varying = unique(long_ds$Item),
                              # v.names = c("Value"),
                              sep = ".",
                              timevar = "Item",
                              idvar = "country_name", 
                              direction = "wide",
                              new.row.names = NULL)  
    
    vars_slct <- grepl("Value.", names(wide_ds))
    
    # those variables that have been reshaped, give the Element identifier to their names
    names(wide_ds)[vars_slct] <- paste0(elmt,"_",names(wide_ds)[vars_slct])
    
    # remove "Value." part in names
    names(wide_ds)[vars_slct] <- gsub("Value.", "", 
                                      x = names(wide_ds)[vars_slct])
    
    elmt_wide_ds_list[[elmt]] <- wide_ds
  }
  rm(wide_ds)
  
  # and then join them back based on country_name key
  wide_cb <- elmt_wide_ds_list[[1]]
  for(i in 2:length(elmt_wide_ds_list)){
    wide_cb <- left_join(wide_cb, elmt_wide_ds_list[[i]], by = "country_name")
  }
  # at this point, 175 rows, one for each country_name, and, if no Item has been removed, 808 columns, one for each type Element*Item 
  
  length(all_items_cb)*7 == ncol(wide_cb) - 1
  # not all selected items are available for every 7 elements. 
  # in particular, grand total and vegetable and animal product totals are available only in nutrient, not in weights 
  
  ### Clean some country_name related things
  unique(wide_cb$country_name)
  # Handle China: get Taiwan apart (makes sense in food security context)
  # wide_cb[grepl("China", wide_cb$country_name), c("country_name", "rice_production_ktonnes")]
  
  wide_cb$country_name[wide_cb$country_name == "China, Taiwan Province of"] <- "Taiwan"
  
  # and remove Hong Kong and Macao 
  wide_cb <- dplyr::filter(wide_cb, country_name != "China, Hong Kong SAR")
  wide_cb <- dplyr::filter(wide_cb, country_name != "China, Macao SAR") 
  
  # Remove China (which aggregates China mainland and Taiwan), and keep only China mainland 
  wide_cb <- dplyr::filter(wide_cb, country_name != "China")

  # and call it just China to match main_data
  wide_cb$country_name[wide_cb$country_name == "China, mainland"] <- "China"
  
  # Remove oversea territories
  # wide_cb[grepl("Fr", wide_cb$country_name), c("country_name", "rice_production_ktonnes")]
  wide_cb <- dplyr::filter(wide_cb, country_name != "French Polynesia")
  wide_cb <- dplyr::filter(wide_cb, country_name != "Netherlands Antilles (former)")
  wide_cb <- dplyr::filter(wide_cb, country_name != "Bermuda")
  wide_cb <- dplyr::filter(wide_cb, country_name != "New Caledonia")
  
  # Handle some weird names 
  # "TÃ¼rkiye" and "CÃ´te d'Ivoire" 
  wide_cb$country_name[wide_cb$country_name == "Türkiye"] <- "Turkey"
  wide_cb$country_name[wide_cb$country_name == "TÃ¼rkiye"] <- "Turkey"
  wide_cb$country_name[wide_cb$country_name == "T?rkiye"] <- "Turkey"
  # wide_cb[grepl("Tur", wide_cb$country_name), c("country_name", "rice_production_ktonnes")]
  # wide_cb[grepl("?", wide_cb$country_name), c("country_name", "rice_production_ktonnes")]
  wide_cb$country_name[wide_cb$country_name == "Côte d'Ivoire"] <- "Côte d'Ivoire"
  wide_cb$country_name[wide_cb$country_name == "CÃ´te d'Ivoire"] <- "Côte d'Ivoire"
  wide_cb$country_name[wide_cb$country_name == "C?te d'Ivoire"] <- "Côte d'Ivoire"
  
  # take only Serbia, after it splitted with Montenegro, and call it as it was prior splitting
  # grep(pattern = "erbia", x = unique(wide_cb$country_name), value = TRUE) 
  wide_cb$country_name[wide_cb$country_name=="Serbia"] <- "Serbia and Montenegro" # Montenegro independent since 2006, (almost only) after our data period
  wide_cb <- dplyr::filter(wide_cb, country_name != "Montenegro")
  
  # just change sudan so that it matches main_data, don't bother more
  wide_cb$country_name[wide_cb$country_name == "Sudan (former)"] <- "South Sudan"
  
  wide_cb$country_name[wide_cb$country_name == "Gambia"] <- "Gambia, The"
  
  # Keep track of the year 
  wide_cb$year <- year
  
  wide_cb_list[[match(year, pretreatment_years)]] <- wide_cb
} # Close loop over years here, because we will need imputations from other years for NAs in the next steps

# stack to make a panel of food balance data  
pcb <- bind_rows(wide_cb_list)

# remove commodities that are already in food balances under same name (palm_kernels, cottonseed and rapeseed) (and preserve country_name and year variables)
pcb <- pcb[,!(names(pcb) %in% names(pfb)[!(names(pfb)%in%c("country_name", "year"))])]

# merge commodity and food balances 
pfb <- full_join(pfb, pcb, by = c("country_name", "year"))

# remove food elements
pfb <- pfb[, !grepl("fsupply_", x = names(pfb))]

# for convenience
row.names(pfb) <- dplyr::mutate(pfb, country_year = paste0(country_name, "_", year))$country_year


pfb_save <- pfb

#### AGGREGATE COMMODITIES -----------------------------------------------------------------------------------------
all_items[selected_items]

# see sources and discussion of oil/cake content in prepare_loss_commodity.R 
# Remove NAs, implies to assume NAs are ignorable by their small sizes 

for(elmt in c("production_ktonnes_", "export_ktonnes_", "import_ktonnes_", "stock_var_ktonnes_")){
  pfb[, paste0(elmt,"seed_eq_soybean_oil")] <-   pfb[, paste0(elmt,"soybean_oil")]*(1/0.18)
  pfb[, paste0(elmt,"seed_eq_soybean_cake")] <-   pfb[, paste0(elmt,"soybean_cake")]*(1/0.82)
  pfb <- mutate(pfb, !!as.symbol(paste0(elmt,"soy_compo")) := rowSums(across(.cols = any_of(paste0(elmt,c("soybeans", 
                                                                                                           "seed_eq_soybean_oil",
                                                                                                           "seed_eq_soybean_cake")))), 
                                                                      na.rm = TRUE))

  pfb[, paste0(elmt,"seed_eq_groundnut_oil")] <-   pfb[, paste0(elmt,"groundnut_oil")]*(1/0.31)
  pfb[, paste0(elmt,"seed_eq_groundnut_cake")] <-   pfb[, paste0(elmt,"groundnut_cake")]*(1/(1-0.31))  
  pfb <- mutate(pfb, !!as.symbol(paste0(elmt,"groundnut_compo")) := rowSums(across(.cols = any_of(paste0(elmt,c("groundnuts", 
                                                                                                          "seed_eq_groundnut_oil",
                                                                                                          "seed_eq_groundnut_cake")))), 
                                                                      na.rm = TRUE))

  
  pfb[, paste0(elmt,"seed_eq_rapeseed_oil")] <-   pfb[, paste0(elmt,"rapeseed_oil")]*(1/0.41)
  pfb[, paste0(elmt,"seed_eq_rapeseed_cake")] <-   pfb[, paste0(elmt,"rapeseed_cake")]*(1/(1-0.41))
  pfb <- mutate(pfb, !!as.symbol(paste0(elmt,"rapeseed_compo")) := rowSums(across(.cols = any_of(paste0(elmt,c("rapeseed", 
                                                                                                          "seed_eq_rapeseed_oil",
                                                                                                          "seed_eq_rapeseed_cake")))), 
                                                                      na.rm = TRUE))
  
  
  pfb[, paste0(elmt,"seed_eq_sunflowerseed_oil")] <-   pfb[, paste0(elmt,"sunflowerseed_oil")]*(1/0.42)
  pfb[, paste0(elmt,"seed_eq_sunflowerseed_cake")] <-   pfb[, paste0(elmt,"sunflowerseed_cake")]*(1/(1-0.42))
  pfb <- mutate(pfb, !!as.symbol(paste0(elmt,"sunflowerseed_compo")) := rowSums(across(.cols = any_of(paste0(elmt,c("sunflowerseed", 
                                                                                                               "seed_eq_sunflowerseed_oil",
                                                                                                               "seed_eq_sunflowerseed_cake")))), 
                                                                           na.rm = TRUE))
  
  # and add oil crops up 
  pfb <- mutate(pfb, !!as.symbol(paste0(elmt,"oilfeed_crops")) := rowSums(across(.cols = any_of(paste0(elmt,c("groundnut_compo", 
                                                                                                              "rapeseed_compo",
                                                                                                              "sunflowerseed_compo")))), 
                                                                                na.rm = TRUE))
  
  # source for oil and cake is https://cottonaustralia.com.au/uses-of-cotton#:~:text=single%20bed%20sheets-,Cotton%20seed,meal%20and%20300kg%20of%20hulls.
  pfb[, paste0(elmt,"seed_eq_cottonseed_oil")] <-   pfb[, paste0(elmt,"cottonseed_oil")]*(1/0.2)
  pfb[, paste0(elmt,"seed_eq_cottonseed_cake")] <-   pfb[, paste0(elmt,"cottonseed_cake")]*(1/(0.5))
  pfb[, paste0(elmt,"seed_eq_cotton_lint")] <-   pfb[, paste0(elmt,"cotton_lint")]*(1/(0.33)) # this is the conversion factor from GAEZ
  pfb <- mutate(pfb, !!as.symbol(paste0(elmt,"cotton")) := rowSums(across(.cols = any_of(paste0(elmt,c("cottonseed", 
                                                                                                      "seed_eq_cottonseed_oil",
                                                                                                      "seed_eq_cottonseed_cake", 
                                                                                                      "seed_eq_cotton_lint")))), 
                                                                                na.rm = TRUE))
  
  # And for cereals and cocoa and coffee, weadd up directly
  pfb <- mutate(pfb, !!as.symbol(paste0(elmt,"cereals")) := rowSums(across(.cols = any_of(paste0(elmt,c("barley", 
                                                                                                        "sorghum",
                                                                                                        "wheat")))), 
                                                                    na.rm = TRUE))
  
  
  pfb <- mutate(pfb, !!as.symbol(paste0(elmt,"cocoa_coffee")) := rowSums(across(.cols = any_of(paste0(elmt,c("cocoa", "coffee")))), 
                                                                    na.rm = TRUE))

}


#### MAKE STATISTICS OF INTEREST ---------------------------------------------------------------------------------

# Exposure of selected crops
# (making sure that outputs are NA, not Inf or NaN)
for(item in final_items){

  pfb <- dplyr::mutate(pfb, 
                       !!as.symbol(paste0("trade_expo_",item)) := ( !!as.symbol(paste0("export_ktonnes_", item)) + !!as.symbol(paste0("import_ktonnes_", item)) ) /
                                                                    !!as.symbol(paste0("production_ktonnes_", item)),  
                                                                      
                       !!as.symbol(paste0("export_expo_",item)) := ( !!as.symbol(paste0("export_ktonnes_", item)) ) / !!as.symbol(paste0("production_ktonnes_", item)) ) 
  
  # handle cases where production is 0
  pfb[, paste0("trade_expo_",item)][!is.finite(pfb[, paste0("trade_expo_",item)])] <- NA
  pfb[, paste0("export_expo_",item)][!is.finite(pfb[, paste0("export_expo_",item)])] <- NA
   # this can produce annual NA but not a pb, as na.rm = TRUE handles it in ddply below. 

}

### AVERAGE OVER YEARS 
# We average over different sets of years, as it is not clear in advance what is the most appropriate period (there are trade-offs)
pretreat_year_sets <- list(`2001_2007` = c(2001:2007), 
                           `2006_2007` = c(2006:2007))


csfb <- list()
for(pretreat_period in pretreat_year_sets){
  
  csfb <- ddply(pfb[pfb$year %in% c(pretreat_period), ], 
                "country_name", summarise, 
                
                # Exposure of specific crops, either from all trade or from exports only
                !!as.symbol("trade_expo_fodder") := mean(!!as.symbol("trade_expo_fodder"), na.rm = TRUE), 
                !!as.symbol("trade_expo_maizegrain") := mean(!!as.symbol("trade_expo_maizegrain"), na.rm = TRUE), 
                !!as.symbol("trade_expo_cereals") := mean(!!as.symbol("trade_expo_cereals"), na.rm = TRUE), 
                !!as.symbol("trade_expo_rice") := mean(!!as.symbol("trade_expo_rice"), na.rm = TRUE), 
                !!as.symbol("trade_expo_soy_compo") := mean(!!as.symbol("trade_expo_soy_compo"), na.rm = TRUE), 
                !!as.symbol("trade_expo_cotton") := mean(!!as.symbol("trade_expo_cotton"), na.rm = TRUE), 
                !!as.symbol("trade_expo_oilfeed_crops") := mean(!!as.symbol("trade_expo_oilfeed_crops"), na.rm = TRUE), 
                !!as.symbol("trade_expo_sugarcane") := mean(!!as.symbol("trade_expo_sugarcane"), na.rm = TRUE), 
                !!as.symbol("trade_expo_tobacco") := mean(!!as.symbol("trade_expo_tobacco"), na.rm = TRUE), 
                !!as.symbol("trade_expo_citrus") := mean(!!as.symbol("trade_expo_citrus"), na.rm = TRUE), 
                !!as.symbol("trade_expo_oilpalm") := mean(!!as.symbol("trade_expo_oilpalm"), na.rm = TRUE),
                !!as.symbol("trade_expo_coconut") := mean(!!as.symbol("trade_expo_coconut"), na.rm = TRUE),
                !!as.symbol("trade_expo_cocoa_coffee") := mean(!!as.symbol("trade_expo_cocoa_coffee"), na.rm = TRUE),
                !!as.symbol("trade_expo_tea") := mean(!!as.symbol("trade_expo_tea"), na.rm = TRUE), 
                !!as.symbol("trade_expo_banana") := mean(!!as.symbol("trade_expo_banana"), na.rm = TRUE), 
                !!as.symbol("trade_expo_rubber") := mean(!!as.symbol("trade_expo_rubber"), na.rm = TRUE), 
                
                !!as.symbol("export_expo_fodder") := mean(!!as.symbol("export_expo_fodder"), na.rm = TRUE), 
                !!as.symbol("export_expo_maizegrain") := mean(!!as.symbol("export_expo_maizegrain"), na.rm = TRUE), 
                !!as.symbol("export_expo_cereals") := mean(!!as.symbol("export_expo_cereals"), na.rm = TRUE), 
                !!as.symbol("export_expo_rice") := mean(!!as.symbol("export_expo_rice"), na.rm = TRUE), 
                !!as.symbol("export_expo_soy_compo") := mean(!!as.symbol("export_expo_soy_compo"), na.rm = TRUE), 
                !!as.symbol("export_expo_cotton") := mean(!!as.symbol("export_expo_cotton"), na.rm = TRUE), 
                !!as.symbol("export_expo_oilfeed_crops") := mean(!!as.symbol("export_expo_oilfeed_crops"), na.rm = TRUE), 
                !!as.symbol("export_expo_sugarcane") := mean(!!as.symbol("export_expo_sugarcane"), na.rm = TRUE), 
                !!as.symbol("export_expo_tobacco") := mean(!!as.symbol("export_expo_tobacco"), na.rm = TRUE), 
                !!as.symbol("export_expo_citrus") := mean(!!as.symbol("export_expo_citrus"), na.rm = TRUE), 
                !!as.symbol("export_expo_oilpalm") := mean(!!as.symbol("export_expo_oilpalm"), na.rm = TRUE),
                !!as.symbol("export_expo_coconut") := mean(!!as.symbol("export_expo_coconut"), na.rm = TRUE),
                !!as.symbol("export_expo_cocoa_coffee") := mean(!!as.symbol("export_expo_cocoa_coffee"), na.rm = TRUE),
                !!as.symbol("export_expo_tea") := mean(!!as.symbol("export_expo_tea"), na.rm = TRUE), 
                !!as.symbol("export_expo_banana") := mean(!!as.symbol("export_expo_banana"), na.rm = TRUE), 
                !!as.symbol("export_expo_rubber") := mean(!!as.symbol("export_expo_rubber"), na.rm = TRUE)
  )
  
  # make another metric, that imputes the international average to country-crops where all information is missing (i.e. production, export, and import are NA every year)
  for(item in final_items){
    
    csfb[,paste0("trade_expo_imp_",item)] <- csfb[,paste0("trade_expo_",item)]  
    csfb[is.nan(csfb[,paste0("trade_expo_imp_",item)]),paste0("trade_expo_imp_",item)] <- mean(csfb[,paste0("trade_expo_",item)], na.rm = TRUE)
    
    csfb[,paste0("export_expo_imp_",item)] <- csfb[,paste0("export_expo_",item)]  
    csfb[is.nan(csfb[,paste0("export_expo_imp_",item)]),paste0("export_expo_imp_",item)] <- mean(csfb[,paste0("export_expo_",item)], na.rm = TRUE)
    
  }
  
  # Finally, simply add columns for biomass, just to match main_data without error
  csfb$trade_expo_biomass <- NA
  csfb$trade_expo_imp_biomass <- NA
  csfb$export_expo_biomass <- NA
  csfb$export_expo_imp_biomass <- NA
  
  saveRDS(csfb, file = here("temp_data", "processed_trade_exposures", paste0("trade_exposures_",
                                                                       min(pretreat_period),
                                                                       max(pretreat_period),".Rdata")))
  
}



rm(fb, pfb, pfb_save, cb, pcb, wide_fb, wide_fb_list,wide_cb, wide_cb_list, csfb, long_ds, elmt_wide_ds_list)



# # used to make some country names match with main_data 
# un_csfb <- unique(csfb$country_name)
# un_md <- unique(main_data$country_name)
# 
# un_csfb[!(un_csfb %in% un_md)] #  many, normal puisqu'on n'a que les pays tropicaux dans un_md
# 
# un_md[!(un_md %in% un_csfb)]
# 
# grep(pattern="urund", x=un_csfb, value = TRUE)








