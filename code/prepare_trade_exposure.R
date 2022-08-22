

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


# adding up nutritional contents of these categories yields the grand total. 
categories <- c("Cereals - Excluding Beer" = "all_cereals", 
                "Starchy Roots" = "roots", 
                "Sugar Crops" = "sugar_crops",
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
             "Maize and products" = "maize",
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
                    
                    "Palm Oil" = "palm_oil", 
                    "Coconut Oil" = "coconut_oil",
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
                 "maize",
                 "cereals",
                 "rice", 
                 "soy_compo",
                 "cotton",
                 "oilfeed_crops",
                 "sugar_crops",
                 "tobacco", 
                 "citrus", 
                 "palm_oil", 
                 "coconut_oil", 
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
  
  names(fb)[names(fb)=="Area"] <- "country"
  
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
  
  # Items with several item codes are duplicates (within the same country and Element)
  fb[duplicated(fb[,c("country", "Element", "Item")]), ] %>% nrow() # 2951 obs. in 2001
  # we want those that are not duplicates
  fb <- fb[!duplicated(fb[,c("country", "Element", "Item")]), ]
  
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
    long_ds <- fb[fb$Element==elmt, c("country", "Item", "Value")]
    
    wide_ds <- stats::reshape(long_ds,
                              # varying = unique(long_ds$Item),
                              # v.names = c("Value"),
                              sep = ".",
                              timevar = "Item",
                              idvar = "country", 
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
  
  # and then join them back based on country key
  wide_fb <- elmt_wide_ds_list[[1]]
  for(i in 2:length(elmt_wide_ds_list)){
    wide_fb <- left_join(wide_fb, elmt_wide_ds_list[[i]], by = "country")
  }
  # at this point, 175 rows, one for each country, and, if no Item has been removed, 808 columns, one for each type Element*Item 
  
  length(all_items_fb)*7 == ncol(wide_fb) - 1
  # not all selected items are available for every 7 elements. 
  # in particular, grand total and vegetable and animal product totals are available only in nutrient, not in weights 
  
  ### Clean some country related things
  unique(wide_fb$country)
  # Handle China: get Taiwan apart (makes sense in food security context)
  # wide_fb[grepl("China", wide_fb$country), c("country", "rice_production_ktonnes")]
  
  wide_fb$country[wide_fb$country == "China, Taiwan Province of"] <- "Taiwan"
  
  # Remove China (which aggregates China mainland and Taiwan), and keep only China mainland 
  wide_fb <- dplyr::filter(wide_fb, country != "China")
  
  # and remove Hong Kong and Macao 
  wide_fb <- dplyr::filter(wide_fb, country != "China, Hong Kong SAR")
  wide_fb <- dplyr::filter(wide_fb, country != "China, Macao SAR")
  
  
  # Remove oversea territories
  # wide_fb[grepl("Fr", wide_fb$country), c("country", "rice_production_ktonnes")]
  wide_fb <- dplyr::filter(wide_fb, country != "French Polynesia")
  wide_fb <- dplyr::filter(wide_fb, country != "Netherlands Antilles (former)")
  wide_fb <- dplyr::filter(wide_fb, country != "Bermuda")
  wide_fb <- dplyr::filter(wide_fb, country != "New Caledonia")
  
  # Handle some weird names 
  # "TÃ¼rkiye" and "CÃ´te d'Ivoire" 
  wide_fb$country[wide_fb$country == "Türkiye"] <- "Turkey"
  wide_fb$country[wide_fb$country == "TÃ¼rkiye"] <- "Turkey"
  wide_fb$country[wide_fb$country == "T?rkiye"] <- "Turkey"
  # wide_fb[grepl("Tur", wide_fb$country), c("country", "rice_production_ktonnes")]
  # wide_fb[grepl("?", wide_fb$country), c("country", "rice_production_ktonnes")]
  wide_fb$country[wide_fb$country == "Côte d'Ivoire"] <- "Ivory Coast"
  wide_fb$country[wide_fb$country == "CÃ´te d'Ivoire"] <- "Ivory Coast"
  wide_fb$country[wide_fb$country == "C?te d'Ivoire"] <- "Ivory Coast"
  
  # CONGO - this is special: there is only one Congo in the data, named simply "Congo". 
  # Checking by the population size in the data, it is the Republic of the Congo (i.e. Congo Brazzaville)
  wide_fb$country[wide_fb$country == "Congo"] <- "Republic of the Congo"
  
  # take only Serbia, after it splitted with Montenegro, and call it as it was prior splitting
  # grep(pattern = "erbia", x = unique(wide_fb$country), value = TRUE) 
  wide_fb$country[wide_fb$country=="Serbia"] <- "Serbia and Montenegro" # Montenegro independent since 2006, (almost only) after our data period
  wide_fb <- dplyr::filter(wide_fb, country != "Montenegro")
  # The problem does not occur for sudan and south sudan, as the  distinction is not made in the data untile 2009 at least. 
  # thus, only Sudan (former) appears here, and need not be adjusted. 
  # grep(pattern = "udan", x = unique(wide_fb$country), value = TRUE) 
  
  
  # Keep track of the year 
  wide_fb$year <- year
  
  wide_fb_list[[match(year, pretreatment_years)]] <- wide_fb
} # Close loop over years here, because we will need imputations from other years for NAs in the next steps

# stack to make a panel of food balance data  
pfb <- bind_rows(wide_fb_list)
# rm(wide_fb_list, wide_fb)
# unique(pfb$country)




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
  
  names(cb)[names(cb)=="Area"] <- "country"
  
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
  
  # Items with several item codes are duplicates (within the same country and Element)
  cb[duplicated(cb[,c("country", "Element", "Item")]), ] %>% nrow() # 2951 obs. in 2001
  # we want those that are not duplicates
  cb <- cb[!duplicated(cb[,c("country", "Element", "Item")]), ]
  
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
    long_ds <- cb[cb$Element==elmt, c("country", "Item", "Value")]
    
    wide_ds <- stats::reshape(long_ds,
                              # varying = unique(long_ds$Item),
                              # v.names = c("Value"),
                              sep = ".",
                              timevar = "Item",
                              idvar = "country", 
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
  
  # and then join them back based on country key
  wide_cb <- elmt_wide_ds_list[[1]]
  for(i in 2:length(elmt_wide_ds_list)){
    wide_cb <- left_join(wide_cb, elmt_wide_ds_list[[i]], by = "country")
  }
  # at this point, 175 rows, one for each country, and, if no Item has been removed, 808 columns, one for each type Element*Item 
  
  length(all_items_cb)*7 == ncol(wide_cb) - 1
  # not all selected items are available for every 7 elements. 
  # in particular, grand total and vegetable and animal product totals are available only in nutrient, not in weights 
  
  ### Clean some country related things
  unique(wide_cb$country)
  # Handle China: get Taiwan apart (makes sense in food security context)
  # wide_cb[grepl("China", wide_cb$country), c("country", "rice_production_ktonnes")]
  
  wide_cb$country[wide_cb$country == "China, Taiwan Province of"] <- "Taiwan"
  
  # Remove China (which aggregates China mainland and Taiwan), and keep only China mainland 
  wide_cb <- dplyr::filter(wide_cb, country != "China")
  
  # and remove Hong Kong and Macao 
  wide_cb <- dplyr::filter(wide_cb, country != "China, Hong Kong SAR")
  wide_cb <- dplyr::filter(wide_cb, country != "China, Macao SAR")
  
  
  # Remove oversea territories
  # wide_cb[grepl("Fr", wide_cb$country), c("country", "rice_production_ktonnes")]
  wide_cb <- dplyr::filter(wide_cb, country != "French Polynesia")
  wide_cb <- dplyr::filter(wide_cb, country != "Netherlands Antilles (former)")
  wide_cb <- dplyr::filter(wide_cb, country != "Bermuda")
  wide_cb <- dplyr::filter(wide_cb, country != "New Caledonia")
  
  # Handle some weird names 
  # "TÃ¼rkiye" and "CÃ´te d'Ivoire" 
  wide_cb$country[wide_cb$country == "Türkiye"] <- "Turkey"
  wide_cb$country[wide_cb$country == "TÃ¼rkiye"] <- "Turkey"
  wide_cb$country[wide_cb$country == "T?rkiye"] <- "Turkey"
  # wide_cb[grepl("Tur", wide_cb$country), c("country", "rice_production_ktonnes")]
  # wide_cb[grepl("?", wide_cb$country), c("country", "rice_production_ktonnes")]
  wide_cb$country[wide_cb$country == "Côte d'Ivoire"] <- "Ivory Coast"
  wide_cb$country[wide_cb$country == "CÃ´te d'Ivoire"] <- "Ivory Coast"
  wide_cb$country[wide_cb$country == "C?te d'Ivoire"] <- "Ivory Coast"
  
  # CONGO - this is special: there is only one Congo in the data, named simply "Congo". 
  # Checking by the population size in the data, it is the Republic of the Congo (i.e. Congo Brazzaville)
  wide_cb$country[wide_cb$country == "Congo"] <- "Republic of the Congo"
  
  # take only Serbia, after it splitted with Montenegro, and call it as it was prior splitting
  # grep(pattern = "erbia", x = unique(wide_cb$country), value = TRUE) 
  wide_cb$country[wide_cb$country=="Serbia"] <- "Serbia and Montenegro" # Montenegro independent since 2006, (almost only) after our data period
  wide_cb <- dplyr::filter(wide_cb, country != "Montenegro")
  # The problem does not occur for sudan and south sudan, as the  distinction is not made in the data untile 2009 at least. 
  # thus, only Sudan (former) appears here, and need not be adjusted. 
  # grep(pattern = "udan", x = unique(wide_cb$country), value = TRUE) 
  
  
  # Keep track of the year 
  wide_cb$year <- year
  
  wide_cb_list[[match(year, pretreatment_years)]] <- wide_cb
} # Close loop over years here, because we will need imputations from other years for NAs in the next steps

# stack to make a panel of food balance data  
pcb <- bind_rows(wide_cb_list)

# remove commodities that are already in food balances under same name (palm_kernels, cottonseed and rapeseed) (and preserve country and year variables)
pcb <- pcb[,!(names(pcb) %in% names(pfb)[!(names(pfb)%in%c("country", "year"))])]

# merge commodity and food balances 
pfb <- full_join(pfb, pcb, by = c("country", "year"))

# remove food elements
pfb <- pfb[, !grepl("fsupply_", x = names(pfb))]

# for convenience
row.names(pfb) <- dplyr::mutate(pfb, country_year = paste0(country, "_", year))$country_year


pfb_save <- pfb

#### AGGREGATE COMMODITIES -----------------------------------------------------------------------------------------
all_items[selected_items]

# see sources and discussion of oil/cake content in prepare_loss_commodity.R 
# Remove NAs, implies to assume NAs are ignorable by their small sizes 

for(elmt in c("production_ktonnes_", "export_ktonnes_", "import_ktonnes_", "stock_var_ktonnes_")){
  pfb[, paste0(elmt,"seed_eq_soybean_oil")] <-   pfb[, paste0(elmt,"soybean_oil")]*(1/0.18)
  pfb[, paste0(elmt,"seed_eq_soybean_cake")] <-   pfb[, paste0(elmt,"soybean_cake")]*(1/0.82)
  pfb <- mutate(pfb, !!as.symbol(paste0(elmt,"soy_compo")) := sum(c(!!as.symbol(paste0(elmt,"soybeans")), 
                                                                        !!as.symbol(paste0(elmt,"seed_eq_soybean_oil")), 
                                                                        !!as.symbol(paste0(elmt,"seed_eq_soybean_cake"))), na.rm = TRUE) )

  pfb[, paste0(elmt,"seed_eq_groundnut_oil")] <-   pfb[, paste0(elmt,"groundnut_oil")]*(1/0.31)
  pfb[, paste0(elmt,"seed_eq_groundnut_cake")] <-   pfb[, paste0(elmt,"groundnut_cake")]*(1/(1-0.31))  
  pfb <- mutate(pfb, !!as.symbol(paste0(elmt,"groundnut_compo")) := sum(c(!!as.symbol(paste0(elmt,"groundnuts")), 
                                                                        !!as.symbol(paste0(elmt,"seed_eq_groundnut_oil")), 
                                                                        !!as.symbol(paste0(elmt,"seed_eq_groundnut_cake"))), na.rm = TRUE) )
  
  pfb[, paste0(elmt,"seed_eq_rapeseed_oil")] <-   pfb[, paste0(elmt,"rapeseed_oil")]*(1/0.41)
  pfb[, paste0(elmt,"seed_eq_rapeseed_cake")] <-   pfb[, paste0(elmt,"rapeseed_cake")]*(1/(1-0.41))
  pfb <- mutate(pfb, !!as.symbol(paste0(elmt,"rapeseed_compo")) := sum(c(!!as.symbol(paste0(elmt,"rapeseed")), 
                                                                        !!as.symbol(paste0(elmt,"seed_eq_rapeseed_oil")), 
                                                                        !!as.symbol(paste0(elmt,"seed_eq_rapeseed_cake"))), na.rm = TRUE) )
  
  pfb[, paste0(elmt,"seed_eq_sunflowerseed_oil")] <-   pfb[, paste0(elmt,"sunflowerseed_oil")]*(1/0.42)
  pfb[, paste0(elmt,"seed_eq_sunflowerseed_cake")] <-   pfb[, paste0(elmt,"sunflowerseed_cake")]*(1/(1-0.42))
  pfb <- mutate(pfb, !!as.symbol(paste0(elmt,"sunflowerseed_compo")) := sum(c(!!as.symbol(paste0(elmt,"sunflowerseed")), 
                                                                        !!as.symbol(paste0(elmt,"seed_eq_sunflowerseed_oil")), 
                                                                        !!as.symbol(paste0(elmt,"seed_eq_sunflowerseed_cake"))), na.rm = TRUE) )
  
  # and add oil crops up 
  pfb <- mutate(pfb, !!as.symbol(paste0(elmt,"oilfeed_crops")) := sum(c(!!as.symbol(paste0(elmt,"groundnut_compo")), 
                                                                        !!as.symbol(paste0(elmt,"rapeseed_compo")), 
                                                                        !!as.symbol(paste0(elmt,"sunflowerseed_compo"))), na.rm = TRUE) )
  
  # source for oil and cake is https://cottonaustralia.com.au/uses-of-cotton#:~:text=single%20bed%20sheets-,Cotton%20seed,meal%20and%20300kg%20of%20hulls.
  pfb[, paste0(elmt,"seed_eq_cottonseed_oil")] <-   pfb[, paste0(elmt,"cottonseed_oil")]*(1/0.2)
  pfb[, paste0(elmt,"seed_eq_cottonseed_cake")] <-   pfb[, paste0(elmt,"cottonseed_cake")]*(1/(0.5))
  pfb[, paste0(elmt,"seed_eq_cotton_lint")] <-   pfb[, paste0(elmt,"cotton_lint")]*(1/(0.33)) # this is the conversion factor from GAEZ
  pfb <- mutate(pfb, !!as.symbol(paste0(elmt,"cotton")) := sum(c(!!as.symbol(paste0(elmt,"cottonseed")), 
                                                                 !!as.symbol(paste0(elmt,"seed_eq_cottonseed_oil")), 
                                                                 !!as.symbol(paste0(elmt,"seed_eq_cottonseed_cake")), 
                                                                 !!as.symbol(paste0(elmt,"seed_eq_cotton_lint"))), na.rm = TRUE) )
  
  # And for cereals and cocoa and coffee, weadd up directly
  pfb <- mutate(pfb, !!as.symbol(paste0(elmt,"cereals")) := sum(c(!!as.symbol(paste0(elmt,"barley")), 
                                                                  !!as.symbol(paste0(elmt,"sorghum")), 
                                                                  !!as.symbol(paste0(elmt,"wheat"))), na.rm = TRUE) )
  
  pfb <- mutate(pfb, !!as.symbol(paste0(elmt,"cocoa_coffee")) := sum(c(!!as.symbol(paste0(elmt,"cocoa")), 
                                                                       !!as.symbol(paste0(elmt,"coffee"))), na.rm = TRUE) )
}


#### MAKE STATISTICS OF INTEREST ---------------------------------------------------------------------------------

# Exposure of selected crops
for(item in all_items[final_items]){
  pfb[, paste0("trade_expo_",item)] <- ( pfb[, paste0("export_ktonnes_",item)] + pfb[, paste0("import_ktonnes_",item)] ) /
                                                          pfb[, paste0("production_ktonnes_",item)]
  
  pfb[, paste0("export_expo_",item)] <- ( pfb[, paste0("export_ktonnes_",item)] ) /
                                          pfb[, paste0("production_ktonnes_",item)]
  
}




### AVERAGE OVER YEARS 
# We average over different sets of years, as it is not clear in advance what is the most appropriate period (there are trade-offs)
pretreat_year_sets <- list(`2001_2007` = c(2001:2007), 
                           `2004_2007` = c(2004:2007),
                           `2006_2007` = c(2006:2007))


csfb <- list()
for(pretreat_period in pretreat_year_sets){
  
  csfb <- ddply(pfb[pfb$year %in% c(pretreat_period, 2010, 2011), ], # the 2010-2011 part is to include countries missing in old FAO data
                "country", summarise, 
                

                
                # Main dependency variables: 
                dependency_calorie_total = mean(dependency_calorie_total, na.rm = TRUE), 
                dependency_protein_total = mean(dependency_protein_total, na.rm = TRUE), 
                dependency_fat_total = mean(dependency_fat_total, na.rm = TRUE), 

                # import, total: 
                import_kcal_total = mean(import_kcal_total, na.rm = TRUE), 
                import_gprot_total = mean(import_gprot_total, na.rm = TRUE), 
                import_gfat_total = mean(import_gfat_total, na.rm = TRUE), 
                # gross supply, total: 
                gross_supply_kcal_total = mean(gross_supply_kcal_total, na.rm = TRUE), 
                gross_supply_gprot_total = mean(gross_supply_gprot_total, na.rm = TRUE), 
                gross_supply_gfat_total = mean(gross_supply_gfat_total, na.rm = TRUE), 
                
                # Dependency through specific crops, for raw weight, calorie, protein, and fat
                !!as.symbol("dependency_cereals") := mean(!!as.symbol("dependency_cereals"), na.rm = TRUE), 
                !!as.symbol("dependency_wheat") := mean(!!as.symbol("dependency_wheat"), na.rm = TRUE), 
                !!as.symbol("dependency_rice") := mean(!!as.symbol("dependency_rice"), na.rm = TRUE), 
                !!as.symbol("dependency_barley") := mean(!!as.symbol("dependency_barley"), na.rm = TRUE), 
                !!as.symbol("dependency_maize") := mean(!!as.symbol("dependency_maize"), na.rm = TRUE), 
                !!as.symbol("dependency_oilcrops") := mean(!!as.symbol("dependency_oilcrops"), na.rm = TRUE), 
                !!as.symbol("dependency_soybeans") := mean(!!as.symbol("dependency_soybeans"), na.rm = TRUE), 
                !!as.symbol("dependency_vegetable_oils") := mean(!!as.symbol("dependency_vegetable_oils"), na.rm = TRUE), 
                !!as.symbol("dependency_soybean_oil") := mean(!!as.symbol("dependency_soybean_oil"), na.rm = TRUE), 
                !!as.symbol("dependency_sunflowerseed_oil") := mean(!!as.symbol("dependency_sunflowerseed_oil"), na.rm = TRUE), 
                !!as.symbol("dependency_rapeseed_oil") := mean(!!as.symbol("dependency_rapeseed_oil"), na.rm = TRUE), 
                !!as.symbol("dependency_palm_oil") := mean(!!as.symbol("dependency_palm_oil"), na.rm = TRUE),
                !!as.symbol("dependency_meat") := mean(!!as.symbol("dependency_meat"), na.rm = TRUE),
                !!as.symbol("dependency_fish") := mean(!!as.symbol("dependency_fish"), na.rm = TRUE),
                
                !!as.symbol("dependency_calorie_cereals") := mean(!!as.symbol("dependency_calorie_cereals"), na.rm = TRUE), 
                !!as.symbol("dependency_calorie_wheat") := mean(!!as.symbol("dependency_calorie_wheat"), na.rm = TRUE), 
                !!as.symbol("dependency_calorie_rice") := mean(!!as.symbol("dependency_calorie_rice"), na.rm = TRUE), 
                !!as.symbol("dependency_calorie_barley") := mean(!!as.symbol("dependency_calorie_barley"), na.rm = TRUE), 
                !!as.symbol("dependency_calorie_maize") := mean(!!as.symbol("dependency_calorie_maize"), na.rm = TRUE), 
                !!as.symbol("dependency_calorie_oilcrops") := mean(!!as.symbol("dependency_calorie_oilcrops"), na.rm = TRUE), 
                !!as.symbol("dependency_calorie_soybeans") := mean(!!as.symbol("dependency_calorie_soybeans"), na.rm = TRUE), 
                !!as.symbol("dependency_calorie_vegetable_oils") := mean(!!as.symbol("dependency_calorie_vegetable_oils"), na.rm = TRUE), 
                !!as.symbol("dependency_calorie_soybean_oil") := mean(!!as.symbol("dependency_calorie_soybean_oil"), na.rm = TRUE), 
                !!as.symbol("dependency_calorie_sunflowerseed_oil") := mean(!!as.symbol("dependency_calorie_sunflowerseed_oil"), na.rm = TRUE), 
                !!as.symbol("dependency_calorie_rapeseed_oil") := mean(!!as.symbol("dependency_calorie_rapeseed_oil"), na.rm = TRUE), 
                !!as.symbol("dependency_calorie_palm_oil") := mean(!!as.symbol("dependency_calorie_palm_oil"), na.rm = TRUE),
                !!as.symbol("dependency_calorie_meat") := mean(!!as.symbol("dependency_calorie_meat"), na.rm = TRUE),
                !!as.symbol("dependency_calorie_fish") := mean(!!as.symbol("dependency_calorie_fish"), na.rm = TRUE),
                
                !!as.symbol("dependency_protein_cereals") := mean(!!as.symbol("dependency_protein_cereals"), na.rm = TRUE), 
                !!as.symbol("dependency_protein_wheat") := mean(!!as.symbol("dependency_protein_wheat"), na.rm = TRUE), 
                !!as.symbol("dependency_protein_rice") := mean(!!as.symbol("dependency_protein_rice"), na.rm = TRUE), 
                !!as.symbol("dependency_protein_barley") := mean(!!as.symbol("dependency_protein_barley"), na.rm = TRUE), 
                !!as.symbol("dependency_protein_maize") := mean(!!as.symbol("dependency_protein_maize"), na.rm = TRUE), 
                !!as.symbol("dependency_protein_oilcrops") := mean(!!as.symbol("dependency_protein_oilcrops"), na.rm = TRUE), 
                !!as.symbol("dependency_protein_soybeans") := mean(!!as.symbol("dependency_protein_soybeans"), na.rm = TRUE), 
                !!as.symbol("dependency_protein_vegetable_oils") := mean(!!as.symbol("dependency_protein_vegetable_oils"), na.rm = TRUE), 
                !!as.symbol("dependency_protein_soybean_oil") := mean(!!as.symbol("dependency_protein_soybean_oil"), na.rm = TRUE), 
                !!as.symbol("dependency_protein_sunflowerseed_oil") := mean(!!as.symbol("dependency_protein_sunflowerseed_oil"), na.rm = TRUE), 
                !!as.symbol("dependency_protein_rapeseed_oil") := mean(!!as.symbol("dependency_protein_rapeseed_oil"), na.rm = TRUE), 
                !!as.symbol("dependency_protein_palm_oil") := mean(!!as.symbol("dependency_protein_palm_oil"), na.rm = TRUE),
                !!as.symbol("dependency_protein_meat") := mean(!!as.symbol("dependency_protein_meat"), na.rm = TRUE),
                !!as.symbol("dependency_protein_fish") := mean(!!as.symbol("dependency_protein_fish"), na.rm = TRUE),
                
                !!as.symbol("dependency_fat_cereals") := mean(!!as.symbol("dependency_fat_cereals"), na.rm = TRUE), 
                !!as.symbol("dependency_fat_wheat") := mean(!!as.symbol("dependency_fat_wheat"), na.rm = TRUE), 
                !!as.symbol("dependency_fat_rice") := mean(!!as.symbol("dependency_fat_rice"), na.rm = TRUE), 
                !!as.symbol("dependency_fat_barley") := mean(!!as.symbol("dependency_fat_barley"), na.rm = TRUE), 
                !!as.symbol("dependency_fat_maize") := mean(!!as.symbol("dependency_fat_maize"), na.rm = TRUE), 
                !!as.symbol("dependency_fat_oilcrops") := mean(!!as.symbol("dependency_fat_oilcrops"), na.rm = TRUE), 
                !!as.symbol("dependency_fat_soybeans") := mean(!!as.symbol("dependency_fat_soybeans"), na.rm = TRUE), 
                !!as.symbol("dependency_fat_vegetable_oils") := mean(!!as.symbol("dependency_fat_vegetable_oils"), na.rm = TRUE), 
                !!as.symbol("dependency_fat_soybean_oil") := mean(!!as.symbol("dependency_fat_soybean_oil"), na.rm = TRUE), 
                !!as.symbol("dependency_fat_sunflowerseed_oil") := mean(!!as.symbol("dependency_fat_sunflowerseed_oil"), na.rm = TRUE), 
                !!as.symbol("dependency_fat_rapeseed_oil") := mean(!!as.symbol("dependency_fat_rapeseed_oil"), na.rm = TRUE), 
                !!as.symbol("dependency_fat_palm_oil") := mean(!!as.symbol("dependency_fat_palm_oil"), na.rm = TRUE),
                !!as.symbol("dependency_fat_meat") := mean(!!as.symbol("dependency_fat_meat"), na.rm = TRUE),
                !!as.symbol("dependency_fat_fish") := mean(!!as.symbol("dependency_fat_fish"), na.rm = TRUE),
                
                # import of specific commodities 
                !!as.symbol("import_kcal_cereals") := mean(!!as.symbol("import_kcal_cereals"), na.rm = TRUE), 
                !!as.symbol("import_kcal_wheat") := mean(!!as.symbol("import_kcal_wheat"), na.rm = TRUE), 
                !!as.symbol("import_kcal_rice") := mean(!!as.symbol("import_kcal_rice"), na.rm = TRUE), 
                !!as.symbol("import_kcal_barley") := mean(!!as.symbol("import_kcal_barley"), na.rm = TRUE), 
                !!as.symbol("import_kcal_maize") := mean(!!as.symbol("import_kcal_maize"), na.rm = TRUE), 
                !!as.symbol("import_kcal_oilcrops") := mean(!!as.symbol("import_kcal_oilcrops"), na.rm = TRUE), 
                !!as.symbol("import_kcal_soybeans") := mean(!!as.symbol("import_kcal_soybeans"), na.rm = TRUE), 
                !!as.symbol("import_kcal_vegetable_oils") := mean(!!as.symbol("import_kcal_vegetable_oils"), na.rm = TRUE), 
                !!as.symbol("import_kcal_soybean_oil") := mean(!!as.symbol("import_kcal_soybean_oil"), na.rm = TRUE), 
                !!as.symbol("import_kcal_sunflowerseed_oil") := mean(!!as.symbol("import_kcal_sunflowerseed_oil"), na.rm = TRUE), 
                !!as.symbol("import_kcal_rapeseed_oil") := mean(!!as.symbol("import_kcal_rapeseed_oil"), na.rm = TRUE), 
                !!as.symbol("import_kcal_palm_oil") := mean(!!as.symbol("import_kcal_palm_oil"), na.rm = TRUE),
                !!as.symbol("import_kcal_meat") := mean(!!as.symbol("import_kcal_meat"), na.rm = TRUE),
                !!as.symbol("import_kcal_fish") := mean(!!as.symbol("import_kcal_fish"), na.rm = TRUE),
                
                # gross supply of specific commodities 
                !!as.symbol("gross_supply_kcal_cereals") := mean(!!as.symbol("gross_supply_kcal_cereals"), na.rm = TRUE), 
                !!as.symbol("gross_supply_kcal_wheat") := mean(!!as.symbol("gross_supply_kcal_wheat"), na.rm = TRUE), 
                !!as.symbol("gross_supply_kcal_rice") := mean(!!as.symbol("gross_supply_kcal_rice"), na.rm = TRUE), 
                !!as.symbol("gross_supply_kcal_barley") := mean(!!as.symbol("gross_supply_kcal_barley"), na.rm = TRUE), 
                !!as.symbol("gross_supply_kcal_maize") := mean(!!as.symbol("gross_supply_kcal_maize"), na.rm = TRUE), 
                !!as.symbol("gross_supply_kcal_oilcrops") := mean(!!as.symbol("gross_supply_kcal_oilcrops"), na.rm = TRUE), 
                !!as.symbol("gross_supply_kcal_soybeans") := mean(!!as.symbol("gross_supply_kcal_soybeans"), na.rm = TRUE), 
                !!as.symbol("gross_supply_kcal_vegetable_oils") := mean(!!as.symbol("gross_supply_kcal_vegetable_oils"), na.rm = TRUE), 
                !!as.symbol("gross_supply_kcal_soybean_oil") := mean(!!as.symbol("gross_supply_kcal_soybean_oil"), na.rm = TRUE), 
                !!as.symbol("gross_supply_kcal_sunflowerseed_oil") := mean(!!as.symbol("gross_supply_kcal_sunflowerseed_oil"), na.rm = TRUE), 
                !!as.symbol("gross_supply_kcal_rapeseed_oil") := mean(!!as.symbol("gross_supply_kcal_rapeseed_oil"), na.rm = TRUE), 
                !!as.symbol("gross_supply_kcal_palm_oil") := mean(!!as.symbol("gross_supply_kcal_palm_oil"), na.rm = TRUE),
                !!as.symbol("gross_supply_kcal_meat") := mean(!!as.symbol("gross_supply_kcal_meat"), na.rm = TRUE),
                !!as.symbol("gross_supply_kcal_fish") := mean(!!as.symbol("gross_supply_kcal_fish"), na.rm = TRUE),
                
                # conversion factors, for some major items: 
                !!as.symbol("kcal_per_kg_cereals") := mean(!!as.symbol("kcal_per_kg_cereals"), na.rm = TRUE), 
                !!as.symbol("gprot_per_kg_cereals") := mean(!!as.symbol("gprot_per_kg_cereals"), na.rm = TRUE), 
                !!as.symbol("gfat_per_kg_cereals") := mean(!!as.symbol("gfat_per_kg_cereals"), na.rm = TRUE), 
                !!as.symbol("kcal_per_kg_oilcrops") := mean(!!as.symbol("kcal_per_kg_oilcrops"), na.rm = TRUE), 
                !!as.symbol("kcal_per_kg_vegetable_oils") := mean(!!as.symbol("kcal_per_kg_vegetable_oils"), na.rm = TRUE), 
                !!as.symbol("gfat_per_kg_vegetable_oils") := mean(!!as.symbol("gfat_per_kg_vegetable_oils"), na.rm = TRUE),
                !!as.symbol("gprot_per_kg_meat") := mean(!!as.symbol("gprot_per_kg_meat"), na.rm = TRUE), 
                !!as.symbol("gprot_per_kg_fish") := mean(!!as.symbol("gprot_per_kg_fish"), na.rm = TRUE) 
  )
  
  
  saveRDS(csfb, file = here("temp_data", "exposures", paste0("dependency_",
                                                             min(pretreat_period),
                                                             max(pretreat_period),".Rdata")))
  
}



rm(fb, pfb, pfbsave, rdc, prdc, wide_fb, wide_fb_list, wide_rdc, wide_rdc_list, csfb, long_ds, avg_convfact, missing_cols, elmt_wide_ds_list)




#### Data exploration  -----------------------------------------------------------------------------------------------
# Don't erase this part, it is not as complete in / it's different in anaylses_rfsFOOD.R
sffb <- st_read(here("input_data", "Global_LSIB_Polygons_Detailed"))

# necessary for simplifying below, and better for unioning 
sffb <- st_transform(sffb, crs = 4088)

names(sffb)[names(sffb) == "COUNTRY_NA"] <- "country"
# sffb[grepl("ongo", sffb$country), ]
# pfb[grepl("ongo", pfb$country), "country"]

### Match country names to those from FAOSTAT 
## /!\ NEXT TIME CONSIDER USING A FAOSTAT ADMIN BOUNDARIES MAP DIRECTLY ! 
# (but some advantages to adjust manually some things here)

# Note that some FAOSTAT names have been modified above already (they are on top here)
# Note also that we do not associate islands to their main country, for simplicity (by commenting their section out)
sffb$country[sffb$country=="Cote d'Ivoire"] <- "Ivory Coast" # NOTICE THIS (FAOSTAT name is Côte d'Ivoire)
sffb$country[sffb$country=="Congo, Rep of the"] <- "Republic of the Congo" # This name is not from FAOSTAT, I changed it above

# Political discrepancies
sffb$country[sffb$country=="China"] <- "China, mainland" # NOTICE THIS 
sffb$country[sffb$country=="Hong Kong (Ch)"] <- "China, Hong Kong SAR"
sffb$country[sffb$country=="Macau (Ch)"] <- "China, Macao SAR"
# sffb$country[sffb$country=="Taiwan"] 

# for Sudan and South Sudan, they will be splitted in the anaysis period. They will have outcome data separately. 
# But they will have exactly the same treatment. So we can treat them as within the same cluster, 
# or merge them into a single country, with average outcome. 
sffb[sffb$country=="Sudan", "geometry"] <- st_union(sffb[sffb$country == "Sudan", "geometry"], 
                                                    sffb[sffb$country == "South Sudan", "geometry"]) # %>% dplyr::select(geometry)
sffb$country[sffb$country=="Sudan"] <- "Sudan (former)" # South Sudan independent since 2005, after our data period

sffb[sffb$country=="Serbia", "geometry"] <- st_union(sffb[sffb$country == "Serbia", "geometry"], 
                                                     sffb[sffb$country == "Montenegro", "geometry"]) # %>% dplyr::select(geometry)
sffb$country[sffb$country=="Serbia"] <- "Serbia and Montenegro" # Montenegro independent since 2006, after our data period

# Countries for which it's simply a matter of different way to write the name down
sffb$country[sffb$country=="Congo, Dem Rep of the"] <- "Democratic Republic of the Congo"
sffb$country[sffb$country=="Macedonia"] <- "North Macedonia"
sffb$country[sffb$country=="United Kingdom"] <- "United Kingdom of Great Britain and Northern Ireland"
sffb$country[sffb$country=="United States"] <- "United States of America"
sffb$country[sffb$country=="Russia"] <- "Russian Federation"
sffb$country[sffb$country=="Syria"] <- "Syrian Arab Republic"
sffb$country[sffb$country=="Bosnia & Herzegovina"] <- "Bosnia and Herzegovina"
sffb$country[sffb$country=="Korea, North"] <- "Democratic People's Republic of Korea"
sffb$country[sffb$country=="Korea, South"] <- "Republic of Korea"
sffb$country[sffb$country=="Moldova"] <- "Republic of Moldova"
sffb$country[sffb$country=="Gambia, The"] <- "Gambia"
sffb$country[sffb$country=="Swaziland"] <- "Eswatini"
sffb$country[sffb$country=="Iran"] <- "Iran (Islamic Republic of)"
sffb$country[sffb$country=="Burma"] <- "Myanmar"
sffb$country[sffb$country=="Bahamas, The"] <- "Bahamas"
sffb$country[sffb$country=="Vietnam"] <- "Viet Nam"
sffb$country[sffb$country=="Laos"] <- "Lao People's Democratic Republic"
sffb$country[sffb$country=="Antigua & Barbuda"] <- "Antigua and Barbuda"
sffb$country[sffb$country=="St Kitts & Nevis"] <- "Saint Kitts and Nevis"
sffb$country[sffb$country=="St Lucia"] <- "Saint Lucia"
sffb$country[sffb$country=="St Vincent & the Grenadines"] <- "Saint Vincent and the Grenadines"
sffb$country[sffb$country=="Venezuela"] <- "Venezuela (Bolivarian Republic of)"
sffb$country[sffb$country=="Trinidad & Tobago"] <- "Trinidad and Tobago"
sffb$country[sffb$country=="Central African Rep"] <- "Central African Republic"
sffb$country[sffb$country=="Brunei"] <- "Brunei Darussalam"
sffb$country[sffb$country=="Sao Tome & Principe"] <- "Sao Tome and Principe"
sffb$country[sffb$country=="Tanzania"] <- "United Republic of Tanzania"
sffb$country[sffb$country=="Solomon Is"] <- "Solomon Islands"
sffb$country[sffb$country=="Bolivia"] <- "Bolivia (Plurinational State of)"
# sffb$country[sffb$country=="Micronesia, Fed States of"]

# ISLANDS
# sffb$country[sffb$country=="Niue (NZ)"] <- "New Zealand"
# sffb$country[sffb$country=="Turks & Caicos Is (UK)"] <- "United Kingdom of Great Britain and Northern Ireland"
# sffb$country[sffb$country=="Montserrat (UK)"] <- "United Kingdom of Great Britain and Northern Ireland"
# sffb$country[sffb$country=="Br Virgin Is (UK)"] <- "United Kingdom of Great Britain and Northern Ireland"
# sffb$country[sffb$country=="Pitcairn Is (UK)"] <- "United Kingdom of Great Britain and Northern Ireland"
# sffb$country[sffb$country=="Cayman Is (UK)"] <- "United Kingdom of Great Britain and Northern Ireland"
# sffb$country[sffb$country=="Puerto Rico (US)"] <- "United States of America"
# sffb$country[sffb$country=="Northern Mariana Is (US)"] <- "United States of America"
# sffb$country[sffb$country=="American Samoa (US)"] <- "United States of America"
# sffb$country[sffb$country=="Guam (US)"] <- "United States of America"
# sffb$country[sffb$country=="US Virgin Is (US)"] <- "United States of America"
# sffb$country[sffb$country=="Christmas I (Aus)"] <- "Australia"
# sffb$country[sffb$country=="Spain [Canary Is]"] <- "Spain"
# sffb$country[sffb$country=="Guadeloupe (Fr)"] <- "France"
# sffb$country[sffb$country=="Martinique (Fr)"] <- "France"
# sffb$country[sffb$country=="French Guiana (Fr)"] <- "France"
# sffb$country[sffb$country=="French Polynesia (Fr)"] <- "France"
# sffb$country[sffb$country=="Mayotte (Fr)"] <- "France"
# sffb$country[sffb$country=="Wallis & Futuna (Fr)"] <- "France"
# sffb$country[sffb$country=="New Caledonia (Fr)"] <- "France"
# sffb$country[sffb$country=="Reunion (Fr)"] <- "France"

# length(unique(sffb$country))
# sffb$country[duplicated(sffb$country)]
# length(unique(csfb5$country))
# sfctry <- unique(sffb$country)
# fbctry <- unique(csfb5$country)
# extra_names <- sfctry[!(sfctry %in% fbctry)]
# # (all names in rdc data are in global data)
# fbctry[!(fbctry %in% sfctry)]
# sfctry[grepl("yri", sfctry)]
# fbctry[grepl("apua", fbctry)]
pfbnames[grepl("uantity", pfbnames)]


# simplify before plotting
sffb$geometry <- st_simplify(sffb$geometry, dTolerance = 1000)

sffb <- st_transform(sffb, crs = 4326)

csfb5_sf <- left_join(csfb5, sffb, by = "country") %>% st_as_sf()

csfb5_sf$country[st_is_empty(csfb5_sf$geometry)]

#### PLOT TOTAL DEPENDENCY #### 
pal_dep <- colorNumeric("viridis", # "viridis" (green-purple), "magma" (yellow-purple), "inferno" (like magma), or "plasma", "BuPu", "Greens"
                        domain = st_drop_geometry(csfb5_sf[,"dependency_calorie_total"]),
                        #bins = 4, 
                        na.color = "transparent", 
                        reverse = F)


# popup
csfb5_sf$popup_total <- paste0(csfb5_sf$country, "<br/>",
                               "Imports: ", formatC(csfb5_sf$import_kcal_total/1e6, format = "e", digits = 2), " bn cal", "<br/>",
                               "Supply + Exports: ", formatC(csfb5_sf$gross_supply_kcal_total/1e6,format = "e", digits = 2), " bn cal", "<br/>"
)



leaflet() %>% 
  addTiles()%>%
  addProviderTiles(providers$Esri.WorldGrayCanvas, group ="ESRI") %>%
  setView(lat = 0, 
          lng = 0, 
          zoom = 1) %>% 
  addPolygons(data = csfb5_sf, 
              opacity = 0, color = "black", weight = 2, 
              fill = TRUE, fillColor = ~pal_dep(csfb5_sf$dependency_calorie_total), fillOpacity = 0.5,
              popup = ~csfb5_sf$popup_total, 
              popupOptions = popupOptions(riseOnHover = TRUE, 
                                          bringToFront = TRUE),
              highlightOptions = highlightOptions(bringToFront = TRUE)
  ) %>% 
  # addMarkers(data = csfb5_sf, 
  #            popup = ~csfb5_sf$popup_total,
  #            options = markerOptions(riseOnHover = TRUE)) %>% 
  addLegend(pal = pal_dep,  
            values = csfb5_sf$dependency_calorie_total, 
            bins = 5, opacity = 0.4,
            title = "Import dependency <br/> for calories, total",
            position = "bottomright") 



#### PLOT CONVERSION FACTORS #### 
pal_dep <- colorNumeric("viridis", # "viridis" (green-purple), "magma" (yellow-purple), "inferno" (like magma), or "plasma", "BuPu", "Greens"
                        domain = st_drop_geometry(csfb5_sf[,"kcal_per_kg_cereals"]),
                        #bins = 4, 
                        na.color = "transparent", 
                        reverse = F)


# popup
csfb5_sf$popup_convfact <- paste0(csfb5_sf$country, "<br/>",
                                  "kcal per kg cereals (excl. beer): ", round(csfb5_sf$kcal_per_kg_cereals,1), "<br/>",
                                  "g protein per kg cereals (excl. beer): ", round(csfb5_sf$gprot_per_kg_cereals,1), "<br/>",
                                  "g fat per kg cereals (excl. beer): ", round(csfb5_sf$gfat_per_kg_cereals,1), "<br/>",
                                  "kcal per kg vegetable oils: ", round(csfb5_sf$kcal_per_kg_vegetable_oils,1), "<br/>",
                                  "gprot per kg meat: ", round(csfb5_sf$gprot_per_kg_meat,1), "<br/>",
                                  "gprot per kg fish & seafood: ", round(csfb5_sf$gprot_per_kg_fish,1), "<br/>"
)

leaflet() %>% 
  addTiles()%>%
  addProviderTiles(providers$Esri.WorldGrayCanvas, group ="ESRI") %>%
  setView(lat = 0, 
          lng = 0, 
          zoom = 1) %>% 
  addPolygons(data = csfb5_sf, 
              opacity = 0, color = "black", weight = 2, 
              fill = TRUE, fillColor = ~pal_dep(csfb5_sf$kcal_per_kg_cereals), fillOpacity = 0.5,
              popup = ~csfb5_sf$popup_convfact, 
              popupOptions = popupOptions(riseOnHover = TRUE, 
                                          bringToFront = TRUE),
              highlightOptions = highlightOptions(bringToFront = TRUE)) %>% 
  # addMarkers(data = csfb5_sf, 
  #            popup = ~csfb5_sf$popup_convfact,
  #            options = markerOptions(riseOnHover = TRUE)) %>% 
  addLegend(pal = pal_dep,  
            values = csfb5_sf$kcal_per_kg_cereals, 
            bins = 5, opacity = 0.4,
            title = "Nutritional density <br/> of cereals (kcal/kg) <br/> (2001-2005)",
            position = "bottomright") 

plot(csfb5_sf[,"dependency_calorie_total"])
plot(csfb5_sf[,"dependency_protein_total"])
plot(csfb5_sf[,"gross_supply_kcal_total"])
plot(csfb5_sf[,"import_kcal_total"])


csfb5_sf[csfb5_sf$country=="Niger", c("dependency_calorie_total", "geometry")] %>% plot()










