### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages = c("plyr", "dplyr", "readxl", "foreign", "here", 
                   "DataCombine") 
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

### Script steps description ### 
# Import and prepare Pink Sheet data. Price annual time series of most crops are taken from there.   
# Import time series from IMF, "manually" extract series for missing crops (those in GAEZ but not in Pink Sheet) and prepare
# 
# 
# 



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# PREPARE INTERNATIONAL PRICE TIME SERIES 

#### PINK SHEET DATA (WORLD BANK) ####
# "Annual Prices" link at: https://www.worldbank.org/en/research/commodity-markets (in April 2021)
ps <- read_xlsx(here("input_data", "price_data", "CMO-Historical-Data-Annual.xlsx"), 
                sheet = "Annual Prices (Real)", 
                range = "A7:BS69", 
                col_names =TRUE)

names(ps)[1] <- "year"

ps <- as.data.frame(ps)

# Convert $/kg in $/mt (the warnings are for columns that have some NAs)
ps[-1,] <- ps[-1,] %>% mutate(across(.cols = colnames(ps[,ps[1,]=="($/kg)" & !is.na(ps[1,])]),
                           .fns = function(c){as.numeric(c)*1000}))
                                               
# drop first row (price unit) for manipulation convenience. 
ps <- ps[-1,]
# Convert to numeric
ps <- summarise(ps, across(.fns = as.numeric))


# Select commodities of interest 
ps_commo <- c("Banana, US", 
              "Barley",
              "Beef", 
              "Meat, chicken", 
              "Meat, sheep",
              "Crude oil, average",
              "Orange", 
              "Cocoa", 
              "Coconut oil", 
              "Coffee, Arabica",
              "Coffee, Robusta", 
              "Cotton, A Index",
              "Rice, Thai 5%",
              "Groundnuts", 
              "Maize",
              "Palm oil", 
              "Rubber, SGP/MYS",
              "Sorghum", 
              "Soybeans", 
              "Soybean oil",
              "Soybean meal",
              "Sugar, world", 
              "Tea, avg 3 auctions", 
              "Tobacco, US import u.v.", 
              "Wheat, US HRW")

ps2 <- ps[, c("year", ps_commo, "MUV Index")] 



# # Average prices of the two coffee types arabica and robusta, to better match coffee price from FAOSTAT
# ps2 <-  dplyr::mutate(ps2, Coffee = rowMeans(across(.cols = all_of(c("Coffee, Arabica", "Coffee, Robusta")))))
# ps2 <- dplyr::select(ps2, -'Coffee, Arabica', - 'Coffee, Robusta')

# Change (simplify) some names
names(ps2)[names(ps2) == "Banana, US"] <- "Banana"
names(ps2)[names(ps2) == "Coconut oil"] <- "Coconut_oil"
names(ps2)[names(ps2) == "Coffee, Arabica"] <- "Coffee"
names(ps2)[names(ps2) == "Cotton, A Index"] <- "Cotton"
names(ps2)[names(ps2) == "Crude oil, average"] <- "Crude_oil"
names(ps2)[names(ps2) == "Groundnuts"] <- "Groundnut"
names(ps2)[names(ps2) == "Meat, chicken"] <- "Chicken"
names(ps2)[names(ps2) == "Meat, sheep"] <- "Sheep"
names(ps2)[names(ps2) == "Rice, Thai 5%"] <- "Rice"
names(ps2)[names(ps2) == "Palm oil"] <- "Palm_oil"
names(ps2)[names(ps2) == "Rubber, SGP/MYS"] <- "Rubber"
names(ps2)[names(ps2) == "Soybeans"] <- "Soybean"
names(ps2)[names(ps2) == "Soybean oil"] <- "Soybean_oil"
names(ps2)[names(ps2) == "Soybean meal"] <- "Soybean_meal"
names(ps2)[names(ps2) == "Sugar, world"] <- "Sugar"
names(ps2)[names(ps2) == "Tea, avg 3 auctions"] <- "Tea"
names(ps2)[names(ps2) == "Tobacco, US import u.v."] <- "Tobacco"
names(ps2)[names(ps2) == "Wheat, US HRW"] <- "Wheat"

names(ps2)[names(ps2) == "MUV Index"] <- "MUV_index"

head(ps2)

# Convert from base 2010 to base 2014-2016
# First construct index for 2014-2016, to match the FAOSTAT index
# So this is the MUV index for 2014-2016 in base 2010=100
muv_index_2014_16 <- mean(c(as.numeric(ps2[ps2$year >= 2014 & ps2$year <= 2016, "MUV_index"])))

# Important to not modify MUV_index, we will need it again as such
ps2 <- dplyr::mutate(ps2, across(.cols = (!contains("year") & !contains("MUV_index")), 
                                 .fns = ~.*muv_index_2014_16/100)) 




#### IMF DATA #### 
# "Excel Database" link at: webpage: https://www.imf.org/en/Research/commodity-prices (in April 2021)
imf <- read_xls(here("input_data", "price_data", "external-dataAPR.xls"), col_names = TRUE)

length(unique(colnames(imf))) == ncol(imf)
names(imf)[1] <- "year"
# Select columns we need 
colnames(imf)
needed_imf_col <- c("year", 
                    colnames(imf)[grepl(pattern = "oat", x = imf[1,], ignore.case = TRUE)], 
                    colnames(imf)[grepl(pattern = "olive", x = imf[1,], ignore.case = TRUE)], 
                    colnames(imf)[grepl(pattern = "milk", x = imf[1,], ignore.case = TRUE)], 
                    colnames(imf)[grepl(pattern = "rapeseed", x = imf[1,], ignore.case = TRUE)], 
                    colnames(imf)[grepl(pattern = "sunflower", x = imf[1,], ignore.case = TRUE)],
                    # colnames(imf)[grepl(pattern = "arabica", x = imf[1,], ignore.case = TRUE)],
                    colnames(imf)[grepl(pattern = "pork", x = imf[1,], ignore.case = TRUE)],
                    # PFOOD is the name of the column for the food price index 
                    colnames(imf)[grepl(pattern = "pfood", x = imf[1,], ignore.case = TRUE)]) 
imf <- imf[, needed_imf_col]

# drop first, unnecessary rows
imf <- imf[-c(1,2,3),]

imf <- as.data.frame(imf)

# Rename /!\ ORDER MATTERS HERE ! 
names(imf)
names(imf) <- c("year", "Oat", "Olive_oil", "FPI", "Milk", "Rapeseed_oil", "Sunflower_oil",  "Pork")#"Coffee",

# Convert in $/mt 
# Oat is in USD/bushel (ratio from https://www.sagis.org.za/conversion_table.html)
imf[,"Oat"] <- as.numeric(imf[,"Oat"]) / 0.014515 
# Milk is in USD/cwt i.e. hundredweight or 45.36kg  https://en.wikipedia.org/wiki/Hundredweight
imf[,"Milk"] <- as.numeric(imf[,"Milk"]) * 45.36 / 1000 
# Pork is in USD cents / pound (ratio from https://www.rapidtables.com/convert/weight/pound-to-kg.html?x=1&x2=&x3=)
imf[,"Pork"] <- as.numeric(imf[,"Pork"])*0.01/0.000453592 

# imf[,"Coffee"] <- as.numeric(imf[,"Coffee"])*0.01/0.000453592 

# Average monthly prices to annual 
imf$year <- gsub("M.*", "", imf$year) %>% as.numeric()


imf <- ddply(imf, "year", summarise, 
                    Oat = Oat %>% as.numeric() %>% mean(), 
                    Olive_oil = Olive_oil %>% as.numeric() %>% mean(), 
                    Milk = Milk %>% as.numeric() %>% mean(), 
                    Rapeseed_oil = Rapeseed_oil %>% as.numeric() %>% mean(), 
                    Sunflower_oil = Sunflower_oil %>% as.numeric() %>% mean(), 
                    # Coffee = mean(Coffee), 
                    Pork = Pork %>% as.numeric() %>% mean(), 
                    FPI = FPI %>% as.numeric() %>% mean())

# From 1992 on, use the Food price index (FPI) from the IMF, in base 2016=100, to deflate the nominal time series.
imf[imf$year >=1992,] <- mutate(imf[imf$year >=1992,], Oat = as.numeric(Oat)/(FPI/100), 
                                                       Olive_oil = as.numeric(Olive_oil)/(FPI/100), 
                                                       Rapeseed_oil = as.numeric(Rapeseed_oil)/(FPI/100), 
                                                       Sunflower_oil = as.numeric(Sunflower_oil)/(FPI/100),
                                                       # Coffee = as.numeric(Coffee)/muv_index_2016,
                                                       Pork = as.numeric(Pork)/(FPI/100))

# Then change the base from 2016 to 2014-2016, to match FAOSTAT
fpi_2014_16 <- mean(c(as.numeric(imf[imf$year >= 2014 & imf$year <= 2016, "FPI"])))

imf <- dplyr::mutate(imf, across(.cols = (!contains("year")), 
                                 .fns = ~.*fpi_2014_16/100)) 

# For years prior 1992, we deflate with the MUV index base 2014-2016 from the Pink Sheet
imf$MUV_index_10 <- NA
imf$MUV_index_10[imf$year<1992] <- ps2$MUV_index[ps2$year >= min(imf$year) & ps2$year < 1992]

imf$MUV_index_2014_16 <- NA
imf$MUV_index_2014_16 <- imf$MUV_index_10/(muv_index_2014_16/100)

imf[imf$year < 1992,] <- mutate(imf[imf$year < 1992,], Oat = as.numeric(Oat)/(MUV_index_2014_16/100), 
                                Olive_oil = as.numeric(Olive_oil)/(MUV_index_2014_16/100), 
                                Rapeseed_oil = as.numeric(Rapeseed_oil)/(MUV_index_2014_16/100), 
                                Sunflower_oil = as.numeric(Sunflower_oil)/(MUV_index_2014_16/100),
                                # Coffee = as.numeric(Coffee)/(MUV_index_2014_16/100),
                                Pork = as.numeric(Pork)/(MUV_index_2014_16/100))


class(imf[,2])



#### Consolidate international price annual time series data frame #### 
# keep all years from each sources
ip <- full_join(x = ps2, y = imf, by = "year") 

ip <- ip %>% dplyr::select(-MUV_index, -FPI, -MUV_index_10, - MUV_index_2014_16)

### PREPARE ADDITIONAL INTERNATIONAL PRICE VARIABLES #### 

# ### Group some commodity prices (but not used currently)
# ip <- ip %>% rowwise() %>% mutate(cereal_crops = mean(c(Barley, Maize, Sorghum, Wheat, Oat), na.rm = T), # excluding rice as not comparable enough
#                                   oil_crops = mean(c(Palm_oil, Rapeseed_oil, Soybean_oil, Sunflower_oil), na.rm = T), # using only "unflavored" oils
# ) %>% as.data.frame() # other commodities are not comparable enough to be grouped.
# 
# head(ip)


# unique(ip$year) %>% length()

### Lags 
# Lagging here is made in preparation of USD international price time series for final use in analysis, not for merging with FAOSTAT. 
inter_prices <- ip
## Lag international prices over whole period 
ip_variables <- names(inter_prices)[names(inter_prices)!="year"]

for(voi in ip_variables){
  
  ## short to long lags
  for(lag in c(1:5)){
    inter_prices <- dplyr::arrange(inter_prices, year)
    inter_prices <- DataCombine::slide(inter_prices,
                                 Var = voi, 
                                 TimeVar = "year",
                                 NewVar = paste0(voi,"_lag",lag),
                                 slideBy = -lag, 
                                 keepInvalid = FALSE)
    inter_prices <- dplyr::arrange(inter_prices, year)
    
  }
  
  for(py in c(2:5)){
    ## Past-year averages (2, 3 and 4 years)  
    # note that we DON'T add voi column (not lagged) in the row mean
    inter_prices$newv <- rowMeans(x = inter_prices[,paste0(voi,"_lag",c(1:py))], na.rm = FALSE)
    inter_prices[is.nan(inter_prices$newv),"newv"] <- NA
    colnames(inter_prices)[colnames(inter_prices)=="newv"] <- paste0(voi,"_",py,"pya")
  }
}

# remove some variables that were only temporarily necessary
# (we want to keep lag1)
vars_torm <- names(inter_prices)[(grepl(pattern = "_lag2", x = names(inter_prices)) |
                              grepl(pattern = "_lag3", x = names(inter_prices)) |
                              grepl(pattern = "_lag4", x = names(inter_prices)) |
                              grepl(pattern = "_lag5", x = names(inter_prices)))]

inter_prices <- inter_prices[,!(names(inter_prices) %in% vars_torm)]

### Make logarithms
inter_prices <- dplyr::mutate(inter_prices, across(.cols = !c("year"),
                                       .fns = log,
                                       .names = paste0("ln_", "{.col}")))

saveRDS(inter_prices, here("temp_data", "prepared_international_prices.Rdata"))

rm(inter_prices)






#### FAOSTAT ####

# fao <- read.csv(file = here("input_data", "price_data", "FAOSTAT_data_6-11-2021.csv"))

fao <- read.csv(file = here("input_data", "price_data", "FAOSTAT_producerprice_slc_6-14-2021.csv"))
fao <- fao %>% dplyr::arrange(Area, Year, Item)
names(fao)[names(fao)=="Value"] <- "ppslc"

# Deflate PP-SLC from nominal values to real, base 2014-2016, values.
idx <- read.csv(file = here("input_data", "price_data", "FAOSTAT_producerprice_index_6-14-2021.csv"))
idx <- idx %>% dplyr::arrange(Area, Year, Item)
names(idx)[names(idx)=="Value"] <- "ppi"

unique(idx$Flag.Description)
idx <- dplyr::select(idx, Area, Year, Item, ppi)

fao <- inner_join(fao, idx, by = c("Area", "Year", "Item"))

fao <- dplyr::mutate(fao, ppslc = ppslc/(ppi/100))

rm(idx)


# change (simplify) some commodity names
unique(fao$Item)
# http://www.fao.org/faostat/en/#data/PP définitions and standards > Items for more detail

fao$Item[fao$Item == "Bananas"] <- "Banana" 
fao$Item[fao$Item == "Cocoa, beans"] <- "Cocoa" 
fao$Item[fao$Item == "Coconuts"] <- "Coconut" # the commodity is the coconutS, not the derived oil.
fao$Item[fao$Item == "Coffee, green"] <- "Coffee" # Raw, arabica, robusta, liberica
fao$Item[fao$Item == "Cotton lint"] <- "Cotton"
fao$Item[fao$Item == "Milk, whole fresh cow"] <- "Milk"
fao$Item[fao$Item == "Meat live weight, cattle"] <- "Beef"
fao$Item[fao$Item == "Meat live weight, pig"] <- "Pork"
fao$Item[fao$Item == "Meat live weight, chicken"] <- "Chicken"
fao$Item[fao$Item == "Meat live weight, sheep"] <- "Sheep"
fao$Item[fao$Item == "Oats"] <- "Oat"
fao$Item[fao$Item == "Olives"] <- "Olive" # the commodity is oliveS, incl. those for oil, not the derived oil.
fao$Item[fao$Item == "Oil, palm"] <- "Palm_oil" # the oil, not the fruits (likely crude, although not explicit)
fao$Item[fao$Item == "Oil palm fruit"] <- "FFB"
fao$Item[fao$Item == "Oranges"] <- "Orange" 
fao$Item[fao$Item == "Rice, paddy"] <- "Rice"
fao$Item[fao$Item == "Palm kernels"] <- "Palm_kernel"
fao$Item[fao$Item == "Rubber, natural"] <- "Rubber"
fao$Item[fao$Item == "Soybeans"] <- "Soybean"# the beans, not the oil
fao$Item[fao$Item == "Sugar cane"] <- "Sugarcane"
fao$Item[fao$Item == "Sugar beet"] <- "Sugarbeet"
fao$Item[fao$Item == "Sunflower seed"] <- "Sunflower" # the seeds, not the oil, although mostly values for its oil.
fao$Item[fao$Item == "Groundnuts, with shell"] <- "Groundnut"
fao$Item[fao$Item == "Tobacco, unmanufactured"] <- "Tobacco"

unique(fao$Item)

### Reshape commodities from long to wide 
fao <- dplyr::select(fao, Area, Item, Year, ppslc)
fao <- stats::reshape(data = fao, 
                       direction = "wide", 
                       v.names = "ppslc", # names of variables in the long format that correspond to multiple variables in the wide format
                       timevar = "Item", # the variable in long format that differentiates multiple records from the same group or individual
                       idvar = c("Area", "Year")) 

### Attribute exchange rates to be able to convert international prices 
exr <- read.csv(file = here("input_data", "price_data", "FAOSTAT_exchangerate_6-14-2021.csv"))
names(exr)[names(exr)=="Value"] <- "exr"
fao <- full_join(fao, exr[,c("Area", "Year", "exr")], by = c("Area", "Year"))
fao <- dplyr::arrange(fao, Area, Year)
rm(exr)

fao[fao$exr==0 & !is.na(fao$exr),"exr"] <- NA

### Merge international prices 
names(fao)[names(fao) == "Year"] <- "year"
names(fao)[names(fao) == "Area"] <- "country_name"
names(ip)[names(ip) != "year"] <- paste0("ip_", names(ip)[names(ip) != "year"])
prices <- left_join(fao, ip, by = "year")



### Convert international prices in USD to local currency units
prices <- dplyr::mutate(prices, across(.cols = contains("ip_"),#
                                       .fns = ~.*exr)) # exr is standard local currency units per USD


### Lags

## Lag international prices over whole period 
ip_variables <- names(prices)[grepl(pattern = "ip_", x = names(prices))]

for(voi in ip_variables){
  
  ## short to long lags
  for(lag in c(1:5)){
    prices <- dplyr::arrange(prices, year)
    prices <- DataCombine::slide(prices,
                             Var = voi, 
                             TimeVar = "year",
                             NewVar = paste0(voi,"_lag",lag),
                             slideBy = -lag, 
                             keepInvalid = FALSE)
    prices <- dplyr::arrange(prices, year)
    
  }
  
  for(py in c(2:5)){
    ## Past-year averages (2, 3 and 4 years)  
    # note that we DON'T add voi column (not lagged) in the row mean
    prices$newv <- rowMeans(x = prices[,paste0(voi,"_lag",c(1:py))], na.rm = FALSE)
    prices[is.nan(prices$newv),"newv"] <- NA
    colnames(prices)[colnames(prices)=="newv"] <- paste0(voi,"_",py,"pya")
  }
}

# remove some variables that were only temporarily necessary
# (we want to keep lag1)
vars_torm <- names(prices)[(grepl(pattern = "_lag2", x = names(prices)) |
                          grepl(pattern = "_lag3", x = names(prices)) |
                          grepl(pattern = "_lag4", x = names(prices)) |
                          grepl(pattern = "_lag5", x = names(prices)))]

prices <- prices[,!(names(prices) %in% vars_torm)]


## Lag ppslc prices over period with available data --> first year non missing is 1992 (with lag1 only) 
ppslc_variables <- names(prices)[grepl(pattern = "ppslc.", x = names(prices))]

for(voi in ppslc_variables){
  
  ## short to long lags
  for(lag in c(1:5)){
    prices <- dplyr::arrange(prices, year)
    prices <- DataCombine::slide(prices,
                                 Var = voi, 
                                 TimeVar = "year",
                                 NewVar = paste0(voi,"_lag",lag),
                                 slideBy = -lag, 
                                 keepInvalid = FALSE)
    prices <- dplyr::arrange(prices, year)
    
  }
  
  for(py in c(2:5)){
    ## Past-year averages (2, 3 and 4 years)  
    # note that we DON'T add voi column (not lagged) in the row mean
    prices$newv <- rowMeans(x = prices[,paste0(voi,"_lag",c(1:py))], na.rm = FALSE)
    prices[is.nan(prices$newv),"newv"] <- NA
    colnames(prices)[colnames(prices)=="newv"] <- paste0(voi,"_",py,"pya")
  }
}

# remove some variables that were only temporarily necessary
# (we want to keep lag1)
vars_torm <- names(prices)[(grepl(pattern = "_lag2", x = names(prices)) |
                              grepl(pattern = "_lag3", x = names(prices)) |
                              grepl(pattern = "_lag4", x = names(prices)) |
                              grepl(pattern = "_lag5", x = names(prices)))]

prices <- prices[,!(names(prices) %in% vars_torm)]


### For each lag type and commodity, replace missing values from starting year to corresponding year. 
# FAOSTAT crops
crops <- sub(pattern = ".*ppslc.", x = names(fao)[grepl("ppslc.", names(fao))], replacement = "")

# FAOSTAT producer price time series starts in 1991. Therefore: 
# lag1 variables are avalaible only as of 1992, 2pya as of 1993, 3pya as of 1994, etc. 
lag_types <- c("", "_lag1", "_2pya", "_3pya", "_4pya", "_5pya")
maxyear <- 1991
for(type in lag_types){
  period <- prices$year < maxyear 
  variables <- paste0(crops, type)
  for(voi in variables){
    ip_var <- paste0("ip_", voi)
    ppslc_var <- paste0("ppslc.",voi)
    if(ip_var %in% names(prices)){
      prices[period, ppslc_var] <- prices[period, ip_var] 
    }
  }
  # do the same thing for slightly different commodities:
  prices[period, paste0("ppslc.Rapeseed", type)] <- prices[period, paste0("ip_Rapeseed_oil", type)] 
  prices[period, paste0("ppslc.Sunflower", type)] <- prices[period, paste0("ip_Sunflower_oil", type)] 
  prices[period, paste0("ppslc.Coconut", type)] <- prices[period, paste0("ip_Coconut_oil", type)] 
  prices[period, paste0("ppslc.Olive", type)] <- prices[period, paste0("ip_Olive_oil", type)] 
  prices[period, paste0("ppslc.Sugarbeet", type)] <- prices[period, paste0("ip_Sugar", type)] 
  prices[period, paste0("ppslc.Sugarcane", type)] <- prices[period, paste0("ip_Sugar", type)] 
  
  maxyear <- maxyear + 1
}

# remove ip variables
prices <- dplyr::select(prices, !contains("ip_"))
prices <- dplyr::select(prices, -exr)

names(prices) <- sub(pattern = ".*ppslc.", x = names(prices), replacement = "")

### Make logarithms
prices <- dplyr::mutate(prices, across(.cols = !c("country_name", "year"),
                                       .fns = log,
                                       .names = paste0("ln_", "{.col}")))

# View(ip[,c("Banana", "Banana_lag1", "Banana_2pya", "Banana_3pya", "Banana_4pya")])

saveRDS(prices, here("temp_data", "prepared_producer_prices.Rdata"))



rm(prices, fao, imf, ip, ps, ps2, crops, fpi_2014_16, muv_index_2014_16, period, ppslc_var, ppslc_variables, py, ip_var, ip_variables, lag_types, type, maxyear, ps_commo, needed_imf_col, lag, voi, variables, vars_torm)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


# # Map Pink Sheet and IMF commodities with GAEZ crops
# # Rename layers (will be lost when writing the masked_gaez in the current code, so useless here and we rename later)
# gaez_crops <- list.files(path = here("temp_data", "GAEZ", "AES_index_value", "Rain-fed", "High-input"), 
#                          pattern = "", 
#                          full.names = FALSE)
# gaez_crops <- gsub(pattern = ".tif", replacement = "", x = gaez_crops)
# 
# 
# ps_gaez_map <- data.frame("pink_sheet" = "",
#                           "IMF" = "",
#                           "gaez" = gaez_crops)
# 
# ps_gaez_map[ps_gaez_map$gaez == "Alfalfa", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Banana", "pink_sheet"] <- "Banana, US"
# ps_gaez_map[ps_gaez_map$gaez == "Barley", "pink_sheet"] <- "Barley"
# ps_gaez_map[ps_gaez_map$gaez == "Buckweat", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Cabbage", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Carrot", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Cassava", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Chickpea", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Citrus", "pink_sheet"] <- "Orange"
# ps_gaez_map[ps_gaez_map$gaez == "Cocoa", "pink_sheet"] <- "Cocoa"
# ps_gaez_map[ps_gaez_map$gaez == "Coconut", "pink_sheet"] <- "Coconut oil"
# ps_gaez_map[ps_gaez_map$gaez == "Coffee", "pink_sheet"] <- "Coffee, Arabica" # According to Table A4-3 in GAEZ v3 User guide, coffee is Arabica.
# ps_gaez_map[ps_gaez_map$gaez == "Cotton", "pink_sheet"] <- "Cotton, A Index"
# ps_gaez_map[ps_gaez_map$gaez == "Cowpea", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Dryland_rice", "pink_sheet"] <- "Rice, Thai 5%" # As it seems to be the most consistent series
# ps_gaez_map[ps_gaez_map$gaez == "Drypea", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Flax", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Foxtailmillet", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Gram", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Groundnut", "pink_sheet"] <- "Groundnuts" # take the rawest product available
# ps_gaez_map[ps_gaez_map$gaez == "Jatropha", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Maize", "pink_sheet"] <- "Maize"
# ps_gaez_map[ps_gaez_map$gaez == "Miscanthus", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Oat", "IMF"] <- colnames(imf)[grepl(pattern = "oat", x = imf[1,], ignore.case = TRUE)] # Oats, Generic 1st 'O ' Future, USD/bushel
# ps_gaez_map[ps_gaez_map$gaez == "Oilpalm", "pink_sheet"] <- "Palm oil"
# ps_gaez_map[ps_gaez_map$gaez == "Olive", "IMF"] <- colnames(imf)[grepl(pattern = "olive", x = imf[1,], ignore.case = TRUE)] # Olive Oil, extra virgin less than 1% free fatty acid, ex-tanker price U.K., US$ per metric ton
# ps_gaez_map[ps_gaez_map$gaez == "Onion", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Pearlmillet", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Phaseolusbean", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Pigeonpea", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Rapeseed", "IMF"] <- colnames(imf)[grepl(pattern = "rapeseed", x = imf[1,], ignore.case = TRUE)] # Rapeseed oil, crude, fob Rotterdam, US$ per metric ton
# ps_gaez_map[ps_gaez_map$gaez == "Reedcanarygrass", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Rye", "pink_sheet"] <- "" # cereal grain
# ps_gaez_map[ps_gaez_map$gaez == "Sorghum", "pink_sheet"] <- "Sorghum"
# ps_gaez_map[ps_gaez_map$gaez == "Soybean", "pink_sheet"] <- "Soybeans" # take the rawest product available
# ps_gaez_map[ps_gaez_map$gaez == "Sugarbeet", "pink_sheet"] <- "Sugar, world"
# ps_gaez_map[ps_gaez_map$gaez == "Sugarcane", "pink_sheet"] <- "Sugar, world"
# ps_gaez_map[ps_gaez_map$gaez == "Sunflower", "IMF"] <- colnames(imf)[grepl(pattern = "sunflower", x = imf[1,], ignore.case = TRUE)] # Sunflower oil, Sunflower Oil, US export price from Gulf of Mexico, US$ per metric ton
# ps_gaez_map[ps_gaez_map$gaez == "Sweetpotato", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Switchgrass", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Tea", "pink_sheet"] <- "Tea, avg 3 auctions"
# ps_gaez_map[ps_gaez_map$gaez == "Tobacco", "pink_sheet"] <- "Tobacco, US import u.v."
# ps_gaez_map[ps_gaez_map$gaez == "Tomato", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Wetland_rice", "pink_sheet"] <- "Rice, Thai 5%" # same price as dryland rice, as output is the same...
# ps_gaez_map[ps_gaez_map$gaez == "Wheat", "pink_sheet"] <- "Wheat, US HRW" # hard red winter wheat is more common
# ps_gaez_map[ps_gaez_map$gaez == "Whitepotato", "pink_sheet"] <- ""
# ps_gaez_map[ps_gaez_map$gaez == "Yam", "pink_sheet"] <- ""
# 
# 
# 
# 
# 
# data <- readRDS(here("temp_data", "merged_datasets", "tropical_aoi", "glass_gaez_long.Rdata"))
# 
# ex <- data[1, 1:47 ]
# std <- ex/sum(ex)
# 
# names(data)
