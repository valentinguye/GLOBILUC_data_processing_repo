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

years_oi <- seq(1983, 2020, 1) 


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


#### PINK SHEET DATA (WORLD BANK) ####
# "Annual Prices" link at: https://www.worldbank.org/en/research/commodity-markets (in April 2021)
ps <- read_xlsx(here("input_data", "price_data", "CMO-Historical-Data-Annual.xlsx"), 
                sheet = "Annual Prices (Real)", 
                range = "A7:BS69", 
                col_names =TRUE)

names(ps)[1] <- "year"

ps <- as.data.frame(ps)

# drop first row (price unit) for manipulation convenience. 
ps <- ps[-1,]

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

ps2 <- ps[, c("year", ps_commo)] 

# Select years of interest 
ps2 <- dplyr::filter(ps2, ps2$year %in% years_oi)


# Convert to numeric
ps2 <- summarise(ps2, across(.fns = as.numeric))

# Change (simplify) some names
names(ps2)[names(ps2) == "Banana, US"] <- "Banana"
names(ps2)[names(ps2) == "Coconut oil"] <- "Coconut_oil"
names(ps2)[names(ps2) == "Coffee, Arabica"] <- "Coffee"
names(ps2)[names(ps2) == "Cotton, A Index"] <- "Cotton"
names(ps2)[names(ps2) == "Crude oil, average"] <- "Crude_oil"
names(ps2)[names(ps2) == "Meat, chicken"] <- "Chicken"
names(ps2)[names(ps2) == "Meat, sheep"] <- "Sheep"
names(ps2)[names(ps2) == "Rice, Thai 5%"] <- "Rice"
names(ps2)[names(ps2) == "Palm oil"] <- "Palm_oil"
names(ps2)[names(ps2) == "Rubber, SGP/MYS"] <- "Rubber"
names(ps2)[names(ps2) == "Soybean oil"] <- "Soybean_oil"
names(ps2)[names(ps2) == "Soybean meal"] <- "Soybean_meal"
names(ps2)[names(ps2) == "Sugar, world"] <- "Sugar"
names(ps2)[names(ps2) == "Tea, avg 3 auctions"] <- "Tea"
names(ps2)[names(ps2) == "Tobacco, US import u.v."] <- "Tobacco"
names(ps2)[names(ps2) == "Wheat, US HRW"] <- "Wheat"

head(ps2)


#### IMF DATA #### 
# "Excel Database" link at: webpage: https://www.imf.org/en/Research/commodity-prices (in April 2021)
imf <- read_xls(here("input_data", "price_data", "external-dataAPR.xls"), 
                col_names = TRUE)
length(unique(colnames(imf))) == ncol(imf)
names(imf)[1] <- "year"
# Select columns we need 
colnames(imf)
needed_imf_col <- c("year", 
                    colnames(imf)[grepl(pattern = "oat", x = imf[1,], ignore.case = TRUE)], 
                    colnames(imf)[grepl(pattern = "olive", x = imf[1,], ignore.case = TRUE)], 
                    colnames(imf)[grepl(pattern = "rapeseed", x = imf[1,], ignore.case = TRUE)], 
                    colnames(imf)[grepl(pattern = "sunflower", x = imf[1,], ignore.case = TRUE)],
                    colnames(imf)[grepl(pattern = "pork", x = imf[1,], ignore.case = TRUE)])
    
imf <- imf[, needed_imf_col]

# drop first, unnecessary rows
imf <- imf[-c(1,2,3),]

imf <- as.data.frame(imf)

# Rename
names(imf) <- c("year", "Oat", "Olive_oil", "Rapeseed_oil", "Sunflower_oil", "Pork")

# Convert in 2010 real value (instead of 2016), using MUV index from Pink Sheet
muv_index_2016 <- ps[ps$year == 2016, "MUV Index"] %>% as.numeric()/100

imf <- mutate(imf, Oat = as.numeric(Oat)/muv_index_2016, 
              Olive_oil = as.numeric(Olive_oil)/muv_index_2016, 
              Rapeseed_oil = as.numeric(Rapeseed_oil)/muv_index_2016, 
              Sunflower_oil = as.numeric(Sunflower_oil)/muv_index_2016,
              Pork = as.numeric(Pork)/muv_index_2016)
class(imf[,2])

# Average monthly prices to annual 
imf$year <- gsub("M.*", "", imf$year) %>% as.numeric()

annual_imf <- ddply(imf, "year", summarise, 
                    Oat = mean(Oat), 
                    Olive_oil = mean(Olive_oil), 
                    Rapeseed_oil = mean(Rapeseed_oil), 
                    Sunflower_oil = mean(Sunflower_oil), 
                    Pork = mean(Pork))


# Restrict to years of interest
annual_imf <- dplyr::filter(annual_imf, annual_imf$year %in% years_oi)



#### Consolidate price annual time series data frame #### 

prices <- merge(x = ps2, y = annual_imf, by = "year") 
head(prices)

#### PREPARE ADDITIONAL PRICE VARIABLES #### 

### Make logarithms
logs <- mutate(prices[,names(prices)[names(prices) != "year"]], 
               across(.fns = log))

names(logs) <- paste0("ln_", names(logs))

prices <- cbind(prices, logs)
head(prices)

### Lags 
variables <- names(prices)[grepl(pattern = "ln_", x = names(prices))]

for(voi in variables){
  for(lag in c(1:2)){
    prices <- dplyr::arrange(prices, year)
    prices <- DataCombine::slide(prices,
                                  Var = voi, 
                                  TimeVar = "year",
                                  NewVar = paste0(voi,"_lag", lag),
                                  slideBy = -lag, 
                                  keepInvalid = FALSE)
    prices <- dplyr::arrange(prices, year)
  }
}

# bananas <- prices[, c("year", "Banana", "ln_Banana", "ln_Banana_lag1", "ln_Banana_lag2")] 


saveRDS(prices, here("temp_data", "prepared_prices.Rdata"))



rm(annual_imf, imf, logs, prices, ps, ps2, muv_index_2016, ps_commo, needed_imf_col, years_oi)



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
ps_gaez_map[ps_gaez_map$gaez == "Alfalfa", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Banana", "pink_sheet"] <- "Banana, US"
ps_gaez_map[ps_gaez_map$gaez == "Barley", "pink_sheet"] <- "Barley"
ps_gaez_map[ps_gaez_map$gaez == "Buckweat", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Cabbage", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Carrot", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Cassava", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Chickpea", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Citrus", "pink_sheet"] <- "Orange"
ps_gaez_map[ps_gaez_map$gaez == "Cocoa", "pink_sheet"] <- "Cocoa"
ps_gaez_map[ps_gaez_map$gaez == "Coconut", "pink_sheet"] <- "Coconut oil"
ps_gaez_map[ps_gaez_map$gaez == "Coffee", "pink_sheet"] <- "Coffee, Arabica" # According to Table A4-3 in GAEZ v3 User guide, coffee is Arabica.
ps_gaez_map[ps_gaez_map$gaez == "Cotton", "pink_sheet"] <- "Cotton, A Index"
ps_gaez_map[ps_gaez_map$gaez == "Cowpea", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Dryland_rice", "pink_sheet"] <- "Rice, Thai 5%" # As it seems to be the most consistent series
ps_gaez_map[ps_gaez_map$gaez == "Drypea", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Flax", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Foxtailmillet", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Gram", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Groundnut", "pink_sheet"] <- "Groundnuts" # take the rawest product available
ps_gaez_map[ps_gaez_map$gaez == "Jatropha", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Maize", "pink_sheet"] <- "Maize"
ps_gaez_map[ps_gaez_map$gaez == "Miscanthus", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Oat", "IMF"] <- colnames(imf)[grepl(pattern = "oat", x = imf[1,], ignore.case = TRUE)] # Oats, Generic 1st 'O ' Future, USD/bushel
ps_gaez_map[ps_gaez_map$gaez == "Oilpalm", "pink_sheet"] <- "Palm oil"
ps_gaez_map[ps_gaez_map$gaez == "Olive", "IMF"] <- colnames(imf)[grepl(pattern = "olive", x = imf[1,], ignore.case = TRUE)] # Olive Oil, extra virgin less than 1% free fatty acid, ex-tanker price U.K., US$ per metric ton
ps_gaez_map[ps_gaez_map$gaez == "Onion", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Pearlmillet", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Phaseolusbean", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Pigeonpea", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Rapeseed", "IMF"] <- colnames(imf)[grepl(pattern = "rapeseed", x = imf[1,], ignore.case = TRUE)] # Rapeseed oil, crude, fob Rotterdam, US$ per metric ton
ps_gaez_map[ps_gaez_map$gaez == "Reedcanarygrass", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Rye", "pink_sheet"] <- "" # cereal grain
ps_gaez_map[ps_gaez_map$gaez == "Sorghum", "pink_sheet"] <- "Sorghum"
ps_gaez_map[ps_gaez_map$gaez == "Soybean", "pink_sheet"] <- "Soybeans" # take the rawest product available
ps_gaez_map[ps_gaez_map$gaez == "Sugarbeet", "pink_sheet"] <- "Sugar, world"
ps_gaez_map[ps_gaez_map$gaez == "Sugarcane", "pink_sheet"] <- "Sugar, world"
ps_gaez_map[ps_gaez_map$gaez == "Sunflower", "IMF"] <- colnames(imf)[grepl(pattern = "sunflower", x = imf[1,], ignore.case = TRUE)] # Sunflower oil, Sunflower Oil, US export price from Gulf of Mexico, US$ per metric ton
ps_gaez_map[ps_gaez_map$gaez == "Sweetpotato", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Switchgrass", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Tea", "pink_sheet"] <- "Tea, avg 3 auctions"
ps_gaez_map[ps_gaez_map$gaez == "Tobacco", "pink_sheet"] <- "Tobacco, US import u.v."
ps_gaez_map[ps_gaez_map$gaez == "Tomato", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Wetland_rice", "pink_sheet"] <- "Rice, Thai 5%" # same price as dryland rice, as output is the same...
ps_gaez_map[ps_gaez_map$gaez == "Wheat", "pink_sheet"] <- "Wheat, US HRW" # hard red winter wheat is more common
ps_gaez_map[ps_gaez_map$gaez == "Whitepotato", "pink_sheet"] <- ""
ps_gaez_map[ps_gaez_map$gaez == "Yam", "pink_sheet"] <- ""
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
