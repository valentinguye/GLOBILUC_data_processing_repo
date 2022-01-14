
### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
neededPackages = c("Matrix",
                   "plyr", "dplyr", "here",#"tibble", "data.table",
                   "foreign", "readxl",
                   "raster", "rgdal",  "sp", "sf", # "spdep",
                   "DataCombine",
                   "knitr", "kableExtra",
                   "fixest", "boot",#,"msm", "car",  "sandwich", "lmtest",  "multcomp",
                   "ggplot2", "dotwhisker", #"tmap",# "leaflet", "htmltools"
                   "foreach", "parallel"
)
# "pglm", "multiwayvcov", "clusterSEs", "alpaca", "clubSandwich",

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
              "Groundnut oil", 
              "Maize",
              "Palm oil", 
              "Palm kernel oil",
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
names(ps2)[names(ps2) == "Coffee, Arabica"] <- "Coffee" # for the sake of simplicity, because this is the one we will use with GAEZ. 
names(ps2)[names(ps2) == "Cotton, A Index"] <- "Cotton"
names(ps2)[names(ps2) == "Crude oil, average"] <- "Crude_oil"
names(ps2)[names(ps2) == "Groundnuts"] <- "Groundnuts"
names(ps2)[names(ps2) == "Groundnut oil"] <- "Groundnut_oil"
names(ps2)[names(ps2) == "Meat, chicken"] <- "Chicken"
names(ps2)[names(ps2) == "Meat, sheep"] <- "Sheep"
names(ps2)[names(ps2) == "Rice, Thai 5%"] <- "Rice"
names(ps2)[names(ps2) == "Palm oil"] <- "Palm_oil"
names(ps2)[names(ps2) == "Palm kernel oil"] <- "Palm_kernel_oil"
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
imf <- as.data.frame(imf)
# Select columns we need 
imf <- imf[,(colnames(imf) %in% c("year", "POATS", "POLVOIL", "PFOOD", "PBEVE", "PMILK",   "PROIL", "PSUNO", "PPORK", "PFERT"))]

# drop first, unnecessary rows
imf <- imf[-c(1,2,3),]

head(imf)

# Rename /!\ ORDER MATTERS HERE it must match the observed order, not the one in column selection above
names(imf)
names(imf) <- c("year", "FPI", "BPI", "Fertilizer", "Rapeseed_oil", "Olive_oil",  "Pork", "Sunflower_oil", "Oat", "Milk")#"Coffee",

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
                    Oat = Oat %>% as.numeric() %>% mean(na.rm = TRUE), 
                    Olive_oil = Olive_oil %>% as.numeric() %>% mean(na.rm = TRUE), 
                    Milk = Milk %>% as.numeric() %>% mean(na.rm = TRUE), 
                    Rapeseed_oil = Rapeseed_oil %>% as.numeric() %>% mean(na.rm = TRUE), 
                    Sunflower_oil = Sunflower_oil %>% as.numeric() %>% mean(na.rm = TRUE), 
                    # Coffee = mean(Coffee), 
                    Pork = Pork %>% as.numeric() %>% mean(na.rm = TRUE), 
                    Fertilizer = Fertilizer %>% as.numeric() %>% mean(na.rm = TRUE), 
                    FPI = FPI %>% as.numeric() %>% mean(na.rm = TRUE), 
                    BPI = BPI %>% as.numeric() %>% mean(na.rm = TRUE))

# From 1992 on, use the Food price index (FPI) from the IMF, in base 2016=100, to deflate the nominal time series.
# Except for fertilizers, which is already an index base 100 in 2016
imf[imf$year >=1992,] <- mutate(imf[imf$year >=1992,], Oat = as.numeric(Oat)/(FPI/100), 
                                                       Olive_oil = as.numeric(Olive_oil)/(FPI/100), 
                                                       Rapeseed_oil = as.numeric(Rapeseed_oil)/(FPI/100), 
                                                       Sunflower_oil = as.numeric(Sunflower_oil)/(FPI/100),
                                                       # Coffee = as.numeric(Coffee)/muv_index_2016,
                                                       Pork = as.numeric(Pork)/(FPI/100))

# Then change the base from 2016 to 2014-2016, to match FAOSTAT (including fertilizers)
fpi_2014_16 <- mean(c(as.numeric(imf[imf$year >= 2014 & imf$year <= 2016, "FPI"])))
# do it for the Beverage Price Index too
bpi_2014_16 <- mean(c(as.numeric(imf[imf$year >= 2014 & imf$year <= 2016, "BPI"])))

imf <- dplyr::mutate(imf, across(.cols = (!contains("year") & !contains("BPI")), 
                                 .fns = ~.*fpi_2014_16/100)) 

# For years prior 1992, we deflate with the MUV index base 2014-2016 from the Pink Sheet (do not include fertilizer index)
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


#### ICO DATA for COFFEE ####
# data are downloaded here
# https://www.ico.org/new_historical.asp

ico <- read_excel(here("input_data", "price_data", "ICO_indicator_prices.xlsx")) %>% as.data.frame()
# two first columns are time and ico composite indicator 
ico <- ico[,1:2]
colnames(ico) <- c("year", "Coffee")
ico <- ico[(grepl(pattern = "19", ico$year) | grepl(pattern = "20", ico$year)), ]
ico$year <- as.numeric(ico$year)
ico$Coffee <- as.numeric(ico$Coffee)
ico$Coffee <- round(ico$Coffee, 2)
ico$Coffee <- ico$Coffee*0.01/0.000453592 

# Deflate using IMF beverage price index
# keep only prices 
ico <- left_join(ico, imf[,c("year", "BPI")], by = "year")
ico <- mutate(ico, Coffee = Coffee*(bpi_2014_16/100) / (BPI/100) )


# compare with arabicas from PS
summary(ps2[ps2$year >= 2000, "Coffee"])
summary(ico[ico$year >= 2000, "Coffee"])

# the distribution is not the same. 2 things explain the difference: first and mostly, it's that ico ts is deflated with BPI, not pinksheet index 
# second, ico ts is composite of arabica and robusta
# for both reasons, we prefer to use ico time series
ps2 <- dplyr::select(ps2, -Coffee)


#### WEIGHTS ####

### VEGETABLE OILS WEIGHTS AND CEREALS-FEEDS WEIGHTS
# here, use faostat data
faocer <- read.csv(file = here("input_data", "FAOSTAT_data_12-10-2021.csv"))
unique(faocer$Item)
template <- faocer[!duplicated(faocer$Year), c("Year", "Unit")]
cereals <- c("Barley", "Maize", "Maize, green", "Oats", "Rice, paddy", "Sorghum", "Soybeans", "Wheat") # they are all expressed in tonnes
for(cer in cereals){
  cer_df <- faocer[faocer$Item == cer, c("Year", "Value")]
  names(cer_df)[names(cer_df)=="Value"] <- cer
  template <- full_join(template, cer_df, by = "Year")
}
rm(cer_df)
# Group maize 
faocer <- mutate(template, Maize := Maize + !!as.symbol("Maize, green"))
cereals <- cereals[cereals != "Maize, green"]
faocer <- faocer[,c("Year", cereals)]

# Weight Soybeans to not double factor in soybean_oil used for food or industrial use (17% of the total according to OurWorldinData) 
faocer <- mutate(faocer, Soybeans = Soybeans*0.83)

# Make shares of total supply for each crop
faocer <- mutate(faocer, cereal_supply = rowSums(across(.cols = !contains("Year"))))
for(cer in cereals){
  faocer <- mutate(faocer, !!as.symbol(cer) := !!as.symbol(cer)/cereal_supply)
}

# change names 
names(faocer) <- c("year", "share_Barley", "share_Maize", "share_Oat", "share_Rice", "share_Sorghum", "share_Soybean_meal", "share_Wheat", "cereal_supply")

### SOY WEIGHTS
# here, use USDA PSD data
# here we make a price index for soy. 
psd <- read.csv(here("input_data", "USDA", "psd_alldata_csv", "psd_alldata.csv"))

meal <- filter(psd, Country_Name %in% c("Brazil","China", "Argentina", "India", "United States"), #"Brazil", "China", "Argentina", "India", 
               Commodity_Description %in% c("Meal, Soybean"), 
               Attribute_Description == "Production") %>% ddply(c("Market_Year"), summarise, 
                                                                Soybean_meal = sum(Value, na.rm = TRUE))

oil <- filter(psd, Country_Name %in% c("Brazil","China", "Argentina", "India", "United States"), 
              Commodity_Description %in% c("Oil, Soybean"), 
              Attribute_Description == "Production") %>% ddply("Market_Year", summarise, 
                                                               Soybean_oil = sum(Value, na.rm = TRUE))
# this is the total soybean production. 
# total <- filter(psd, Country_Name %in% c("Brazil","China", "Argentina", "India", "United States"),
#                Commodity_Description %in% c("Oilseed, Soybean"),
#                Attribute_Description == "Production") %>% ddply("Market_Year", summarise,
#                                                                   Soybeans = sum(Value, na.rm = TRUE))

soy <- full_join(oil, meal, "Market_Year")
# soy <- full_join(soy, total, "Market_Year")

# give them a different name (soyshare) because share_Soybean_meal is already taken
soy <- mutate(soy, soyshare_Soybean_meal = Soybean_meal/(Soybean_meal + Soybean_oil),
              soyshare_Soybean_oil = Soybean_oil/(Soybean_meal + Soybean_oil))

names(soy)[names(soy)=="Market_Year"] <- "year"

#### Consolidate international price annual time series data frame #### 
# keep all years from each sources
ip <- full_join(x = ps2, y = imf, by = "year") 
ip <- full_join(x = ip, y = ico[,c("year", "Coffee")],  by = "year") # necessary to not merge BPI column again
ip <- ip %>% dplyr::select(-MUV_index, -FPI, -BPI, -MUV_index_10, - MUV_index_2014_16)

### Make calorie prices

### CALORIE CONTENT OF CROPS *before log*

# values are in calories per kg, from https://www.fao.org/3/x9892e/X9892e05.htm
# or kcalories per ton
cal_content <- c(
  "Banana" = 600,
  "Barley" = 3320,
  "Beef" = 1500, 
  "Orange" = 340, 
  "Cocoa" = 4140,
  "Coconut_oil" = 8840,
  "Coffee" = 470,
  "Cotton" = 2530, # we have prices not on cotton oil but on cotton raw from Cotlook 'A' Index
  "Groundnut_oil" = 8840,
  "Maize" = 3560,
  "Oat" = 3850,
  "Olive_oil" = 8840,
  "Palm_oil" = 8840,
  "Palm_kernel_oil" = 8840,
  "Rapeseed_oil" = 8840,
  "Rice" = 3600, # it's not rice paddy, but 5% broken, so we take the value for rice milled from calorie content source
  "Sorghum" = 3430,
  "Soybean_oil" = 8840,
  "Soybean_meal" = 2610,
  "Sugar" = 3730, # it's raw in Pink Sheet
  "Sunflower_oil" = 8840,
  "Tea" = 400,
  "Wheat" = 3340)

for(commo in names(cal_content)){
  ip <- mutate(ip, !!as.symbol(paste0("kca_",commo)) := !!as.symbol(commo) / cal_content[commo])  
}

### Make logarithms
ip <- dplyr::mutate(ip, across(.cols = !c("year"),
                         .fns = log,
                         .names = paste0("ln_", "{.col}")))

### Make weighted prices
ip <- left_join(x = ip, y = faocer, by = "year")

# Weight only log and kca prices

## Cereals_feeds
for(wc in c("Barley", "Maize", "Oat", "Rice", "Sorghum", "Soybean_meal", "Wheat")){
  # ip <- mutate(ip, 
  #              !!as.symbol(paste0("w_ln_",wc)) := !!as.symbol(paste0("share_",wc)) * (!!as.symbol(paste0("ln_",wc))) )
  # intervert kca_ and ln_ in variable names for convenience in make_main_reg function
  ip <- mutate(ip, 
               !!as.symbol(paste0("w_kca_ln_",wc)) := !!as.symbol(paste0("share_",wc)) * (!!as.symbol(paste0("ln_kca_",wc))) )
}

## Vegetable_oils
# This is a share of world supply, as per the USDA WorldSupplyUseOilseedandProducts, Table 42–World vegetable oils supply and distribution, 2013/14–2020/21
# Quantities are averaged from year 2013/14 to 2018/19 
oil_sha <- c(0.019101432,0.026418347,0.355059925,0.041790248,0.030875014,0.150684719,0.283082617,0.092987698)
# weights with olive oil
# oil_sha <- c(0.018795739,0.025995556,0.01600368,0.34937766,0.04112145,0.0303809,0.148273209,0.278552254,0.091499553)
names(oil_sha) <- c("Coconut_oil","Cotton","Palm_oil","Palm_kernel_oil","Groundnut_oil","Rapeseed_oil","Soybean_oil","Sunflower_oil")
oil_sha <- oil_sha[base::order(names(oil_sha))]

for(wc in c(names(oil_sha))){
  ip <- mutate(ip, 
               !!as.symbol(paste0("w_kca_ln_",wc)) := oil_sha[wc] * (!!as.symbol(paste0("ln_kca_",wc))) )
}

# do not average weighted prices now, because this will be done more flexibly in the analysis directly. 

## Soy

ip <- left_join(ip, soy[,c("year", "soyshare_Soybean_meal", "soyshare_Soybean_oil")], "year")

ip <- mutate(ip, Soy_index = soyshare_Soybean_meal * Soybean_meal + soyshare_Soybean_oil * Soybean_oil)
ip <- mutate(ip, ln_Soy_index = log(Soy_index))


## IF WE ARE TO MERGE MACRO FROM USDA
# ip <- full_join(x = ip, y = psd, by = "year")

# macro_vars <- ip
# for(country in c("United States", "EU-25", "China")){
#   for(attribute_des in c("Area Harvested", "Total Supply", "Domestic Consumption")){
#     macro_vars <- full_join(macro_vars, psd_list[[paste0(country, " ", attribute_des)]], by = "year")
#   }
# }



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
  
  # past year averages
  for(py in c(2:5)){
    ## Past-year averages (2, 3 and 4 years)  
    # note that we DON'T add voi column (not lagged) in the row mean
    inter_prices$newv <- rowMeans(x = inter_prices[,paste0(voi,"_lag",c(1:py))], na.rm = FALSE)
    inter_prices[is.nan(inter_prices$newv),"newv"] <- NA
    colnames(inter_prices)[colnames(inter_prices)=="newv"] <- paste0(voi,"_",py,"pya")
  }
  
  # leads
  for(lead in 1){
    inter_prices <- dplyr::arrange(inter_prices, year)
    inter_prices <- DataCombine::slide(inter_prices,
                                       Var = voi, 
                                       TimeVar = "year",
                                       NewVar = paste0(voi,"_lead",lead),
                                       slideBy = lead, 
                                       keepInvalid = FALSE)
    inter_prices <- dplyr::arrange(inter_prices, year)
  }
}

# remove some variables that were only temporarily necessary
# (we want to keep lag1)
vars_torm <- names(inter_prices)[(#grepl(pattern = "_lag2", x = names(inter_prices)) |
                              grepl(pattern = "_lag3", x = names(inter_prices)) |
                              grepl(pattern = "_lag4", x = names(inter_prices)) |
                              grepl(pattern = "_lag5", x = names(inter_prices)))]

inter_prices <- inter_prices[,!(names(inter_prices) %in% vars_torm)]

saveRDS(inter_prices, here("temp_data", "prepared_international_prices.Rdata"))

rm(inter_prices)





#### WEIGHTS ####
# this is the code if we want to compute a data frame of weights on its own. Currently not needed, as weighting is done in this script, above. 
#fg <- read_xlsx("input_data", "USDA", "FeedGrainsYearbookTables-AllYears.xlsx", sheet = "FGYearbookTable02-Full") 
# 
# faocer <- read.csv(file = here("input_data", "FAOSTAT_data_12-10-2021.csv"))
# unique(faocer$Item)
# template <- faocer[!duplicated(faocer$Year), c("Year", "Unit")]
# cereals <- c("Barley", "Maize", "Maize, green", "Oats", "Rice, paddy", "Sorghum", "Soybeans", "Wheat") # they are all expressed in tonnes
# for(cer in cereals){
#   cer_df <- faocer[faocer$Item == cer, c("Year", "Value")]
#   names(cer_df)[names(cer_df)=="Value"] <- cer
#   template <- full_join(template, cer_df, by = "Year")
# }
# rm(cer_df)
# # Group maize 
# faocer <- mutate(template, Maize := Maize + !!as.symbol("Maize, green"))
# cereals <- cereals[cereals != "Maize, green"]
# faocer <- faocer[,c("Year", cereals)]
# 
# 
# # Make shares of total supply for each crop
# faocer <- mutate(faocer, cereal_supply = rowSums(across(.cols = !contains("Year"))))
# for(cer in cereals){
#   faocer <- mutate(faocer, !!as.symbol(cer) := !!as.symbol(cer)/cereal_supply)
# }
# 
# # change names 
# names(faocer) <- c("year", "share_Barley", "share_Maize", "share_Oat", "share_Rice", "share_Sorghum", "share_Soybeans", "share_Wheat", "cereal_supply")
# 
# # lag 
# ip_variables <- names(faocer)[names(faocer)!="year"]
# for(voi in ip_variables){
#   ## short to long lags
#   for(lag in c(1:5)){
#     faocer <- dplyr::arrange(faocer, year)
#     faocer <- DataCombine::slide(faocer,
#                                  Var = voi, 
#                                  TimeVar = "year",
#                                  NewVar = paste0(voi,"_lag",lag),
#                                  slideBy = -lag, 
#                                  keepInvalid = FALSE)
#     faocer <- dplyr::arrange(faocer, year)
#     
#   }
#   for(py in c(2:5)){
#     ## Past-year averages (2, 3 and 4 years)  
#     # note that we DON'T add voi column (not lagged) in the row mean
#     faocer$newv <- rowMeans(x = faocer[,paste0(voi,"_lag",c(1:py))], na.rm = FALSE)
#     faocer[is.nan(faocer$newv),"newv"] <- NA
#     colnames(faocer)[colnames(faocer)=="newv"] <- paste0(voi,"_",py,"pya")
#   }
# }
# 
# # remove some variables that were only temporarily necessary
# # (we want to keep lag1)
# vars_torm <- names(faocer)[(#grepl(pattern = "_lag2", x = names(faocer)) |
#   grepl(pattern = "_lag3", x = names(faocer)) |
#     grepl(pattern = "_lag4", x = names(faocer)) |
#     grepl(pattern = "_lag5", x = names(faocer)))]
# 
# faocer <- faocer[,!(names(faocer) %in% vars_torm)]
# 
# 
# saveRDS(faocer, here("temp_data", "fao_cereals_shares_global_supply.Rdata"))
# 
# rm(template, cereals)




#### USDA MACRO ####

psd <- read.csv(here("input_data", "USDA", "psd_alldata_csv", "psd_alldata.csv"))

## What countries are there
sort(unique(psd$Country_Name))

uspsd <- psd[psd$Country_Name == "United States",]

## What agregate to use? 
unique(uspsd$Attribute_Description)
# Total Use = domestic consumption + exports + ending stocks
uspsd_TU <- uspsd[uspsd$Attribute_Description == "Total Use",]
uspsd_TU$Calendar_Year%>%unique()%>% sort()# - misses years 2001-2007
uspsd_TU$Market_Year%>%unique()%>% sort()
unique(uspsd_TU$Commodity_Description) # only very few commodities
# Total Supply = beginning stocks + domestic production + imports
uspsd_TS <- uspsd[uspsd$Attribute_Description == "Total Supply",]
uspsd_TS$Calendar_Year%>%unique()%>% sort()
uspsd_TS$Market_Year%>%unique()%>% sort()
unique(uspsd_TS$Commodity_Description)
# Domestic Consumption = all possible uses of the commodity: food, feed, seed, waste, and industrial processing
uspsd_DC <- uspsd[uspsd$Attribute_Description == "Domestic Consumption",]
uspsd_DC$Calendar_Year%>%unique()%>% sort()# - misses years 2001-2005
uspsd_DC$Market_Year%>%unique()%>% sort()
unique(uspsd_DC$Commodity_Description)

# Area Harvested
uspsd_AH <- uspsd[uspsd$Attribute_Description == "Area Harvested",]
uspsd_AH$Calendar_Year%>%unique()%>% sort() # only since 2006
uspsd_AH$Market_Year%>%unique()%>% sort()
unique(uspsd_AH$Commodity_Description) # main crops (by crop, not commodity, for instance only one category for soy)

# they all have all years of information, conditional on using marketing year and not calendar year. 

## What commodities?
unique(uspsd$Commodity_Description)




# select those we want
focal_commodities <- c("Animal Numbers, Cattle", "Corn",
                       "Meal, Rapeseed", "Meal, Soybean", "Meal, Sunflowerseed", # these categories are only in TS and DC, not in AH
                       "Oil, Rapeseed", "Oil, Soybean", "Oil, Sunflowerseed", # these categories are only in TS and DC, not in AH
                       "Oilseed, Rapeseed", "Oilseed, Soybean", "Oilseed, Sunflowerseed", 
                       "Rice, Milled", "Sugar, Centrifugal", "Wheat")
## format data 
psd_list <- list()
elm <- 1
for(country in c("United States", "European Union", "China")){
  for(attribute_des in c("Area Harvested", "Total Supply", "Domestic Consumption")){
    
    # select observations and variables we want
    tmpd<- psd %>% dplyr::filter(Country_Name == country &
                                   Attribute_Description == attribute_des & 
                                   Commodity_Description %in% focal_commodities & 
                                   Market_Year >= 1980 & 
                                   Market_Year <= 2021)
    
    print(unique(tmpd$Unit_Description))
    print(tmpd[tmpd$Unit_Description=="(1000 HEAD)", "Commodity_Description"]%>% unique())
    
    tmpd <- dplyr::select(tmpd, Commodity_Description, Market_Year, Value)
    
    # rename some commodities and variables
    tmpd[tmpd$Commodity_Description == "Animal Numbers, Cattle","Commodity_Description"] <- "Cattle"
    tmpd[tmpd$Commodity_Description == "Corn","Commodity_Description"] <- "Maize" 
    tmpd[tmpd$Commodity_Description == "Meal, Rapeseed","Commodity_Description"] <- "Rapeseed_meal"
    tmpd[tmpd$Commodity_Description == "Meal, Soybean","Commodity_Description"] <- "Soybean_meal"
    tmpd[tmpd$Commodity_Description == "Meal, Sunflowerseed","Commodity_Description"] <- "Sunflowerseed_meal"
    tmpd[tmpd$Commodity_Description == "Oil, Rapeseed","Commodity_Description"] <- "Rapeseed_oil"
    tmpd[tmpd$Commodity_Description == "Oil, Soybean","Commodity_Description"] <- "Soybean_oil"
    tmpd[tmpd$Commodity_Description == "Oil, Sunflowerseed","Commodity_Description"] <- "Sunflowerseed_oil"
    tmpd[tmpd$Commodity_Description == "Oilseed, Rapeseed","Commodity_Description"] <- "Rapeseed_oilseed"
    tmpd[tmpd$Commodity_Description == "Oilseed, Soybean","Commodity_Description"] <- "Soybean_oilseed"
    tmpd[tmpd$Commodity_Description == "Oilseed, Sunflowerseed","Commodity_Description"] <- "Sunflowerseed_oilseed"
    tmpd[tmpd$Commodity_Description == "Rice, Milled","Commodity_Description"] <- "Rice"
    tmpd[tmpd$Commodity_Description == "Sugar, Centrifugal","Commodity_Description"] <- "Sugar"
    
    names(tmpd)[names(tmpd)=="Market_Year"] <- "year"
    names(tmpd)[names(tmpd)=="Value"] <-  gsub(" ", "_", attribute_des)
    
    
    tmpd <- reshape(tmpd, direction = "wide", 
                    v.names = gsub(" ", "_", attribute_des),
                    timevar = "Commodity_Description",
                    idvar = "year")
    
    # add the country in the variable name
    names(tmpd)[names(tmpd)!="year"] <- paste0(gsub(" ", "", country), ".", names(tmpd)[names(tmpd)!="year"]) 
    
    
    # and lines if there are not as many as in other time series
    min_year <- min(tmpd$year, na.rm = TRUE) 
    if(min_year>1980){
      tmpd <- rbind(matrix(NA, nrow = min_year - 1980, ncol = ncol(tmpd), dimnames = list(NULL, colnames(tmpd)) ), tmpd)#rep(NA, ncol(tmpd))
      tmpd[is.na(tmpd$year),"year"] <- (1980 : (min_year - 1))
    }
    
    # transform transformed commodities back to quantities of a single crop.
    # for(crop in c("Rapeseed", "Soybean", "Sunflower")){
    #   commos <- names(tmpd)[grepl(crop, names(tmpd))]
    #   log_commos <- c()
    #   for(co in commos){
    #     log_co <- paste0("log_", co)
    #     log_commos <- c(log_commos, log_co)
    #     tmpd <- mutate(tmpd, !!as.symbol(log_co) := log(!!as.symbol(co)))
    #   }
    #   varname <- paste0(crop, "_commos")
    #   tmpd <- mutate(tmpd, 
    #                  !!as.symbol(varname) := rowMeans(across(.cols = (any_of(log_commos)))))
    # }
    
    psd_list[[paste0(country, " ", attribute_des)]]  <- tmpd   
    
  }
}

# UNITS
# so from the prints: in AH, unit is 1000ha for all goods, 
# for TS is 1000 head only for Cattle, and 1000 MT for all other goods, 
# and for DC it's 1000 MT for all goods

psd <- bind_cols(psd_list)
names(psd)[1] <- "year"
psd <- psd[, !grepl("year.", names(psd))]


# transform transformed commodities back to quantities of a single crop.

# and quantities of interest (ratios?)






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


coffee_p <- prices[,names(prices)%in%c("year", "country_name", "ppslc.Coffee", "ip_Coffee", "ip_Coffee, Robusta")]
coffee_p <- dplyr::filter(coffee_p, year > 2000, year < 2020)
unique(coffee_p$country_name)
summary(coffee_p$ppslc.Coffee)
summary(coffee_p$ip_Coffee)
summary(coffee_p$`ip_Coffee, Robusta`)
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
