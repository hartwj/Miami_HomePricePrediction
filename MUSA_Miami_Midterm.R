
# --- Setup: Libraries ----
# Regex parsing package
install.packages("stringr", dependencies = TRUE)
library(stringr)

library(tidyverse)
library(ggplot2)
library(sf)
library(rmarkdown)
library(kableExtra)
library(viridis)
library(tidycensus)
library(osmdata)
library(mapview)
library(lubridate)
library(ggcorrplot)

library(spdep)
library(geosphere)
library(caret)
library(ckanr)
library(FNN)
library(geosphere)
library(osmdata)
library(grid)
library(gridExtra)
library(ggstance)
library(jtools)     
library(broom)
# library(tufte)    #excluding for now..weird errors
library(readr)

# --- Setup: Aesthetics & Functions ----
## Aesthetics
mapTheme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black"),
    plot.title = element_text(size = 14,colour = "black"),
    plot.subtitle=element_text(face="italic"),
    plot.caption=element_text(hjust=0),
    axis.ticks = element_blank(),
    panel.background = element_blank(),axis.title = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2)
  )
}

plotTheme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black"),
    plot.title = element_text(size = 14,colour = "black"),
    plot.subtitle = element_text(face="italic"),
    plot.caption = element_text(hjust=0),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_line("grey80", size = 0.1),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    strip.background = element_rect(fill = "grey80", color = "white"),
    strip.text = element_text(size=12),
    axis.title = element_text(size=12),
    axis.text = element_text(size=10),
    plot.background = element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(colour = "black", face = "italic"),
    legend.text = element_text(colour = "black", face = "italic"),
    strip.text.x = element_text(size = 14)
  )
}

palette5 <- c("#25CB10", "#5AB60C", "#8FA108",   "#C48C04", "#FA7800")

qBr <- function(df, variable, rnd) {
  if (missing(rnd)) {
    as.character(quantile(round(df[[variable]],0),
                          c(.01,.2,.4,.6,.8), na.rm=T))
  } else if (rnd == FALSE | rnd == F) {
    as.character(formatC(quantile(df[[variable]]), digits = 3),
                 c(.01,.2,.4,.6,.8), na.rm=T)
  }
}

q5 <- function(variable) {as.factor(ntile(variable, 5))}

## Functions


nn_function <- function(measureFrom,measureTo,k) {
  measureFrom_Matrix <- as.matrix(measureFrom)
  measureTo_Matrix <- as.matrix(measureTo)
  nn <-   
    get.knnx(measureTo, measureFrom, k)$nn.dist
  output <-
    as.data.frame(nn) %>%
    rownames_to_column(var = "thisPoint") %>%
    gather(points, point_distance, V1:ncol(.)) %>%
    arrange(as.numeric(thisPoint)) %>%
    group_by(thisPoint) %>%
    dplyr::summarize(pointDistance = mean(point_distance)) %>%
    arrange(as.numeric(thisPoint)) %>% 
    dplyr::select(-thisPoint) %>%
    pull()
  
  return(output)  
}

multipleRingBuffer <- function(inputPolygon, maxDistance, interval) {
  distances <- seq(0, maxDistance, interval)
  distancesCounter <- 2
  numberOfRings <- floor(maxDistance / interval)
  numberOfRingsCounter <- 1
  allRings <- data.frame()
  
  while (numberOfRingsCounter <= numberOfRings) {
    if(distances[distancesCounter] < 0 & distancesCounter == 2){
      buffer1 <- st_buffer(inputPolygon, distances[distancesCounter])
      buffer1_ <- st_difference(inputPolygon, buffer1)
      thisRing <- st_cast(buffer1_, "POLYGON")
      thisRing <- as.data.frame(thisRing[,ncol(thisRing)])
      thisRing$distance <- distances[distancesCounter]
    }
    
    else if(distances[distancesCounter] < 0 & distancesCounter > 2) {
      buffer1 <- st_buffer(inputPolygon, distances[distancesCounter])
      buffer2 <- st_buffer(inputPolygon, distances[distancesCounter-1])
      thisRing <- st_difference(buffer2,buffer1)
      thisRing <- st_cast(thisRing, "POLYGON")
      thisRing <- as.data.frame(thisRing$geometry)
      thisRing$distance <- distances[distancesCounter]
    }
    
    else {
      buffer1 <- st_buffer(inputPolygon, distances[distancesCounter])
      buffer1_ <- st_buffer(inputPolygon, distances[distancesCounter-1])
      thisRing <- st_difference(buffer1,buffer1_)
      thisRing <- st_cast(thisRing, "POLYGON")
      thisRing <- as.data.frame(thisRing[,ncol(thisRing)])
      thisRing$distance <- distances[distancesCounter]
    }  
    
    allRings <- rbind(allRings, thisRing)
    distancesCounter <- distancesCounter + 1
    numberOfRingsCounter <- numberOfRingsCounter + 1
  }
  
  allRings <- st_as_sf(allRings)
}


## Load census API key
census_api_key("dc04d127e79099d0fa300464507544280121fc3b", overwrite = TRUE)

# --- Part 1: Data Wrangling ----

### Reading in Home Price Data & Base Map

# Julian file path "C:/Users/12156/Documents/GitHub/Miami/studentsData.geojson"
# JZhou file path "/Users/julianazhou/Documents/GitHub/Miami/studentsData.geojson"

miamiHomes <- st_read("/Users/julianazhou/Documents/GitHub/Miami/studentsData.geojson")
miamiHomes.sf    <- miamiHomes %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')

mapview(miamiHomes.sf[1,])
glimpse(miamiHomes.sf)

# Specify saleYear as integer
miamiHomes.sf$saleYear <- as.integer(miamiHomes.sf$saleYear)

### Read in base map
miami.base <- 
  st_read("https://opendata.arcgis.com/datasets/5ece0745e24b4617a49f2e098df8117f_0.geojson") %>%
  st_transform('ESRI:102658') %>%
  filter(NAME == "MIAMI BEACH" | NAME == "MIAMI") %>%
  st_union()

# Create border box around Miami base map to pull data from OSM
xmin = st_bbox(miami.base)[[1]]
ymin = st_bbox(miami.base)[[2]]
xmax = st_bbox(miami.base)[[3]]  
ymax = st_bbox(miami.base)[[4]]

ggplot() +
  geom_sf(data=miami.base, fill="black") +
  geom_sf(data=st_as_sfc(st_bbox(miami.base)), colour="red", fill=NA) 

#bars <- opq(bbox = c(xmin, ymin, xmax, ymax)) %>% #fix
#  add_osm_feature(key = 'amenity', value = c("bar", "pub", "restaurant")) %>%
#  osmdata_sf()

### Joining shoreline distance to MiamiHomes.sf

# Read in shoreline shapefile
shoreline <- st_read('https://opendata.arcgis.com/datasets/58386199cc234518822e5f34f65eb713_0.geojson') %>% 
  st_transform('ESRI:102658')

# Clip shoreline by base map
shoreline <- st_intersection(shoreline, miami.base)

# Transform shoreline to points
shoreline.point <- st_cast(shoreline,"POINT") 

# Creating distance to shore feature
miamiHomes.sf <- miamiHomes.sf %>%
  mutate(Shore1 = nn_function(st_c(st_centroid(miamiHomes.sf)),
                              st_c(st_centroid(shoreline.point)),1))

# Map of Distance to Shore
ggplot() + geom_sf(data=miami.base) + 
  geom_sf(data=miamiHomes.sf, aes(colour=Shore1)) + 
  scale_colour_viridis()

### Joining Neighborhoods to Miami.sf
# Load Census data
tracts <- 
  get_acs(geography = "tract", variables = c("B25026_001E","B02001_002E",
                                             "B19013_001E","B25058_001E",
                                             "B06012_002E"), 
  year=2018, state= 12, county= 086, geometry=T, output="wide") %>%
  st_transform('ESRI:102658') %>%
  rename(TotalPop = B25026_001E, 
         Whites = B02001_002E,
         MedHHInc = B19013_001E, 
         MedRent = B25058_001E,
         TotalPoverty = B06012_002E) %>%
  dplyr::select(-NAME, -starts_with("B")) %>%
  mutate(pctWhite = ifelse(TotalPop > 0, Whites / TotalPop,0),
         pctPoverty = ifelse(TotalPop > 0, TotalPoverty / TotalPop, 0)) %>%
  dplyr::select(-Whites, -TotalPoverty) 

miamiHomes.sf <- st_join(miamiHomes.sf, tracts, join = st_within)

### Parsing XF variables
# Make all XF# text lowercase and combine into one string variable
miamiHomes.sf <- miamiHomes.sf %>%
  mutate(XF1 = tolower(XF1)) %>%
  mutate(XF2 = tolower(XF2)) %>%
  mutate(XF3 = tolower(XF3)) %>%
  mutate(XF_all = paste(XF1,XF2,XF3,sep = " "))

# Use regex to create patio, pool, and fence dummy variables
miamiHomes.sf <- miamiHomes.sf %>%
  mutate(Pool = as.integer(str_detect(XF_all,"pool"))) %>%
  mutate(Fence = as.integer(str_detect(XF_all,"fence"))) %>%
  mutate(Patio = as.integer(str_detect(XF_all,"patio")))
  
## Cleaning miamiHomes.sf
miamiHomesClean.sf <- 
  miamiHomes.sf %>%
  dplyr::select(Folio, SalePrice, Property.Zip, Mailing.Zip, Property.City, AdjustedSqFt,
                LotSize, Bed, Bath, Stories, YearBuilt, EffectiveYearBuilt,
                LivingSqFt, ActualSqFt, toPredict, Shore1, GEOID, TotalPop, MedHHInc, MedRent, pctWhite, 
                pctPoverty, geometry) %>%
  mutate(Age = 2018 - YearBuilt) %>%
  mutate(EffectiveAge = 2018 - EffectiveYearBuilt) # I literally could not get this to work without 2018

## Runing a Correlation Matrix to find interesting variables
miamiHomes.train <- miamiHomesClean.sf %>% 
  st_drop_geometry() %>%
  filter(toPredict == 0)

miamiHomes.test <- miamiHomesClean.sf %>% 
  st_drop_geometry() %>%
  filter(toPredict == 1)


numericVars <- 
  select_if(miamiHomes.train, is.numeric) %>% na.omit()


ggcorrplot(
  round(cor(numericVars), 1), 
  p.mat = cor_pmat(numericVars),
  colors = c("#25CB10", "white", "#FA7800"),
  type="lower",
  insig = "blank") +  
  labs(title = "Correlation across numeric variables") 



# --- Part 3. Feature Engineering ----



cor.test(miamiHomes.train$AdjustedSqFt, miamiHomes.train$SalePrice, method = "pearson")


hist(miamiHomes.train$SalePrice)
ggplot(filter(miamiHomes.train, SalePrice <= 2000000), aes(y=SalePrice, x = AdjustedSqFt)) +
  geom_point() +
  geom_smooth(method = "lm")

## Univarite Regression
Reg1 <- lm(miamiHomes.train$SalePrice ~ ., data = miamiHomes.train %>%
  dplyr::select(Shore1))

Reg2 <- lm(miamiHomes.train$SalePrice ~ ., data = miamiHomes.train %>%
             dplyr::select(-Mailing.Zip, -Property.Zip, -GEOID))

summary(Reg1)
summary(Reg2)
summ(Reg1)

#GEOID R2 = .3, MailingZip =.4, PropertyZip =.9

# -----Part 4 - Train/Test Split -----------


# set random seed
set.seed(3171)

# get index for training sample
inTrain <- caret::createDataPartition(
  y = miamiHomes.train$SalePrice, 
  p = .60, list = FALSE)
# split data into training and test
miami.training <- miamiHomes.train[inTrain,] 
miami.test     <- miamiHomes.train[-inTrain,]  

# Regression  
reg3 <- lm(SalePrice ~ ., data = miami.test %>%
             dplyr::select(-Mailing.Zip, -Property.Zip, -GEOID))
#reg3 regularly gets adjusted R2 of 0.80 to 0.81


summary(reg3)

reg4 <- lm(SalePrice ~ ., data = miami.test %>%
             dplyr::select(-Mailing.Zip, -Property.Zip, -GEOID,
                           -AdjustedSqFt, -YearBuilt, -LivingSqFt,
                           -toPredict, -TotalPop, -MedHHInc,
                           -EffectiveAge,-EffectiveYearBuilt))
summary(reg4)
#reg4 improves to .81

reg5 <- lm(SalePrice ~ ., data = miami.test %>%
             dplyr::select(-Property.Zip, -GEOID,
                           -AdjustedSqFt, -YearBuilt, -LivingSqFt,
                           -toPredict, -TotalPop, -MedHHInc,
                           -EffectiveAge,-EffectiveYearBuilt))

summary(reg5)
#GEOID = .8337, property zip = .99, mailing zip 0.845

#The neighborhood variables clearly make our model better, but could still be overfitting


# --------K Fold Cross Validation --------------------

fitControl <- trainControl(method = "cv", 
                           number = 10,
                           # savePredictions differs from book
                           savePredictions = TRUE)

set.seed(722)
# crimes.buffer feature added
# for k-folds CV
reg.cv1 <- 
  train(SalePrice ~ ., data = miamiHomes.train %>%
          dplyr::select(-Property.Zip, -GEOID, -Mailing.Zip,
                        -AdjustedSqFt, -YearBuilt, -LivingSqFt,
                        -toPredict, -TotalPop, -MedHHInc,
                        -EffectiveAge,-EffectiveYearBuilt), 
        method = "lm", 
        trControl = fitControl, 
        na.action = na.pass)

reg.cv1

#No neighborhoods R2=0.788, MAE=525,277

reg.cv2 <- 
  train(SalePrice ~ ., data = miamiHomes.train %>%
          dplyr::select(-Property.Zip,-Mailing.Zip,
                        -AdjustedSqFt, -YearBuilt, -LivingSqFt,
                        -toPredict, -TotalPop, -MedHHInc,
                        -EffectiveAge,-EffectiveYearBuilt), 
        method = "lm", 
        trControl = fitControl, 
        na.action = na.pass)

reg.cv2

#mailing zip R2=0.735, MAE=537,844    GEOID R2=0.796, MAE=523,302  Property zip broke R lol

#starting from scratch, lets build a model
#our best model was reg.cv2

set.seed(724)
reg.cv3 <- 
  train(SalePrice ~ ., data = miamiHomes.train %>%
          dplyr::select(SalePrice, AdjustedSqFt, LotSize, GEOID, Bed, Bath, stories, Shore1,
                        Age, TotalPop, MedHHInc, MedRent, pctWhite, pctPoverty), 
        method = "lm", 
        trControl = fitControl, 
        na.action = na.pass)

reg.cv3







# finding counts by group
plot1 <- group_by(bostonCrimes, OFFENSE_CODE_GROUP) %>%
  summarize(count = n()) %>%
  arrange(-count) %>% top_n(10)

ggplot(plot1,aes(x=reorder(OFFENSE_CODE_GROUP,count),y=count))+
  geom_bar(stat="identity") +
  coord_flip()+
  theme_bw()

# ggplot, reorder

# Mapping data
ggplot() +
  geom_sf(data = nhoods, fill = "grey40") +
  geom_sf(data = boston.sf, aes(colour = q5(PricePerSq)), 
          show.legend = "point", size = .75) +
  scale_colour_manual(values = palette5,
                      labels=qBr(boston,"PricePerSq"),
                      name="Quintile\nBreaks") +
  labs(title="Price Per Square Foot, Boston") +
  mapTheme()

names(boston.sf)

# Cleaning Crime Data
bostonCrimes.sf <-
  bostonCrimes %>%
  filter(OFFENSE_CODE_GROUP == "Aggravated Assault",
         Lat > -1) %>%
  dplyr::select(Lat, Long) %>%
  na.omit() %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102286') %>%
  distinct()


# Counts of crime per buffer of house sale
boston.sf$crimes.Buffer =
  st_buffer(boston.sf, 660) %>% 
  aggregate(mutate(bostonCrimes.sf, counter = 1),., sum) %>%
  pull(counter)

ggplot() + geom_sf(data = nhoods, fill = "grey40") +
  stat_density2d(data = data.frame(st_coordinates(bostonCrimes.sf)), 
                 aes(X, Y, fill = ..level.., alpha = ..level..),
                 size = 0.01, bins = 40, geom = 'polygon') +
  scale_fill_gradient(low = "#25CB10", high = "#FA7800", name = "Density") +
  scale_alpha(range = c(0.00, 0.35), guide = FALSE) +
  labs(title = "Density of Aggravated Assaults, Boston") +
  mapTheme()

## Nearest Neighbor Feature
st_c <- st_coordinates

## measure from, measure to, k
boston.sf <-
  boston.sf %>% 
  mutate(
    crime_nn1 = nn_function(st_c(boston.sf), st_c(bostonCrimes.sf), 1),
    crime_nn2 = nn_function(st_c(boston.sf), st_c(bostonCrimes.sf), 2), 
    crime_nn3 = nn_function(st_c(boston.sf), st_c(bostonCrimes.sf), 3), 
    crime_nn4 = nn_function(st_c(boston.sf), st_c(bostonCrimes.sf), 4), 
    crime_nn5 = nn_function(st_c(boston.sf), st_c(bostonCrimes.sf), 5)) 

## Plot NN count over space - Should increase or decrease?
boston.sf.plot <- boston.sf %>% 
  st_drop_geometry() %>% 
  dplyr::select(Parcel_No, starts_with("crime_")) %>% 
  tidyr::pivot_longer(cols = -Parcel_No, names_to = "crime_nn")

ggplot(boston.sf.plot, aes(x = crime_nn, y = value, group = Parcel_No)) +
  geom_line(alpha = 0.05, color = "firebrick") +
  theme_bw()

# boston.sf   <- st_read(file.path(data_src,"boston_sf_Ch1_wrangled.geojson"))

## Home Features cor
st_drop_geometry(boston.sf) %>% 
  mutate(Age = 2015 - YR_BUILT) %>%
  dplyr::select(SalePrice, LivingArea, Age, GROSS_AREA) %>%
  filter(SalePrice <= 1000000, Age < 500) %>%
  gather(Variable, Value, -SalePrice) %>% 
  ggplot(aes(Value, SalePrice)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se=F, colour = "#FA7800") +
  facet_wrap(~Variable, ncol = 3, scales = "free") +
  labs(title = "Price as a function of continuous variables") +
  plotTheme()

## Crime cor
st_drop_geometry(boston.sf) %>% 
  mutate(Age = 2015 - YR_BUILT) %>%
  dplyr::select(SalePrice, starts_with("crime_")) %>%
  filter(SalePrice <= 1000000) %>%
  gather(Variable, Value, -SalePrice) %>% 
  ggplot(aes(Value, SalePrice)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se=F, colour = "#FA7800") +
  facet_wrap(~Variable, nrow = 1, scales = "free") +
  labs(title = "Price as a function of continuous variables") +
  plotTheme()

### Corr matrix
numericVars <- 
  select_if(st_drop_geometry(boston.sf), is.numeric) %>% na.omit()

ggcorrplot(
  round(cor(numericVars), 1), 
  p.mat = cor_pmat(numericVars),
  colors = c("#25CB10", "white", "#FA7800"),
  type="lower",
  insig = "blank") +  
  labs(title = "Correlation across numeric variables") 


## Univarite correlation vs. regression
cor.test(boston$LivingArea, boston$SalePrice, method = "pearson")

ggplot(filter(boston, SalePrice <= 2000000), aes(y=SalePrice, x = LivingArea)) +
  geom_point() +
  geom_smooth(method = "lm")

## Univarite Regrssion
livingReg <- lm(SalePrice ~ LivingArea, data = boston)

summary(livingReg)
summ(livingReg)

## Prediction example
new_LivingArea = 4000
# "by hand"
157968.32 + 216.54 * new_LivingArea
# predict() function
predict(livingReg, newdata = data.frame(LivingArea = 4000))

## plot of marginal regression response
effect_plot(livingReg, pred = LivingArea, interval = TRUE, plot.points = TRUE)

## Multivariate Regression
reg1 <- lm(SalePrice ~ ., data = st_drop_geometry(boston.sf) %>% 
             dplyr::select(SalePrice, LivingArea, Style, 
                           GROSS_AREA, R_TOTAL_RM, NUM_FLOORS,
                           R_BDRMS, R_FULL_BTH, R_HALF_BTH, 
                           R_KITCH, R_AC, R_FPLACE))
summ(reg1)
broom::glance(reg1)

## Plot of marginal response
effect_plot(reg1, pred = R_BDRMS, interval = TRUE, plot.points = TRUE)

## Plot coefficients
plot_summs(reg1)

## plot multiple model coeffs
plot_summs(reg1, livingReg)


## What is the Coefficient of LivingArea when Average Distance to 2-nearest crimes are considered?
## Build a regression with LivingArea and crime_nn2
## report regression coefficient for LivingArea
## Is it different? Why?

## External model validation


# Evaluating a model through its predictions on *new* data.

# set random seed
set.seed(31357)

# get index for training sample
inTrain <- caret::createDataPartition(
  y = boston.sf$SalePrice, 
  p = .60, list = FALSE)
# split data into training and test
boston.training <- boston.sf[inTrain,] 
boston.test     <- boston.sf[-inTrain,]  

# Regression  
reg2 <- lm(SalePrice ~ ., data = st_drop_geometry(boston.training) %>% 
             dplyr::select(SalePrice, LivingArea, 
                           GROSS_AREA, R_TOTAL_RM, NUM_FLOORS,
                           R_BDRMS, R_FULL_BTH, R_HALF_BTH, 
                           R_KITCH, R_FPLACE))

# Run this a number of times to see Adjusted R2
summary(reg2)

## predicting on new data
reg2_predict <- predict(reg2, newdata = boston.test)

## Measure generalizability?
## Mean Square Error train and test
rmse.train <- caret::MAE(predict(reg2), boston.training$SalePrice)
rmse.test  <- caret::MAE(reg2_predict, boston.test$SalePrice)

cat("Train MAE: ", as.integer(rmse.train), " \n","Test MAE: ", as.integer(rmse.test))

preds.train <- data.frame(pred   = predict(reg2),
                          actual = boston.training$SalePrice,
                          source = "training data")
preds.test  <- data.frame(pred   = reg2_predict,
                          actual = boston.test$SalePrice,
                          source = "testing data")
preds <- rbind(preds.train, preds.test)

ggplot(preds, aes(x = pred, y = actual, color = source)) +
  geom_point() +
  geom_smooth(method = "lm", color = "green") +
  geom_abline(color = "orange") +
  coord_equal() +
  theme_bw() +
  facet_wrap(~source, ncol = 2) +
  labs(title = "Comparing predictions to actual values",
       x = "Predicted Value",
       y = "Actual Value") +
  theme(
    legend.position = "none"
  )


## Cross-validation

# use caret package cross-validation method
fitControl <- trainControl(method = "cv", 
                           number = 10,
                           # savePredictions differs from book
                           savePredictions = TRUE)

set.seed(717)
# crimes.buffer feature added
# for k-folds CV
reg.cv <- 
  train(SalePrice ~ ., data = st_drop_geometry(boston.sf) %>% 
          dplyr::select(SalePrice, LivingArea,  
                        GROSS_AREA, R_TOTAL_RM, NUM_FLOORS,
                        R_BDRMS, R_FULL_BTH, R_HALF_BTH, 
                        R_KITCH, R_FPLACE, crimes.Buffer), 
        method = "lm", 
        trControl = fitControl, 
        na.action = na.pass)

reg.cv

reg.cv$resample

reg.cv$resample %>% 
  pivot_longer(-Resample) %>% 
  mutate(name = as.factor(name)) %>% 
  ggplot(., aes(x = name, y = value, color = name)) +
  geom_jitter(width = 0.1) +
  facet_wrap(~name, ncol = 3, scales = "free") +
  theme_bw() +
  theme(
    legend.position = "none"
  )

# How do we interpret the spread of values between the folds?
# extract predictions from CV object
cv_preds <- reg.cv$pred
# compare number of observations between data sets
nrow(boston.sf)
nrow(cv_preds)

## Create dataset with "out of fold" predictions and original data
map_preds <- boston.sf %>% 
  rowid_to_column(var = "rowIndex") %>% 
  left_join(cv_preds, by = "rowIndex") %>% 
  mutate(SalePrice.AbsError = abs(pred - SalePrice)) %>% 
  cbind(st_coordinates(.))
# weird CRS fix to boston.sf
st_crs(map_preds) <- st_crs(nhoods)

# plot errors on a map
ggplot() +
  geom_sf(data = nhoods, fill = "grey40") +
  geom_sf(data = map_preds, aes(colour = q5(SalePrice.AbsError)),
          show.legend = "point", size = 1) +
  scale_colour_manual(values = palette5,
                      labels=qBr(map_preds,"SalePrice.AbsError"),
                      name="Quintile\nBreaks") +
  labs(title="Absolute sale price errors on the OOF set",
       subtitle = "OOF = 'Out Of Fold'") +
  mapTheme()

## Spatial Correlation of Errors

# Going back to single random data split to align with book

inTrain <- createDataPartition(
  y = paste(boston.sf$Name, boston.sf$NUM_FLOORS.cat, 
            boston.sf$Style, boston.sf$R_AC), 
  p = .60, list = FALSE)
boston.training <- boston.sf[inTrain,] 
boston.test <- boston.sf[-inTrain,]  

reg.training <- 
  lm(SalePrice ~ ., data = as.data.frame(boston.training) %>% 
       dplyr::select(SalePrice, LivingArea, Style, 
                     GROSS_AREA, NUM_FLOORS.cat,
                     R_BDRMS, R_FULL_BTH, R_HALF_BTH, 
                     R_KITCH, R_AC, R_FPLACE, crimes.Buffer))

boston.test <-
  boston.test %>%
  mutate(Regression = "Baseline Regression",
         SalePrice.Predict = predict(reg.training, boston.test),
         SalePrice.Error = SalePrice.Predict - SalePrice,
         SalePrice.AbsError = abs(SalePrice.Predict - SalePrice),
         SalePrice.APE = (abs(SalePrice.Predict - SalePrice)) / SalePrice.Predict)%>%
  filter(SalePrice < 5000000) 

k_nearest_neighbors = 5
#prices
coords <- st_coordinates(boston.sf) 
# k nearest neighbors
neighborList <- knn2nb(knearneigh(coords, k_nearest_neighbors))
spatialWeights <- nb2listw(neighborList, style="W")
boston.sf$lagPrice <- lag.listw(spatialWeights, boston.sf$SalePrice)

#errors
coords.test <-  st_coordinates(boston.test) 
neighborList.test <- knn2nb(knearneigh(coords.test, k_nearest_neighbors))
spatialWeights.test <- nb2listw(neighborList.test, style="W")
boston.test$lagPriceError <- lag.listw(spatialWeights.test, boston.test$SalePrice.AbsError)

ggplot(boston.sf, aes(x=lagPrice, y=SalePrice)) +
  geom_point(colour = "#FA7800") +
  geom_smooth(method = "lm", se = FALSE, colour = "#25CB10") +
  labs(title = "Price as a function of the spatial lag of price",
       caption = "Public Policy Analytics, Figure 6.6",
       x = "Spatial lag of price (Mean price of 5 nearest neighbors)",
       y = "Sale Price") +
  plotTheme()

ggplot(boston.test, aes(x=lagPriceError, y=SalePrice)) +
  geom_point(colour = "#FA7800") +
  geom_smooth(method = "lm", se = FALSE, colour = "#25CB10") +
  labs(title = "Error as a function of the spatial lag of price",
       caption = "",
       x = "Spatial lag of errors (Mean error of 5 nearest neighbors)",
       y = "Sale Price") +
  plotTheme()

## Moran's I - A measure of spatial correlation

moranTest <- moran.mc(boston.test$SalePrice.AbsError, 
                      spatialWeights.test, nsim = 999)

ggplot(as.data.frame(moranTest$res[c(1:999)]), aes(moranTest$res[c(1:999)])) +
  geom_histogram(binwidth = 0.01) +
  geom_vline(aes(xintercept = moranTest$statistic), colour = "#FA7800",size=1) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(title="Observed and permuted Moran's I",
       subtitle= "Observed Moran's I in orange",
       x="Moran's I",
       y="Count",
       caption="Public Policy Analytics, Figure 6.8") +
  plotTheme()

## Errors by group

nhood_sum <- boston.test %>% 
  group_by(Name) %>%
  summarize(meanPrice = mean(SalePrice, na.rm = T),
            meanPrediction = mean(SalePrice.Predict, na.rm = T),
            meanMAE = mean(SalePrice.AbsError, na.rm = T))

nhood_sum %>% 
  st_drop_geometry %>%
  arrange(desc(meanMAE)) %>% 
  kable() %>% kable_styling()

## Back to CV predictions for the moment. In the split sample predictions we only get predictions for the nhoods in the test sample; the rest are NA. In the CV split, we get a prediction on each nhood. That makes for a more complete looking map.

map_preds_sum <- map_preds %>% 
  group_by(Name) %>% 
  summarise(meanMAE = mean(SalePrice.AbsError))

ggplot() +
  geom_sf(data = nhoods %>% 
            left_join(st_drop_geometry(map_preds_sum), by = "Name"),
          aes(fill = q5(meanMAE))) +
  scale_fill_manual(values = palette5,
                    labels=qBr(nhood_sum,"meanMAE"),
                    name="Quintile\nBreaks") +
  mapTheme() +
  labs(title="Absolute sale price errors on the OOF set by Neighborhood")

## Neighborhood Fixed Effect
reg.nhood <- lm(SalePrice ~ ., data = as.data.frame(boston.training) %>% 
                  dplyr::select(Name, SalePrice, LivingArea, 
                                Style, GROSS_AREA, NUM_FLOORS.cat,
                                R_BDRMS, R_FULL_BTH, R_HALF_BTH, 
                                R_KITCH, R_AC, R_FPLACE,crimes.Buffer))

boston.test.nhood <-
  boston.test %>%
  mutate(Regression = "Neighborhood Effects",
         SalePrice.Predict = predict(reg.nhood, boston.test),
         SalePrice.Error = SalePrice - SalePrice.Predict,
         SalePrice.AbsError = abs(SalePrice - SalePrice.Predict),
         SalePrice.APE = (abs(SalePrice - SalePrice.Predict)) / SalePrice) %>%
  filter(SalePrice < 5000000)

bothRegressions <- 
  rbind(
    dplyr::select(boston.test, starts_with("SalePrice"), Regression, Name) %>%
      mutate(lagPriceError = lag.listw(spatialWeights.test, SalePrice.Error)),
    dplyr::select(boston.test.nhood, starts_with("SalePrice"), Regression, Name) %>%
      mutate(lagPriceError = lag.listw(spatialWeights.test, SalePrice.Error)))    

st_drop_geometry(bothRegressions) %>%
  gather(Variable, Value, -Regression, -Name) %>%
  filter(Variable == "SalePrice.AbsError" | Variable == "SalePrice.APE") %>%
  group_by(Regression, Variable) %>%
  summarize(meanValue = mean(Value, na.rm = T)) %>%
  spread(Variable, meanValue) %>%
  kable() %>%
  kable_styling("striped", full_width = F) %>%
  row_spec(1, color = "black", background = "#25CB10") %>%
  row_spec(2, color = "black", background = "#FA7800")

### The key question to be able to answer from this is:

# Which model better captures the spatial process of home sales prices and why?

# model coefficents for each Neighborhood
tidy(reg.nhood) %>% 
  filter(str_detect(term, "Name")) %>% 
  kable() %>% 
  kable_styling()


# OSM Data

#Study area base
miamiBound <- st_read("/Users/annaduan/Documents/GitHub/2_Miami\ Prediction/Raw\ Data/Municipal_Boundary.geojson") %>%
  filter(NAME == "MIAMI BEACH" | NAME == "MIAMI") %>%
  st_union() %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')

#base for osm (not projected so that it works)
miamiBoundOSM <- st_read("/Users/annaduan/Documents/GitHub/2_Miami\ Prediction/Raw\ Data/Municipal_Boundary.geojson") %>%
  filter(NAME == "MIAMI BEACH" | NAME == "MIAMI") %>%
  st_union()

#OSM bounding box (used the OSM - specific, non-projected, base)
xmin = st_bbox(miamiBoundOSM)[[1]]
ymin = st_bbox(miamiBoundOSM)[[2]]
xmax = st_bbox(miamiBoundOSM)[[3]]  
ymax = st_bbox(miamiBoundOSM)[[4]]

#bars, restaurants, shopping
foodBev <- opq(bbox = c(xmin, ymin, xmax, ymax)) %>% 
  add_osm_feature(key = 'amenity', value = c("bar","pub","restaurant","cafe")) %>%
  osmdata_xml(filename = 'foodBev.osm')

foodBev <- sf::st_read('foodBev.osm', layer = 'points') %>%
  st_as_sf(coords = c("LON", "LAT"), crs = EPSG:3857, agr = "constant") %>% #EPSG:3857 is the projection that most OSM data is in
  st_transform('ESRI:102658')