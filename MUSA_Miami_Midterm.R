
# --- Setup: Libraries ----
# Regex parsing package
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
library(raster)
library(stringr)
library(stargazer)
library(ggpubr)

library(rgeos)
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
s
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

multipleRingBuffer <- function(inputPolygon, maxDistance, interval) 
{
  #create a list of distances that we'll iterate through to create each ring
  distances <- seq(0, maxDistance, interval)
  #we'll start with the second value in that list - the first is '0'
  distancesCounter <- 2
  #total number of rings we're going to create
  numberOfRings <- floor(maxDistance / interval)
  #a counter of number of rings
  numberOfRingsCounter <- 1
  #initialize an otuput data frame (that is not an sf)
  allRings <- data.frame()
  
  #while number of rings  counteris less than the specified nubmer of rings
  while (numberOfRingsCounter <= numberOfRings) 
  {
    #if we're interested in a negative buffer and this is the first buffer
    #(ie. not distance = '0' in the distances list)
    if(distances[distancesCounter] < 0 & distancesCounter == 2)
    {
      #buffer the input by the first distance
      buffer1 <- st_buffer(inputPolygon, distances[distancesCounter])
      #different that buffer from the input polygon to get the first ring
      buffer1_ <- st_difference(inputPolygon, buffer1)
      #cast this sf as a polygon geometry type
      thisRing <- st_cast(buffer1_, "POLYGON")
      #take the last column which is 'geometry'
      thisRing <- as.data.frame(thisRing[,ncol(thisRing)])
      #add a new field, 'distance' so we know how far the distance is for a give ring
      thisRing$distance <- distances[distancesCounter]
    }
    
    
    #otherwise, if this is the second or more ring (and a negative buffer)
    else if(distances[distancesCounter] < 0 & distancesCounter > 2) 
    {
      #buffer by a specific distance
      buffer1 <- st_buffer(inputPolygon, distances[distancesCounter])
      #create the next smallest buffer
      buffer2 <- st_buffer(inputPolygon, distances[distancesCounter-1])
      #This can then be used to difference out a buffer running from 660 to 1320
      #This works because differencing 1320ft by 660ft = a buffer between 660 & 1320.
      #bc the area after 660ft in buffer2 = NA.
      thisRing <- st_difference(buffer2,buffer1)
      #cast as apolygon
      thisRing <- st_cast(thisRing, "POLYGON")
      #get the last field
      thisRing <- as.data.frame(thisRing$geometry)
      #create the distance field
      thisRing$distance <- distances[distancesCounter]
    }
    
    #Otherwise, if its a positive buffer
    else 
    {
      #Create a positive buffer
      buffer1 <- st_buffer(inputPolygon, distances[distancesCounter])
      #create a positive buffer that is one distance smaller. So if its the first buffer
      #distance, buffer1_ will = 0. 
      buffer1_ <- st_buffer(inputPolygon, distances[distancesCounter-1])
      #difference the two buffers
      thisRing <- st_difference(buffer1,buffer1_)
      #cast as a polygon
      thisRing <- st_cast(thisRing, "POLYGON")
      #geometry column as a data frame
      thisRing <- as.data.frame(thisRing[,ncol(thisRing)])
      #add teh distance
      thisRing$distance <- distances[distancesCounter]
    }  
    
    #rbind this ring to the rest of the rings
    allRings <- rbind(allRings, thisRing)
    #iterate the distance counter
    distancesCounter <- distancesCounter + 1
    #iterate the number of rings counter
    numberOfRingsCounter <- numberOfRingsCounter + 1
  }
  
  #convert the allRings data frame to an sf data frame
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

# Specify saleYear as integer
miamiHomes.sf$saleYear <- as.integer(miamiHomes.sf$saleYear)

### Read in base map
miami.base <- 
  st_read("https://opendata.arcgis.com/datasets/5ece0745e24b4617a49f2e098df8117f_0.geojson") %>%
  st_transform('ESRI:102658') %>%
  filter(NAME == "MIAMI BEACH" | NAME == "MIAMI") %>%
  st_union()

### Read in shoreline shapefile
shoreline <- st_read('https://opendata.arcgis.com/datasets/58386199cc234518822e5f34f65eb713_0.geojson') %>% 
  st_transform('ESRI:102658')

# Clip shoreline by base map & transform to points
shoreline <- st_intersection(shoreline, miami.base)
shoreline.point <- st_cast(shoreline,"POINT") 


### 2018 - 5-Year ACS Data by Census Tract
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

### Highways & Major Road Data
roads <- 
  st_read("https://opendata.arcgis.com/datasets/8bc234275dc749329c4e242abcfc5a0f_0.geojson") %>%
  filter(CLASS == c('1','2')) %>%
  st_transform('ESRI:102658') 

miamiRds <- roads[miami.base,]


### Local Park Data
parks <- 
  st_read("https://opendata.arcgis.com/datasets/8c9528d3e1824db3b14ed53188a46291_0.geojson") %>%
  filter(CITY == "Miami Beach" | CITY == "Miami") %>%
  filter(TYPE == "Local") %>%
  st_transform('ESRI:102658') 

parks <- parks[miami.base,]


### Middle School
midschool <- 
  st_read("https://opendata.arcgis.com/datasets/dd2719ff6105463187197165a9c8dd5c_0.geojson") %>%
  filter(CITY == "Miami Beach" | CITY == "Miami") %>% 
  rename(midschoolID = ID)%>%
  st_transform('ESRI:102658') 

midschool <- midschool[miami.base,]

miamiHomes.rings <- st_join(miamiHomes.rings, midschool, join = st_within)


# --- Part 2: Feature Engineering ----

### Parsing XF variables (i.e., Pool/Fence/Patio)
# Make all XF# text lowercase and combine into one string variable
miamiHomes.sf <- miamiHomes.sf %>%
  mutate(XF1 = tolower(XF1)) %>%
  mutate(XF2 = tolower(XF2)) %>%
  mutate(XF3 = tolower(XF3)) %>%
  mutate(XF_all = paste(XF1,XF2,XF3,sep = " "))

# Use stringr package to create patio, pool, and fence dummy variables
miamiHomes.sf <- miamiHomes.sf %>%
  mutate(Pool = as.integer(str_detect(XF_all,"pool"))) %>%
  mutate(Fence = as.integer(str_detect(XF_all,"fence"))) %>%
  mutate(Patio = as.integer(str_detect(XF_all,"patio")))


### Adding shoreline distance
# Creating distance to shore feature
miamiHomes.sf <- miamiHomes.sf %>%
  mutate(Shore1 = nn_function(st_coordinates(st_centroid(miamiHomes.sf)),
                              st_coordinates(st_centroid(shoreline.point)),1))

miamiHomes.sf <- st_join(miamiHomes.sf, tracts, join = st_within)
miamiHomes.sf$Shore.mile <- miamiHomes.sf$Shore1/5280

### Major Roads: 1/8 mile ring buffers
# Create unioned buffer for major roads in Miami & Miami Beach
miamiRds.buffer <- st_union(st_buffer(miamiRds, 660)) %>%
  st_sf() %>%
  mutate(Legend = "Unioned Buffer")
miamiRds.buffer <- filter(miamiRds.buffer, Legend=="Unioned Buffer")

# Create 1/8 mile ring buffers
miami.rings <- multipleRingBuffer(miamiRds.buffer, 660*15, 660) %>%
  rename(road_dist = distance)

# Join home prices with ring buffer -- transform NAs to 0
miamiHomes.rings <- st_join(miamiHomes.sf, miami.rings, join = st_within) %>%
  st_sf() 
miamiHomes.rings[c("road_dist")][is.na(miamiHomes.rings[c("road_dist")])] <- 0


### Parks - 1/8 ring buffers
parks.buffer <- st_union(st_buffer(parks, 660)) %>%
  st_sf() %>%
  mutate(Legend = "Unioned Buffer")

park.rings <- multipleRingBuffer(parks.buffer, 660*8, 660) %>%
  rename(park_dist = distance)

# Join home prices with park ring buffers
miamiHomes.rings <- st_join(miamiHomes.rings, park.rings, join = st_within) %>%
  st_sf() 

miamiHomes.rings[c("park_dist")][is.na(miamiHomes.rings[c("park_dist")])] <- 0


### Cleaning miamiHomes.sf for exploratory analyses
miamiHomesClean.sf <- miamiHomes.rings %>%
  mutate(Age = saleYear - YearBuilt) %>%
  dplyr::select(Folio, SalePrice, Property.City,
                LotSize, Bed, Bath, Stories, Pool, Fence, Patio, ActualSqFt, 
                YearBuilt, EffectiveYearBuilt, Age, toPredict, Shore1, GEOID, TotalPop, 
                MedHHInc, MedRent, pctWhite, pctPoverty, road_dist, park_dist, 
                midschoolID, geometry) 

# --- Markdown: Introduction ----



# --- Markdown: Data ----
trainprice.sf <- miamiHomesClean.sf %>% 
  filter(toPredict == 0)

# Map Home Prices
allprices.sf <- full_join(trainprice.sf %>% as.data.frame(), map_preds %>% as.data.frame())%>%
  st_sf()

allprices.sf$SalePrice[allprices.sf$SalePrice == 0] <- NA   

allprices.sf$SalePrice <- ifelse(is.na(allprices.sf$SalePrice), 
                                 allprices.sf$pred, allprices.sf$SalePrice)

# All Home Prices (Predicted and Training)
ggplot() + 
  geom_sf(data=miami.base) + 
  geom_sf(data=allprices.sf, aes(color=SalePrice),
          show.legend = "line", size = 1) + 
  labs(title = "Home Prices",
       subtitle = "Miami & Miami Beach, FL",
       caption = "Figure 1.0") +
  scale_color_viridis(option="C", 
                      name = "Price ($)", 
                      limits=c(10000,1000000),
                      breaks=c(0, 250000, 500000, 750000, 1000000),
                      direction = -1,
                      begin = 0,
                      end = .95,
                      na.value = .95)

# Map of Predicted Home Prices vs. Training
tst <- ggplot() + 
  geom_sf(data=miami.base) + 
  geom_sf(data=map_preds, aes(color=pred),
          show.legend = "line", size = 1) + 
  labs(title = "Predicted Home Prices",
       caption = "Figure 1.1") +
  scale_color_viridis(option="C", 
                      name = "Price ($)", 
                      limits=c(10000,1000000),
                      breaks=c(0, 250000, 500000, 750000, 1000000),
                      direction = -1,
                      na.value = .97)
trn <- ggplot() + 
  geom_sf(data=miami.base) + 
  geom_sf(data=trainprice.sf, aes(color=SalePrice),
          show.legend = "line", size = 1) + 
  labs(title = "Training Set Sale Prices",
       caption = "Figure 1.2") +
  scale_color_viridis(option="C", 
                      name = "Price ($)", 
                      limits=c(10000,1000000),
                      breaks=c(0, 250000, 500000, 750000, 1000000),
                      direction = -1,
                      na.value = .97) 
ggarrange(tst, trn + rremove("x.text"), ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "right")

# Map of Distance to Shore
ggplot() + 
  geom_sf(data=miami.base) + 
  geom_sf(data=miamiHomes.sf, aes(colour=Shore.mile),
          show.legend = "line", size= 1) + 
  labs(title = "Distance from Shoreline",
       caption = "Figure 2.0") +
  scale_colour_viridis(name = "Distance (miles)")

### Middle School Areas
# Retrieve coordinates for school labels
midschool.pts <- sf::st_point_on_surface(midschool)
midschool.coords <- as.data.frame(sf::st_coordinates(midschool.pts))
midschool.coords$NAME <- midschool$NAME

# Map of Middle School Areas
ggplot() + 
  geom_sf(data = miami.base, fill = "lightgray", lwd = .5) +
  geom_sf(data=allprices.sf, aes(color=SalePrice),
          show.legend = "line", size = 1) + 
  geom_sf(data = midschool, fill = "#d3e9ff", alpha =.4, color = "#498cd3", lwd = .8)+
  geom_text(data = midschool.coords, aes(X, Y, label = NAME), size = 3) +
  labs(title = "Middle Schools & Home Prices",
       caption = "Figure 2.1") +
  scale_color_viridis(option="C", 
                      name = "Price ($)", 
                      limits=c(10000,1000000),
                      breaks=c(0, 250000, 500000, 750000, 1000000),
                      direction = -1,
                      na.value = .97) 

# Plot to check road ring buffers -- Don't think we need, since it wasn't relevant
ggplot() + 
  geom_sf(data = miami.base, fill = "lightgray", lwd = 1) +
  geom_sf(data = miami.rings, fill = "transparent") +
  geom_sf(data=miamiHomes.sf, color='#0000CC', size= .8) +
  geom_sf(data=miamiRds, color = "gold", size= 2) +
  labs(title = "Distance from Major Roads",
       subtitle = "1/8 Mile Increments",
       caption = "Figure 2.1") 

# Plot to check park ring buffers -- Don't think we need, since it wasn't relevant
ggplot() + 
  geom_sf(data = miami.base, fill = "lightgray", lwd = 1) +
  geom_sf(data=miamiHomes.sf, color='#0000CC', size= .8) +
  geom_sf(data = park.rings, fill = "transparent", alpha = 0.1) +
  geom_sf(data=parks, color = "gold", size= 2) +
  labs(title = "Distance from Park",
       subtitle = "1/8 Mile Increments",
       caption = "Figure 2.2")


# --- Markdown: Method ----


# --- Markdown: Results ----


# --- Part 3: Exploratory Analysis ----
## Runing a Correlation Matrix to find interesting variables
miamiHomes.train <- miamiHomesClean.sf %>% 
  filter(toPredict == 0) %>%
  filter(SalePrice <= 1000000)
miamiHomes.test <- miamiHomesClean.sf %>% 
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


cor.test(miamiHomes.train$AdjustedSqFt, miamiHomes.train$SalePrice, method = "pearson")

hist(miamiHomes.train$SalePrice)
ggplot(filter(miamiHomes.train, SalePrice <= 2000000), aes(y=SalePrice, x = AdjustedSqFt)) +
  geom_point() +
  geom_smooth(method = "lm")

# --- END HERE! ----

# --- Part 3: Exploratory Analysis ----
## Runing a Correlation Matrix to find interesting variables
miamiHomes.train <- miamiHomesClean.sf %>% 
  filter(toPredict == 0) %>%
  filter(SalePrice <= 1000000)
miamiHomes.test <- miamiHomesClean.sf %>% 
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



cor.test(miamiHomes.train$AdjustedSqFt, miamiHomes.train$SalePrice, method = "pearson")

hist(miamiHomes.train$SalePrice)
ggplot(filter(miamiHomes.train, SalePrice <= 2000000), aes(y=SalePrice, x = AdjustedSqFt)) +
  geom_point() +
  geom_smooth(method = "lm")

## Univarite Regression
Reg1 <- lm(miamiHomes.train$SalePrice ~ ., data = miamiHomes.train %>%
  dplyr::select(Shore1, TotalPop, MedHHInc, MedRent, pctWhite, pctPoverty, road_dist, park_dist, schoolID))

Reg2 <- lm(miamiHomes.train$SalePrice ~ ., data = miamiHomes.train %>%
             dplyr::select(-GEOID))

summary(Reg1)
summary(Reg2)
summ(Reg1)

stargazer(Reg1, Reg2, title="Training Set LM Results", align=TRUE, type = "html", out = "Regression_Results.htm")



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


reg2_split <- lm(SalePrice ~ ., data = miami.training %>%
             dplyr::select(-GEOID))

summary(reg2_split)


##---------------------- Calculating MAE and MAPE for a single test test

mape <- function(actual,pred){
  mape <- mean(abs((actual - pred)/actual))*100
  return (mape)
}

reg2MAE <- caret::MAE(predict(reg2_split), miami.test$SalePrice)
reg2MAPE <- mape(miami.test$SalePrice,predict(reg2_split))


cat("Train MAE: ", as.integer(reg2MAE), " \n","Test MAE: ", as.integer(reg2MAPE)) %>%
  knitr::kable()



#The neighborhood variables clearly make our model better, but could still be overfitting


# --------K Fold Cross Validation --------------------

fitControl <- trainControl(method = "cv", 
                           number = 100,
                           # savePredictions differs from book
                           savePredictions = TRUE)

set.seed(7350)
#No neighborhoods R2=0.788, MAE=525,277

reg.cv2 <- 
  train(SalePrice ~ ., data = miamiHomes.train %>%
          dplyr::select(-AdjustedSqFt, -LivingSqFt, -YearBuilt, -EffectiveYearBuilt, -TotalPop, 
                        -MedHHInc, -road_dist, -park_dist, -toPredict), 
        method = "lm", 
        trControl = fitControl, 
        na.action = na.omit)

reg.cv2

#Providing Results of Cross Validation Test
reg.cv2$results %>%
  knitr::kable()

#Creating an MAE Histogram
MAE_hist <- as.data.frame(reg.cv2$pred) %>%
  mutate(ABS_error = abs(pred-obs)) %>%
  group_by(Resample) %>%
  summarize(MAE = mean(ABS_error, na.rm = T))

ggplot(MAE_hist, aes(x=MAE))+geom_histogram()+
  labs(title="Histogram of Cross Validation MAE",
       caption="Based on 100 Folds")+
  theme_classic()
  

#Plotting Predicted Prices as a Function of Observed Prices

regCV2plot <- as.data.frame(reg.cv2$pred)


ggplot(regCV2plot,aes(obs,pred)) +
  geom_point() +
  stat_smooth(aes(obs, obs), 
              method = "lm", se = FALSE, size = 1, colour="#FA7800") + 
  stat_smooth(aes(obs, pred), 
              method = "lm", se = FALSE, size = 1, colour="#25CB10") +
  labs(title="Predicted sale price as a function of observed price",
       subtitle="Orange line represents a perfect prediction; Green line represents prediction",
       x="Sale Price",
       y="Predicted Price") +
  plotTheme() + theme(plot.title = element_text(size = 18, colour = "black")) 




#---------Outputting our CV predictions to CV--------------

# Folio, Property.City, LotSize, Bed, Bath, Stories, Pool, Fence, Patio, ActualSqFt, Age, 
# Shore1, GEOID, MedRent, pctWhite, pctPoverty
# GEOID = .8337, property zip = .99, mailing zip 0.845

# How do we interpret the spread of values between the folds?
# extract predictions from CV object
cv_preds <- reg.cv2$pred
# compare number of observations between data sets
nrow(miamiHomes.train)
nrow(cv_preds)

## Create dataset with "out of fold" predictions and original data
map_preds <- miamiHomes.test %>% 
  rowid_to_column(var = "rowIndex") %>% 
  left_join(cv_preds, by = "rowIndex") %>% 
  mutate(SalePrice.AbsError = abs(pred - SalePrice))


#output to CSV

output_preds <- map_preds %>%
  dplyr::select(pred, Folio) %>%
  mutate(team_name = "Florid-API Keys") %>%
  rename(prediction = pred)

write.csv(output_preds, "Florid-API Keys.csv")



#our best model was reg.cv2



#-----Is there a spatial correlation of errors-------

# get index for training sample
inTrain <- caret::createDataPartition(
  y = miamiHomes.train$SalePrice, 
  p = .60, list = FALSE)
# split data into training and test
miami.training <- miamiHomes.train[inTrain,] 
miami.test     <- miamiHomes.train[-inTrain,]  

#this is the dataset and variables we used in regCV2
#however I had to drop GEOID since there was a levels error
finaltraining <- lm(SalePrice ~ ., data = miami.training %>%
                      st_drop_geometry() %>%
                   dplyr::select(-GEOID, -AdjustedSqFt, -LivingSqFt, -YearBuilt, -EffectiveYearBuilt, 
                                 -TotalPop, -MedHHInc, -road_dist, -park_dist, -toPredict))

miami.test <-
  miami.test %>%
  mutate(Regression = "Baseline Regression",
         SalePrice.Predict = predict(finaltraining, miami.test),
         SalePrice.Error = SalePrice.Predict - SalePrice,
         SalePrice.AbsError = abs(SalePrice.Predict - SalePrice),
         SalePrice.APE = (abs(SalePrice.Predict - SalePrice)) / SalePrice.Predict)



k_nearest_neighbors = 5
#prices
coords <- st_coordinates(st_centroid(miamiHomesClean.sf))
# k nearest neighbors
neighborList <- knn2nb(knearneigh(coords, k_nearest_neighbors))
spatialWeights <- nb2listw(neighborList, style="W")
miamiHomesClean.sf$lagPrice <- lag.listw(spatialWeights, miamiHomesClean.sf$SalePrice)


#errors
coords.test <- st_coordinates(st_centroid(miami.test))
neighborList.test <- knn2nb(knearneigh(coords.test, k_nearest_neighbors))
spatialWeights.test <- nb2listw(neighborList.test, style="W")
miami.test$lagPriceError <- lag.listw(spatialWeights.test, miami.test$SalePrice.AbsError)

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



















### The key question to be able to answer from this is:

# Which model better captures the spatial process of home sales prices and why?

# model coefficents for each Neighborhood
tidy(reg.nhood) %>% 
  filter(str_detect(term, "Name")) %>% 
  kable() %>% 
  kable_styling()



