
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
library(caret) 
library(mapview)
library(ggspatial)


library(rgeos)
library(spdep)
library(geosphere)
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

options(scipen = 999)
# --- Setup: Aesthetics ----
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

# --- Setup: Functions ----

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


# --- Part 1: Data Wrangling ----

miamiHomes <- st_read("studentsData.geojson")
miamiHomes.sf <- miamiHomes %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')

xmin = st_bbox(miami.base)[[1]]
ymin = st_bbox(miami.base)[[2]]
xmax = st_bbox(miami.base)[[3]]  
ymax = st_bbox(miami.base)[[4]]

bbox = c(xmin, ymin, xmax, ymax)

# Specify saleYear as integer
miamiHomes.sf$saleYear <- as.integer(miamiHomes.sf$saleYear)
miamiHomes.sf$ID <- seq.int(nrow(miamiHomes.sf))

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
tracts18 <- 
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

tracts18 <- tracts18[miami.base,]



### Middle School
midschool <- 
  st_read("https://opendata.arcgis.com/datasets/dd2719ff6105463187197165a9c8dd5c_0.geojson") %>%
  filter(CITY == "Miami Beach" | CITY == "Miami") %>% 
  rename(midschoolID = ID)%>%
  st_transform('ESRI:102658') 

midschool <- midschool[miami.base,]

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


### Attach tract data to home prices
miamiHomes.sf <- st_join(miamiHomes.sf, tracts18, join = st_within)

### Create dummy variables for middle school using one hot encoding
miamiSchools.sf <- st_join(miamiHomes.sf, midschool, join = st_within)

sch_dmy.df <- miamiSchools.sf %>%
  st_drop_geometry() %>%
  dplyr::select(ID, NAME) 

dmy <- dummyVars(" ~ .", data = sch_dmy.df)
school.dummies<- data.frame(predict(dmy, newdata = sch_dmy.df))

# Join dummies to miamiHomes.sf
miamiHomes.sf <- inner_join(miamiHomes.sf, school.dummies, by = 'ID')

# Rename dummies
miamiHomes.sf <- miamiHomes.sf %>%
  rename(Brownsville.MS = NAMEBrownsville.Middle) %>%
  rename(CitrusGrove.MS = NAMECitrus.Grove.Middle) %>%
  rename(JosedeDiego.MS = NAMEde.Diego..Jose.Middle) %>%
  rename(GeorgiaJA.MS = NAMEJones.Ayers..Georgia.Middle) %>%
  rename(KinlochPk.MS = NAMEKinloch.Park.Middle) %>%
  rename(Madison.MS = NAMEMadison.Middle) %>%
  rename(Nautilus.MS = NAMENautilus.Middle) %>%
  rename(Shenandoah.MS = NAMEShenandoah.Middle) %>%
  rename(WestMiami.MS = NAMEWest.Miami.Middle) 


### Calculating shoreline distance
# Creating distance to shore feature
miamiHomes.sf <- miamiHomes.sf %>%
  mutate(Shore1 = nn_function(st_coordinates(st_centroid(miamiHomes.sf)),
                              st_coordinates(st_centroid(shoreline.point)),1))

miamiHomes.sf$Shore.mile <- miamiHomes.sf$Shore1/5280

### Create home Age variable and clean miamiHomes.sf for exploratory analyses #Julian added medHHInc for a future problem
miamiHomesClean.sf <- miamiHomes.sf %>%
  mutate(Age = saleYear - YearBuilt) %>%
  dplyr::select(ID, Folio, SalePrice, Property.City,
                LotSize, Bed, Bath, Stories, Pool, Fence, Patio, ActualSqFt, 
                Age, toPredict, Shore1, GEOID, MedHHInc, TotalPop, MedRent, pctWhite, pctPoverty, 
                Brownsville.MS, CitrusGrove.MS, JosedeDiego.MS, GeorgiaJA.MS, 
                KinlochPk.MS, Madison.MS, Nautilus.MS, Shenandoah.MS, WestMiami.MS, geometry) 

OLSvars <- c('ID', 'Folio', 'SalePrice','Property.City',
            'LotSize','Bed','Bath','Stories','Pool','Fence','Patio', 'ActualSqFt',
            'Age','toPredict','Shore1','GEOID', 'MedHHInc', 'TotalPop', 'MedRent','pctWhite','pctPoverty', 
            'Brownsville.MS','CitrusGrove.MS','JosedeDiego.MS','GeorgiaJA.MS', 
            'KinlochPk.MS','Madison.MS','Nautilus.MS','Shenandoah.MS','WestMiami.MS','geometry')

# --- Part 3: Exploratory Analysis ----

## Runing a Correlation Matrix to find interesting variables
miamiHomes.train <- miamiHomesClean.sf %>% 
  filter(toPredict == 0) %>%
  filter(SalePrice <= 1000000)
miamiHomes.test <- miamiHomesClean.sf %>% 
  filter(toPredict == 1)

# cor.test(miamiHomes.train$ActualSqFt, miamiHomes.train$SalePrice, method = "pearson")
# Rest of the correlations in Data section


## Methods--------
#-----Markdown The first regression we combined our feature engineering variables to see which were statistically significant
#-----Markdown the second regression includes all of the off-the-shelf features with our custom features. Model improves a lot by R2


Reg1 <- lm(miamiHomes.train$SalePrice ~ ., data = miamiHomes.train %>%
             st_drop_geometry() %>%
  dplyr::select(Shore1, MedHHInc, TotalPop, MedRent, pctWhite, pctPoverty, 
                Brownsville.MS, CitrusGrove.MS, JosedeDiego.MS, GeorgiaJA.MS, 
                KinlochPk.MS, Madison.MS, Nautilus.MS, Shenandoah.MS, WestMiami.MS,))

Reg2 <- lm(miamiHomes.train$SalePrice ~ ., data = miamiHomes.train %>%
             st_drop_geometry()    %>%
             dplyr::select(-GEOID, -ID, -toPredict))

summary(Reg1)
summary(Reg2)
summ(Reg1)

stargazer(Reg1, Reg2, title="Training Set LM Results", align=TRUE, type = "html", out = "Regression_Results.htm")



#GEOID R2 = .3, MailingZip =.4, PropertyZip =.9

# -----Part 4 - Train/Test Split -----------
#-------Markdown: Same features as are 2nd model, but similar R2 showing that its generalizable with a simple test train split

# set random seed
set.seed(31711)

# get index for training sample
inTrain <- caret::createDataPartition(
  y = miamiHomes.train$SalePrice, 
  p = .60, list = FALSE)

# split data into training and test
miami.training <- miamiHomes.train[inTrain,] 
miami.test     <- miamiHomes.train[-inTrain,]  


reg2_split <- lm(SalePrice ~ ., data = miami.training %>%
                   st_drop_geometry() %>%
             dplyr::select(-GEOID, -ID, -toPredict))

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

set.seed(73506)
#No neighborhoods R2=0.788, MAE=525,277

reg.cv2 <- 
  train(SalePrice ~ ., data = miamiHomes.train %>%
          st_drop_geometry() %>%
          dplyr::select(-ActualSqFt, -TotalPop, -ID, 
                        -MedHHInc, -toPredict), 
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

# --- Markdown: Introduction ----



# --- Markdown: Data ----
trainprice.sf <- miamiHomesClean.sf %>% 
  filter(toPredict == 0)

# Map Home Prices
allprices.sf <- full_join(trainprice.sf %>% as.data.frame(), map_preds %>% as.data.frame())%>%
  st_sf()

# Writing predictions into SalePrice variable for mapping
allprices.sf$SalePrice[allprices.sf$SalePrice == 0] <- NA   
allprices.sf$SalePrice <- ifelse(is.na(allprices.sf$SalePrice), 
                                 allprices.sf$pred, allprices.sf$SalePrice)

# All Home Prices (Predicted and Training)
ggplot() + 
  annotation_map_tile("cartolight") +
  geom_sf(data=miami.base, fill = 'transparent') + 
  geom_sf(data=allprices.sf, aes(color=SalePrice),
          show.legend = "line", size = 1) + 
  labs(title = "Home Prices",
       subtitle = "Miami & Miami Beach, FL",
       caption = "Map tiles by Carto, under CC BY 3.0. Data by OpenStreetMap, under ODbL.") +
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
  annotation_map_tile("cartolight") +
  geom_sf(data=miami.base, fill='transparent') + 
  geom_sf(data=map_preds, aes(color=pred),
          show.legend = "line", size = 1) + 
  labs(title = "Predicted Home Prices") +
  scale_color_viridis(option="C", 
                      name = "Price ($)", 
                      limits=c(10000,1000000),
                      breaks=c(0, 250000, 500000, 750000, 1000000),
                      direction = -1,
                      na.value = .97)
trn <- ggplot() + 
  annotation_map_tile("cartolight") +
  geom_sf(data=miami.base, fill='transparent') + 
  geom_sf(data=trainprice.sf, aes(color=SalePrice),
          show.legend = "line", size = 1) + 
  labs(title = "Training Set Sale Prices",
       caption = "Map tiles by Carto, under CC BY 3.0. Data by OpenStreetMap, under ODbL.") +
  scale_color_viridis(option="C", 
                      name = "Price ($)", 
                      limits=c(10000,1000000),
                      breaks=c(0, 250000, 500000, 750000, 1000000),
                      direction = -1,
                      na.value = .97) 
ggarrange(tst, trn, ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "right")

# Map of Distance to Shore
ggplot() + 
  annotation_map_tile("cartolight") +
  geom_sf(data=miami.base, fill='transparent') + 
  geom_sf(data=miamiHomes.sf, aes(colour=Shore.mile),
          show.legend = "line", size= 1) + 
  labs(title = "Distance from Shoreline",
       caption = "\U00a9 OpenStreetMap contributors") +
  scale_colour_viridis(name = "Distance (miles)")

### Middle School Areas
# Retrieve coordinates for school labels
midschool.pts <- sf::st_point_on_surface(midschool)
midschool.coords <- as.data.frame(sf::st_coordinates(midschool.pts))
midschool.coords$NAME <- midschool$NAME

# Map of Middle School Areas
ggplot() + 
  annotation_map_tile(type = 'cartolight', progress = 'none') +
  geom_sf(data=miami.base, fill='transparent') + 
  geom_sf(data=allprices.sf, aes(color=SalePrice),
          show.legend = "line", size = 1) + 
  geom_sf(data = midschool, fill = "#d3e9ff", alpha =.4, color = "#498cd3", lwd = .8)+
  labs(title = "Middle Schools & Home Prices",
       caption = "Map tiles by Carto, under CC BY 3.0. Data by OpenStreetMap, under ODbL.") +
  scale_color_viridis(option="C", 
                      name = "Price ($)", 
                      limits=c(10000,1000000),
                      breaks=c(0, 250000, 500000, 750000, 1000000),
                      direction = -1,
                      na.value = .97) 


#Summary Statistics Table
sumstat.df <- miamiHomes.train %>%
  st_drop_geometry() %>%
  mutate(Miami.dummy = ifelse(Property.City == 'Miami',1,0)) %>%
  dplyr::select(SalePrice, Miami.dummy,GEOID, 
                LotSize, Age, Stories, Bed, Bath, Pool, Fence, Patio, 
                Shore1, MedRent, pctWhite, pctPoverty, 
                Brownsville.MS, CitrusGrove.MS, JosedeDiego.MS, GeorgiaJA.MS, 
                KinlochPk.MS, Madison.MS, Nautilus.MS, Shenandoah.MS, WestMiami.MS) 

summary(sumstat.df)

stargazer(sumstat.df, title="Summary Statistics", type='html', 
          summary.stat = c('n','mean','sd','min','max'), out='Summary_Statistics.htm')

###Correlation Matrix
# Create a new subset of the training data with geometry dropped
numericVars <- miamiHomesClean.sf %>% 
  st_drop_geometry() %>%
  filter(toPredict == 0) %>%
  filter(SalePrice <= 1000000) %>%
  dplyr::select(-ID, -Folio) %>%
  select_if(is.numeric) %>% 
  na.omit() %>%
  rename(ShoreDistance = Shore1) 

# Plot correlation matrix
ggcorrplot(
  round(cor(numericVars), 1), 
  p.mat = cor_pmat(numericVars),
  colors = c("#FA7800", "white", "#25CB10"),
  type="lower",
  insig = "blank") +  
  labs(title = "Correlation across numeric variables") +
  geom_rect(aes(xmin = 0, xmax = 24.5, ymin = 0.25, ymax = 1.75),
            fill = "transparent", color = "red", size = 1)

# Home price correlation scatterplots 
a <- ggplot(trainprice.sf, aes(y=SalePrice, x = ActualSqFt)) +
  geom_point() +
  geom_smooth(method = "lm")+
  labs(title = "Price vs. SqFt (Actual)")+
  xlab(label = 'Square Footage')

b <- ggplot(trainprice.sf, aes(y=SalePrice, x = Shore1/5280)) +
  geom_point() +
  geom_smooth(method = "lm")+
  labs(title = "Price vs. Shore Distance")+
  xlab(label = 'Distance (miles)')

c <- ggplot(trainprice.sf, aes(y=SalePrice, x = MedRent)) +
  geom_point() +
  geom_smooth(method = "lm")+
  labs(title = "Price vs. Avg Rent (Tract)")+
  xlab(label='Median HHI')

d <- ggplot(trainprice.sf, aes(y=SalePrice, x = Age)) +
  geom_point() +
  geom_smooth(method = "lm")+
  labs(title = "Price vs. Home Age")+
  xlab(label='Age')

ggarrange(a,b,c,d, ncol = 2, nrow = 2)

colors <- c('#f1eef6', '#bdc9e1', '#74a9cf', '#2b8cbe', '#045a8d')
ggplot()+
  annotation_map_tile("cartolight") +
  geom_sf(data=tracts18, aes(fill = q5(pctPoverty)), lwd = 0) +
  geom_sf(data=trainprice.sf, aes(color=SalePrice), size=1) + 
  labs(title = "Poverty Rate (2018)") +
  mapTheme() + 
  theme(plot.title = element_text(size=22)) + 
  scale_fill_manual(values = colors,
                    name= "% Poverty\n(Quintile Breaks)", 
                    guide = guide_legend(reverse = TRUE)) +
  scale_color_viridis(option="C", 
                      name = "Price ($)", 
                      limits=c(10000,1000000),
                      breaks=c(0, 250000, 500000, 750000, 1000000),
                      direction = -1,
                      begin = 0,
                      end = .95,
                      na.value = .95)

ggplot() + 
  annotation_map_tile("cartolight") +
  geom_sf(data=miami.base, fill = 'transparent') + 
  geom_sf(data=allprices.sf, aes(color=SalePrice),
          show.legend = "line", size = 1) + 
  labs(title = "Home Prices",
       subtitle = "Miami & Miami Beach, FL",
       caption = "Map tiles by Carto, under CC BY 3.0. Data by OpenStreetMap, under ODbL.") +
  scale_color_viridis(option="C", 
                      name = "Price ($)", 
                      limits=c(10000,1000000),
                      breaks=c(0, 250000, 500000, 750000, 1000000),
                      direction = -1,
                      begin = 0,
                      end = .95,
                      na.value = .95) 
  


# --- Markdown: Method ----


# --- Markdown: Results ----

# Is there a spatial correlation of errors?

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
                   dplyr::select(-GEOID, -ID, 
                                 -TotalPop, -MedHHInc, -toPredict))

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

miami.test.lag <- drop_na(miami.test)
#errors
coords.test <- st_coordinates(st_centroid(miami.test.lag))
neighborList.test <- knn2nb(knearneigh(coords.test, k_nearest_neighbors))
spatialWeights.test <- nb2listw(neighborList.test, style="W")
miami.test$lagPriceError <- lag.listw(spatialWeights.test, miami.test.lag$SalePrice.AbsError)

rmggplot(miamiHomesClean.sf, aes(x=lagPrice, y=SalePrice)) +
  annotation_map_tile("cartolight") +
  geom_sf(data=miami.base, fill='transparent') + 
  geom_point(colour = "#FA7800") +
  geom_smooth(method = "lm", se = FALSE, colour = "#25CB10") +
  labs(title = "Price as a function of the spatial lag of price",
       caption = "Public Policy Analytics, Figure 6.6",
       x = "Spatial lag of price (Mean price of 5 nearest neighbors)",
       y = "Sale Price") +
  plotTheme()

ggplot(miami.test.lag, aes(x=lagPriceError, y=SalePrice)) +
  geom_point(colour = "#FA7800") +
  geom_smooth(method = "lm", se = FALSE, colour = "#25CB10") +
  labs(title = "Error as a function of the spatial lag of price",
       caption = "",
       x = "Spatial lag of errors (Mean error of 5 nearest neighbors)",
       y = "Sale Price") +
  plotTheme()


#Moran's I Test
moranTest <- moran.mc(miami.test$SalePrice.AbsError, 
                      spatialWeights.test, nsim = 999, na.action=na.exclude)

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


#Provide a map of your predicted values for where ‘toPredict’ is both 0 and 1.

PredictMap <-
  miamiHomesClean.sf %>%
  mutate(SalePrice.Predict = predict(finaltraining, miamiHomesClean.sf))

ggplot() +
  geom_sf(data = PredictMap, aes(fill = q5(SalePrice.Predict))) +
  scale_fill_manual(values = palette5,
                    labels=qBr(PredictMap,"SalePrice.Predict"),
                    name="Quintile\nBreaks") +
  mapTheme() +
  labs(title="Precited Prices of All Homes")

#This really needs a basemap or something but it works


# MAPE by Neighborhood

nhood_sum <- miami.test %>% 
  group_by(GEOID) %>%
  summarize(meanPrice = mean(SalePrice, na.rm = T),
            meanPrediction = mean(SalePrice.Predict, na.rm = T),
            meanMAE = mean(SalePrice.AbsError, na.rm = T))

nhood_sum %>% 
  st_drop_geometry %>%
  arrange(desc(meanMAE)) %>% 
  kable() %>% kable_styling()

map_preds_sum <- map_preds %>% 
  group_by(GEOID) %>% 
  summarise(meanMAE = mean(SalePrice.AbsError))

ggplot() +
  geom_sf(data = miamiHomesClean.sf %>% 
            left_join(st_drop_geometry(map_preds_sum), by = "GEOID"),
          aes(fill = q5(meanMAE))) +
  scale_fill_manual(values = palette5,
                    labels=qBr(nhood_sum,"meanMAE"),
                    name="Quintile\nBreaks") +
  mapTheme() +
  labs(title="Absolute sale price errors on the OOF set by Neighborhood")

#Scatter of MAPE & Mean Price by Neighborhood

plot(nhood_sum$meanPrice, nhood_sum$meanMAE, 
     main = "Scatter of Mean Price and Mean MAE by Census Tract", 
     xlab = "Mean Price",
     ylab = "Mean MAE")

#Testing Generalizability of Median Household Income Split

median(as.numeric(miami.test$MedHHInc), na.rm=TRUE)

miamiPoor <- miami.test %>% 
  filter(MedHHInc<= 37117)
miamiRich <- miami.test %>% 
  filter(MedHHInc> 37117)

miamiPoorLM <- 
  train(SalePrice ~ ., data = miamiPoor %>%
          st_drop_geometry() %>%
          dplyr::select(-ActualSqFt, -YearBuilt, -EffectiveYearBuilt, -TotalPop, 
                        -MedHHInc, -toPredict,
                        -Regression, -SalePrice.Predict, -SalePrice.Error, -SalePrice.AbsError, -SalePrice.APE), 
        method = "lm", 
        trControl = fitControl, 
        na.action = na.omit)

miamiPoorLM

miamiRichLM <- 
  train(SalePrice ~ ., data = miamiRich %>%
          st_drop_geometry() %>%
          dplyr::select(-ActualSqFt, -YearBuilt, -EffectiveYearBuilt, -TotalPop, 
                        -MedHHInc, -toPredict,
                        -Regression, -SalePrice.Predict, -SalePrice.Error, -SalePrice.AbsError, -SalePrice.APE), 
        method = "lm", 
        trControl = fitControl, 
        na.action = na.omit)

miamiRichLM


print(kable(miamiPoorLM$results))
print(kable(miamiRichLM$results))


rmarkdown::render("MUSA_Miami_Midterm.R")





