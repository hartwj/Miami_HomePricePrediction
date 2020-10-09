
library(tidyverse)
library(sf)
library(spdep)
library(caret)
library(ckanr)
library(FNN)
library(grid)
library(gridExtra)
library(ggcorrplot)
library(ggstance)
library(jtools)     
library(broom)
# library(tufte)    #excluding for now..weird errors
library(rmarkdown)
library(kableExtra)
library(readr)


# functions
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
    summarize(pointDistance = mean(point_distance)) %>%
    arrange(as.numeric(thisPoint)) %>% 
    dplyr::select(-thisPoint) %>%
    pull()
  
  return(output)  
}



## Reading in Data

miamiHomes <- st_read("C:/Users/12156/Documents/GitHub/Miami/studentsData.geojson")
miamiHomes.sf    <- miamiHomes %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')

plot(miamiHomes.sf[,1])

nhoods <- 
  st_read("http://bostonopendata-boston.opendata.arcgis.com/datasets/3525b0ee6e6b427f9aab5d0a1d0a1a28_0.geojson") %>%
  st_transform('ESRI:102286')


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
```







