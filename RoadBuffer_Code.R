### Calculating distance to major roads: 1/8 mile ring buffers
### Distance from Major Roads buffer

### Read in base map
miami.base <- 
  st_read("https://opendata.arcgis.com/datasets/5ece0745e24b4617a49f2e098df8117f_0.geojson") %>%
  st_transform('ESRI:102658') %>%
  filter(NAME == "MIAMI BEACH" | NAME == "MIAMI") %>%
  st_union()

# Creating buffers for distance to major roads
roads <- 
  st_read("https://opendata.arcgis.com/datasets/8bc234275dc749329c4e242abcfc5a0f_0.geojson") %>%
  filter(CLASS == c('1','2')) %>%
  st_transform('ESRI:102658') 

miamiRds <- roads[miami.base,]

# Create unioned buffer for major roads
miamiRds.buffer <- st_union(st_buffer(miamiRds, 660)) %>%
  st_sf() %>%
  mutate(Legend = "Unioned Buffer")

miamiRds.buffer <- filter(miamiRds.buffer, Legend=="Unioned Buffer")

# Plot to check buffer
ggplot() + geom_sf(data=miami.base) +
  geom_sf(data=st_union(miamiRds.buffer), 
          color = 'red', fill = 'transparent') +
  geom_sf(data=miamiHomes.sf, aes(colour=Shore1),
          show.legend = "line", size= .5) +
  labs(title="Miami Major Roads") +
  scale_colour_viridis()

# Create 1/8 mile ring buffers
miami.rings <- multipleRingBuffer(miamiRds.buffer, 660*15, 660) %>%
  rename(road_dist = distance)

# Join home prices with ring buffer -- transform NAs to 0 
# (inner most ring not showing up as buffer, where dist is 0)
miamiHomes.rings <- st_join(miamiHomes.sf, miami.rings, join = st_within) %>%
  st_sf() 
miamiHomes.rings[c("road_dist")][is.na(miamiHomes.rings[c("road_dist")])] <- 0

miamiHomes.rings <- st_join(miamiHomes.rings, midschool, join = st_within)


