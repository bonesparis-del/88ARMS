# load package 
library(ggplot2)
library(sf)
library(dplyr)
library(automap)
library(gstat)
library(rnaturalearth)
library(viridis)
library(raster)
library(tidyverse)
library(mgcv)
library(vegan)
library(factoextra)
library(car)
library(cowplot)
library(ipdw)


# work D
setwd("~/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#8 modeling/1. data")
rm(list=ls())

# load data & model 
dataAll <- read.csv("260114_allARMS/eData/allArms.ipdw.idp2.csv",row.names = 1)[,-c(3,4,5)]
loaded_model <- readRDS("260114_allARMS/diversityModel/1termModel_glm.rds")
edata <- read.csv("251210/edata.csv",row.names = 1)
sitemeta <- read.csv("sitelocation.csv")
sitemeta <- sitemeta %>% mutate(across(where(is.character), str_trim))# fix formate 
diversity <- read.csv("260114_allARMS/beta/basicTable.plus.ultra.MAX.csv", row.names =1 )


#### 1.0 make hong kong map ####
# let's do mapping 
# Step 1: Load Hong Kong shapefile
hk_shape <- st_read('/Users/moicomputer/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#8 modeling/1. data/hong kong map')
hk_shape_fixed <- st_make_valid(hk_shape)
china_mainland <- ne_countries(scale = "large", country = "China", returnclass = "sf") %>%
  st_transform(4326) # i need the because shenzhen 

# Check
str(hk_shape_fixed)
plot(hk_shape_fixed$geometry)  # Quick visual check
plot(china_mainland$geometry)

# step 2: make a box to include hong kong and part of shenzhen  
pred_grid <- expand.grid(
  lon = seq(113.8, 114.5, length.out = 200),
  lat = seq(22.1, 22.6, length.out = 200)
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Combine Hong Kong and China mainland (shenzhen)
china_cropped <- st_crop(china_mainland, pred_grid)
all_land <- st_union(hk_shape_fixed, china_cropped)
plot(all_land$geometry)

# Step 3: Ensure coordinate systems match
# If your shapefile has a different CRS, transform it to WGS84 (EPSG:4326)
if(st_crs(hk_shape_fixed)!= st_crs(4326)) {
  hk_shape <- st_transform(hk_shape_fixed, crs = 4326)
  print("Transformed shapefile to WGS84")
}

# Step 5: Remove land points using your shapefile
water_grid <- pred_grid[!lengths(st_intersects(pred_grid, all_land)) > 0, ]
print(paste("Water points after masking:", nrow(water_grid)))
plot(water_grid$geometry)  # Quick visual check: worked out well


#### 2.0 make hong kong diversity predition ####

##### 2.2 with annual edata #####
start.day <- "2017-05-01"  
end.day <- "2022-12-31"

# make a bottom water annual mean of five parameters chlA/TIN
edata.d1.1 <- edata[(edata$Depth=="Bottom Water"),] # use bottom because my ARMS are on the bottom 
edata.d1.2 <- edata %>% filter(Depth=="Surface Water") %>% filter(Station=="DM1" |Station== "DM2" |Station=="DM3"|Station== "PT2")
# because sites "DM1" "DM2" "DM3" "PT2" are shallow so only have surface water 
edata.d1 <- rbind(edata.d1.1,edata.d1.2) 
edata.d1$Dates <- as.Date(edata.d1$Dates)

edata.d2 <- edata.d1 %>% filter(Dates>=start.day & Dates <= end.day)

edata.d3 <- as.data.frame(edata.d2 %>% 
                            group_by(Station) %>%
                            summarise(mean_chloA = exp(mean(log(Chlorophyll.a..μg.L.), na.rm = TRUE)),
                                      mean_TIN = exp(mean(log(Total.Inorganic.Nitrogen..mg.L.), na.rm = TRUE)),
                                      n_obs = n(),  # Number of observations per site
                                      n_month=n_distinct(month)) %>%
                            arrange(Station))


edata.d4 <- merge(edata.d3, sitemeta, by.x = "Station", by.y = "Site")

# turn 94 data into all map data (ipdw)
data <- edata.d4 %>%
  dplyr::select(
    Longitude = Longitude..Decimal., 
    Latitude = Latitude..Decimal., 
    chloA = mean_chloA,
    TIN = mean_TIN
  )


data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326)

ipdw_results <- list()

idp.value <- 2

# make the first loop 

# Create cost raster (you'll need coastline_sf)
coastline_sf <- all_land
grid_extent <- st_bbox(pred_grid)
costras <- costrasterGen(
  xymat = as.data.frame(st_coordinates(pred_grid)),
  pols = coastline_sf,
  extent = "pnts",
  projstr = st_crs(data_sf)$proj4string,
  resolution = 0.003 # ~ 300 x 300m squares, truning hk into 40,080 pixels
)
# contras is a raster of 100 x 144 pixel of hong kong wiht water = 1 and land = 10,000


# Two-step IPDW (more efficient)
path_dist_stack <- pathdistGen( # calculate distance from every cell to every sampling point 
  sf_ob = data_sf, # 94 sampling points 
  costras = costras, # the area of 100 x 144 pixel of hk 
  range = 30000, # throw away anything over 25km from the sampling point 
  progressbar = TRUE
)

ipdw_results <- list()
for(var in c('chloA','TIN')) {
  ipdw_results[[var]] <- ipdwInterp(
    sf_ob = data_sf[, var],
    rstack = path_dist_stack,
    paramlist = var,
    overlapped = FALSE,
    dist_power = 2
  )
  names(ipdw_results[[var]]) <- paste0(var, "_pred")
}

# Extract to grid
grid_coords <- st_coordinates(water_grid)
grid_df <- data.frame(
  x = grid_coords[, 1],
  y = grid_coords[, 2]
)

# Now extract
for(var in names(ipdw_results)) {
  water_grid[[paste0(var, "_ipdw")]] <- raster::extract(
    ipdw_results[[var]], 
    grid_df[, c("x", "y")]
  )
}



# Combine results if needed
water_data.d1 <- as.data.frame(water_grid)
water_data.d1$lon <- st_coordinates(water_grid$geometry)[,1]
water_data.d1$lat <- st_coordinates(water_grid$geometry)[,2]


water_data <- water_data.d1[,c(2,4,5)]
names(water_data)[1] <- c("chlA")

water_data$TPCindex_pred <- exp(predict(loaded_model, water_data))


#### 3/4 make plot dataset ####
Strength.mean <- as.data.frame(diversity %>% group_by(siteAct) %>% summarise(Strength.mean=mean(Strength),
                                                                             Degree.mean=mean(Degree),
                                                                             n_obs = n())) %>% arrange(-Strength.mean)

diversity.plot <- merge(Strength.mean,diversity, by.x='siteAct',by.y='siteAct')

#### 4.0 plot them together #### 
combined_plot <- ggplot() +
  # Layer 1: Interpolated surface (geom_tile)
  geom_tile(data = water_data, 
            aes(x = lon, y = lat, fill = TPCindex_pred), 
            alpha = 0.8) +
  
  # Layer 2: Hong Kong land
  geom_sf(data = hk_shape, 
          fill = "lightgray", color = "black", size = 0.3) +
  
  # Layer 3: Diversity points (from your second plot)
  geom_point(data = diversity.plot, 
             aes(x = longitude, y = latitude, 
                 size = Strength.mean, color = Degree.mean),
             alpha = 1) +
  
  # Color scales
  scale_fill_viridis(name = "TPC Index", option = "G") + # D/G is ok
  scale_size_continuous(name = "Strength", range = c(2, 7)) +
  scale_color_gradient(name = "Degree", low = "darkblue", high = "#f23e16") +
  
  # Labels
  labs(title = paste0("Model Biodiversity (", start.day, "~", end.day, ")"),
       subtitle = paste(nrow(diversity), "Community Points"),
       x = "Longitude", y = "Latitude") +
  
  # Coordinate limits
  coord_sf(xlim = c(113.8, 114.5), ylim = c(22.1, 22.6)) +
  
  theme_minimal()


#### 5.0 calculate source index #### 
diversity$region <- 'west'
diversity[diversity$siteAct%in%c('CI', 'CLP'),]$region <- 'dTolo'
diversity[diversity$siteAct%in%c('LM', 'SW','CDA'),]$region <- 'south'
diversity[diversity$siteAct%in%c('TPC', 'BI','PI','NP'),]$region <- 'east'


model1 <- aov(Strength~region, diversity)
summary(model1)
TukeyHSD(model1)

model2 <- aov(Degree~region, diversity)
summary(model2)
TukeyHSD(model2)

p1 <- ggplot(diversity, aes(x=region, y=Strength)) + 
  geom_boxplot()

p2 <- ggplot(diversity, aes(x=region, y=Degree)) + 
  geom_boxplot()

