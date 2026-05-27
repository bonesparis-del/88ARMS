#### 0.0 prep ####

# load package 
library(ggplot2)  # For plotting
library(tidyr)
library(dplyr)
library(stringr)
library(gstat)
library(sf)
library(rnaturalearth)
library(viridis) 
library(ipdw)


# work D
setwd("~/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#8 modeling/1. data")

# clear environment 
rm(list=ls())
par(mfrow = c(1, 1))

# load data 
edata <- read.csv("251210/edata.csv",row.names = 1)
sitemeta <- read.csv("sitelocation.csv")
sitemeta <- sitemeta %>% mutate(across(where(is.character), str_trim)) # fix formate 
diversity <- read.csv("260114_allARMS/hpc/taxAss/alphaTableARMS.csv")
mData <-  read.csv("crf data/sample-metadata.byArms.csv")


#### 1.0 prep data #### 
# make a bottom water annual mean of five parameters chlA/TIN/temp/OP/Ecoli
edata.d1.1 <- edata[(edata$Depth=="Bottom Water"),] # use bottom because my ARMS are on the bottom 
edata.d1.2 <- edata %>% filter(Depth=="Surface Water") %>% filter(Station=="DM1" |Station== "DM2" |Station=="DM3"|Station== "PT2")
# because sites "DM1" "DM2" "DM3" "PT2" are shallow so only have surface water 
edata.d1 <- rbind(edata.d1.1,edata.d1.2) # 33773
edata.d1$Dates <- as.Date(edata.d1$Dates)


# fix time.in
diversity$sampledDate <- as.Date(diversity$sampledDate,"%m/%d/%y")
diversity$deployDate <- as.Date(diversity$deployDate,"%m/%d/%y")

# check if it's right 
temp.d1 <- diversity
temp.d1$days <- temp.d1$sampledDate - temp.d1$deployDate
temp.d1 %>% dplyr::select(siteAct,phase,days) %>% arrange(phase) # ok checked out 


##### 1.1 Geometric mean #####
simplePhase <- unique(diversity$phase)[-3]
i <- simplePhase[1]

## try differernt period 
time.out <- diversity[diversity$phase==i,]$sampledDate[1]
# time.in <- time.out-185 # for six month
# time.in <- diversity[diversity$phase==i,]$deployDate[1] # for soakTime/mix
 time.in <- time.out - 365 # for fix one year soak time 

## bottom/surface 
  edata.d2 <- edata.d1 %>% filter(Dates>time.in & Dates <time.out)  # bottom water
# edata.d2 <- edata.d1s %>% filter(Dates>time.in & Dates <time.out)  # surface water

edata.d3 <- as.data.frame(edata.d2 %>% 
                            group_by(Station) %>%
                            summarise(mean_chloA = exp(mean(log(Chlorophyll.a..μg.L.), na.rm = TRUE)),
                                      mean_TIN = exp(mean(log(Total.Inorganic.Nitrogen..mg.L.), na.rm = TRUE)),
                                      mean_TEMP = mean(Temperature...C., na.rm = TRUE),
                                      mean_P = exp(mean(log(Orthophosphate.Phosphorus..mg.L.), na.rm = TRUE)),
                                      mean_Eco = exp(mean(log(E..coli..cfu.100mL.), na.rm = TRUE)),
                                      phase = i,
                                      n_obs = n(),  # Number of observations per site
                                      n_month=n_distinct(month)) %>%
                            arrange(Station))

  edata.all.gmean <- edata.d3 

 
# loop all the rest 
# turn on/off surface water line
for (i in simplePhase[-1]) {
  
  time.out <- diversity[diversity$phase==i,]$sampledDate[1]
#  time.in <- diversity[diversity$phase==i,]$deployDate[1] # for "actual soak time"
  time.in <- time.out - 365 # for fix one year soak time 
#  time.in <- time.out - 185 # for fix one year soak time 
  
 edata.d2 <- edata.d1 %>% filter(Dates>time.in & Dates <time.out) # bottom water
# edata.d2 <- edata.d1s %>% filter(Dates>time.in & Dates <time.out)  # surface water
  
  edata.d3.temp <- as.data.frame(edata.d2 %>% 
                                   group_by(Station) %>%
                                   summarise(mean_chloA = exp(mean(log(Chlorophyll.a..μg.L.), na.rm = TRUE)),
                                             mean_TIN = exp(mean(log(Total.Inorganic.Nitrogen..mg.L.), na.rm = TRUE)),
                                             mean_TEMP = mean(Temperature...C., na.rm = TRUE),
                                             mean_P = exp(mean(log(Orthophosphate.Phosphorus..mg.L.), na.rm = TRUE)),
                                             mean_Eco = exp(mean(log(E..coli..cfu.100mL.), na.rm = TRUE)),
                                             phase = i,
                                             n_obs = n(),  # Number of observations per site
                                             n_month=n_distinct(month)) %>%
                                   arrange(Station))
  
  edata.all.gmean <- rbind(edata.all.gmean,edata.d3.temp) # bottom water
#  edata.all.gmean.s <- rbind(edata.all.gmean.s,edata.d3.temp) # surface water
  
} 

# 470 (94 x 5) data points, perfect  
# 4 sites using suface water and 90 sites using bottom water 

# now the tricky part, write p3 in 
# make a p3.fin and p3.beg, each with 6 months 
  

## for p3.fin
  time.out <- diversity[diversity$phase=='p3',]$sampledDate[1]
  time.in <- diversity[diversity$phase=='p2',]$sampledDate[1] 

  edata.d2 <- edata.d1 %>% filter(Dates>time.in & Dates <time.out) # bottom water
  
  edata.d3.temp <- as.data.frame(edata.d2 %>% 
                                   group_by(Station) %>%
                                   summarise(mean_chloA = exp(mean(log(Chlorophyll.a..μg.L.), na.rm = TRUE)),
                                             mean_TIN = exp(mean(log(Total.Inorganic.Nitrogen..mg.L.), na.rm = TRUE)),
                                             mean_TEMP = mean(Temperature...C., na.rm = TRUE),
                                             mean_P = exp(mean(log(Orthophosphate.Phosphorus..mg.L.), na.rm = TRUE)),
                                             mean_Eco = exp(mean(log(E..coli..cfu.100mL.), na.rm = TRUE)),
                                             phase = 'p3.fin',
                                             n_obs = n(),  # Number of observations per site
                                             n_month=n_distinct(month)) %>%
                                   arrange(Station))
  
  edata.all.gmean <- rbind(edata.all.gmean,edata.d3.temp) # bottom water


### for p3.beg
  time.out <- diversity[diversity$phase=='p2',]$sampledDate[1] 
  time.in <- time.out - 180
  
  edata.d2 <- edata.d1 %>% filter(Dates>time.in & Dates <time.out) # bottom water
  
  edata.d3.temp <- as.data.frame(edata.d2 %>% 
                                   group_by(Station) %>%
                                   summarise(mean_chloA = exp(mean(log(Chlorophyll.a..μg.L.), na.rm = TRUE)),
                                             mean_TIN = exp(mean(log(Total.Inorganic.Nitrogen..mg.L.), na.rm = TRUE)),
                                             mean_TEMP = mean(Temperature...C., na.rm = TRUE),
                                             mean_P = exp(mean(log(Orthophosphate.Phosphorus..mg.L.), na.rm = TRUE)),
                                             mean_Eco = exp(mean(log(E..coli..cfu.100mL.), na.rm = TRUE)),
                                             phase = 'p3.beg',
                                             n_obs = n(),  # Number of observations per site
                                             n_month=n_distinct(month)) %>%
                                   arrange(Station))
  
  edata.all.gmean <- rbind(edata.all.gmean,edata.d3.temp) # bottom water

  
# 658 obs for (5 + 2) * 94 


# write GPS in 
# let's take gmean 
edata.d4 <- merge(edata.all.gmean, sitemeta, by.x = "Station", by.y = "Site")   # bottom water
# edata.d4 <- merge(edata.all.gmean.s, sitemeta, by.x = "Station", by.y = "Site") # surface water 
# edata.d4 <- merge(edata.all.gmean.m, sitemeta, by.x = "Station", by.y = "Site") # middle water 
hist(log(edata.d4$mean_chloA),breaks =50)
hist(log(edata.d4$mean_TIN), breaks=50)
hist(log(edata.d4$mean_Eco), breaks=50)
hist(log(edata.d4$mean_P), breaks=50)
# should really take the geometric mean for calculation 

#### 2.0 prep map #### 
# make Hong Kong maps 
# Step 1: Load Hong Kong shapefile
hk_shape <- st_read('/Users/moicomputer/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#8 modeling/1. data/hong kong map')
hk_shape_fixed <- st_make_valid(hk_shape)
china_mainland <- ne_countries(scale = "large", country = "China", returnclass = "sf") %>%
  st_transform(4326) # i need the because shenzhen 

# step 2: make a box to include hong kong and part of shenzhen  
pred_grid <- expand.grid(
  lon = seq(113.79, 114.51, length.out = (114.51-113.79)/0.003),
  lat = seq(22.1, 22.6, length.out = (22.6-22.1)/0.003)
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

#### 3.0 write data in the map ####
# k cool, now data are ready, let's plot 
data <- edata.d4 %>%
  dplyr::select(
    Longitude = Longitude..Decimal., 
    Latitude = Latitude..Decimal., 
    chloA = mean_chloA,
    TIN = mean_TIN,
    TEMP = mean_TEMP,
    P = mean_P,
    Eco = mean_Eco,
    phase =phase
  )

# Convert plot data to spatial
data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326)



###############################  ipdw ###############################
# let's try ipdw here 


# Parameters
idp.value <- 2
paramlist <- c("TIN", "chloA", "TEMP", "P", "Eco")

# Filter data
data_sf.temp <- data_sf %>% filter(phase == unique(data_sf$phase)[1])

# Create cost raster (you'll need coastline_sf)
coastline_sf <- all_land
grid_extent <- st_bbox(pred_grid)
costras <- costrasterGen(
  xymat = as.data.frame(st_coordinates(pred_grid)),
  pols = coastline_sf,
  extent = "pnts",
  projstr = st_crs(data_sf.temp)$proj4string,
  resolution = 0.003 # ~ 300 x 300m squares, truning hk into 40,080 pixels
)
# contras is a raster of 100 x 144 pixel of hong kong wiht water = 1 and land = 10,000


# Two-step IPDW (more efficient)
path_dist_stack <- pathdistGen( # calculate distance from every cell to every sampling point 
  sf_ob = data_sf.temp, # 94 sampling points 
  costras = costras, # the area of 100 x 144 pixel of hk 
  range = 30000, # throw away anything over 25km from the sampling point 
  progressbar = TRUE
)

ipdw_results <- list()
for(var in paramlist) {
  ipdw_results[[var]] <- ipdwInterp(
    sf_ob = data_sf.temp[, var],
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

water_grid$phase <- unique(data_sf$phase)[1]
water_grid_all <- water_grid

for (i in unique(data_sf$phase)[-1]) {
  data_sf.temp <- data_sf %>% filter(phase == i)
  ipdw_results <- list()
  
  for(var in paramlist) {
    ipdw_results[[var]] <- ipdwInterp(
      sf_ob = data_sf.temp[, var],
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

water_grid$phase <- i
water_grid_all <-  rbind(water_grid_all,water_grid)
  
}


unique(water_grid_all$phase)




# Step 7: Extract results
water_data <- as.data.frame(water_grid_all)
water_data$lon <- (water_grid_all %>% st_coordinates() %>% as.data.frame() )[,1]
water_data$lat <- (water_grid_all %>% st_coordinates() %>% as.data.frame() )[,2]

water_data.d1 <- water_data[,-7]
names(water_data.d1) <- c("TIN","chlA","TEMP","OP","Eco","phase","lon","lat")
# write.csv(water_data.d1, "260114_allARMS/eData/allArms.ipdw.idp2.csv")

#### 4.0 make some map to see how the modle is like to optimize idp = 1.2/nmax = 8 ####
plot1 <- ggplot() +
  # Interpolated surface
  geom_tile(data = water_data.d1 %>% filter(phase == "p2") %>% 
              dplyr::select(chlA, lon, lat), aes(x = lon, y = lat, fill = chlA), alpha = 0.8) +
  # Hong Kong land from your shapefile
  geom_sf(data = hk_shape, fill = "lightgray", color = "black", size = 0.3) +
  # Sampling points
  geom_point(data = sitemeta, aes(x = Longitude..Decimal., y = Latitude..Decimal.), 
             color = "red", size = 0.5, alpha = 0.5) +
  # Color scale
  scale_fill_viridis(name = "Chlorophyll-A\n(μg/L)", option = "viridis")+ # viridis is the best 
  # Labels
  #  labs(title = paste0(sDate, " ~ ", eDate, " Chlorophyll-A Concentration"), #(", modelSeason, " season)"),
  #       subtitle = paste(nrow(dataReal), "Sampling Points"),
  #       x = "Longitude", y = "Latitude") +
  # Coordinate limits
  coord_sf(xlim = c(113.8, 114.5), ylim = c(22.1, 22.6)) +
  theme_minimal()


###############################




# remark in the end 
# 300m x 300m resolution to balance resolution and computational power 
# distance power = 2, seems to provide best model 


