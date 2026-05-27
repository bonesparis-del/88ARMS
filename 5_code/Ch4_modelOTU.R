#### 0.0 prep ####
# building TPCindex model based on the following dataset 

# mix     : soak time & max 12 month
# bottom  : using bottom eData
# relax   : remove ARMS without full fractions    
# full    : use the full data set 
# ipdw    : use ipdw instead of idw 

## load lib
library(dplyr)      # explore data 
library(ggplot2)    # plot 
library(vegan)      # calculate distance
library(sf)         # for maping 
library(mgcv)       # gam  
library(plotly)     # fancy 3d map     
library(viridis)    # fancy color  

## set WD 
#mac 
setwd("~/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#8 modeling/1. data")
rm(list=ls())
par(mfrow = c(1, 1))

## get data
mTable.d1 <- read.csv("260114_allARMS/beta/basicTable.plus.ultra.csv",row.names = 1) # it also double as the meta data 
fTable <- read.csv("260114_allARMS/fTable/freqTableARMS.csv", row.names = 1) # 88 ARMS with 12169 OTUs
dataAll <- read.csv("260114_allARMS/eData/allArms.ipdw.idp2.csv",row.names = 1) # 12m, bottom, ipdw
mData <- read.csv("crf data/sample-metadata.byArms.csv") 
wDist.d1 <- read.csv("260114_allARMS/beta/FINAL.basicTable.csv",row.names = 1)


##### 0.1 creat phamtom mTable to make p3.fin/p3.end #####
unique(dataAll$phase)

tempTable <- mTable.d1%>% filter(phase=='p3')
mTable.d1[mTable.d1$phase=='p3',]$phase <- 'p3.fin'
tempTable[tempTable$phase=='p3',]$phase <- 'p3.beg'

mTable <- rbind(mTable.d1,tempTable) 
# it becomes 103 obs because there were 15 p3.fin in. 
# then need to change the site of p3.fin because it's not TPC and CDA any more 

mTable%>% filter(phase=='p3.beg') %>% .$siteAct # all tpc/cda 
mData$ARMSno <- paste0('HK',mData$featureid)

# align the order 
p3ARMS <- as.data.frame(mTable[mTable$phase=='p3.beg',]$ARMSno)
names(p3ARMS) <- 'ARMSno'
site.2 <- merge(p3ARMS, mData%>%filter(phase%in% c('resilience','final')), by='ARMSno')

# write it in mTable
mTable[mTable$phase=='p3.beg',]$ARMSno == site.2$ARMSno # true, it's aligned 
mTable[mTable$phase=='p3.beg',]$siteAct <- site.2$site2

# check 
tempA <- mTable[mTable$phase=='p3.beg',1:2]
tempB <- (mData%>% arrange(featureid) %>% filter(phase%in% c('resilience','final')))[,c(27,7)]
all(tempA==tempB) # all aligned 

# ok the name is right but i need to aligh the gps 
gps.temp <- mTable %>% filter(phase=='p3.beg') %>% dplyr::select('ARMSno',"siteAct")
gps.real <- mTable %>% filter(phase=='p2') %>% dplyr::select("siteAct", 'latitude','longitude')
gps.p3.beg <- merge(gps.temp,unique(gps.real),by='siteAct')

# align them again 
gps.replace <- merge(tempA,gps.p3.beg, by='ARMSno')

all(mTable[mTable$phase=='p3.beg',]$ARMSno == gps.replace$ARMSno) # all aligh 
mTable[mTable$phase=='p3.beg',]$latitude <- gps.replace$latitude
mTable[mTable$phase=='p3.beg',]$longitude <- gps.replace$longitude

##### 0.2 write eData into diversity df for modling #####

data.model <- as.data.frame(matrix(0, nrow=0, ncol=17))
for (i in unique(mTable$phase)) {
  
  # use sf to match points 
  df1_sf <- st_as_sf(dataAll %>% filter(phase==i), coords = c("lon", "lat"), crs = 4326)
  df2_sf <- st_as_sf(mTable%>% filter(phase==i), coords = c("longitude", "latitude"), crs = 4326)
  
  # Find nearest point for each site
  nearest_indices <- st_nearest_feature(df2_sf, df1_sf)
  
  # Extract chlorophyll values
  df2_with_t5 <- df2_sf %>%
    mutate(
      nearest_lon = df1_sf$lon[nearest_indices],
      nearest_lat = df1_sf$lat[nearest_indices],
      TIN = df1_sf$TIN[nearest_indices],
      chlA = df1_sf$chlA[nearest_indices],
      TEMP = df1_sf$TEMP[nearest_indices],
      OP = df1_sf$OP[nearest_indices],
      Eco = df1_sf$Eco[nearest_indices],
      distance_m = as.numeric(st_distance(df2_sf, df1_sf[nearest_indices, ], by_element = TRUE)))
  
  data.model <- rbind(data.model,df2_with_t5)
  
}


# rename p3 into p3.end because they are just half of the full set (only 6 months)
data.model.p3 <- data.model %>% filter(phase=='p3.fin')

p3.eData <- data.model %>% 
  filter(phase %in% c('p3.fin','p3.beg')) %>%
  group_by(ARMSno) %>% summarise(TIN = exp(mean(log(TIN), na.rm = TRUE)),
                                 chlA = exp(mean(log(chlA), na.rm = TRUE)),
                                 TEMP = mean(TEMP, na.rm = TRUE),
                                 OP = exp(mean(log(OP), na.rm = TRUE)),
                                 Eco = exp(mean(log(Eco), na.rm = TRUE)),
                                 n_obs = n())

data.model88 <- data.model[1:88,]
data.model88[data.model88$phase=='p3.fin',]$phase <- 'p3'

all(data.model88[data.model88$phase=='p3',]$ARMSno == p3.eData$ARMSno) # all aligned
data.model88[data.model88$phase=='p3',]$TIN <- p3.eData$TIN
data.model88[data.model88$phase=='p3',]$chlA <- p3.eData$chlA
data.model88[data.model88$phase=='p3',]$TEMP <- p3.eData$TEMP
data.model88[data.model88$phase=='p3',]$Eco <- p3.eData$Eco
data.model88[data.model88$phase=='p3',]$OP <- p3.eData$OP

# write it out 
# write.csv(data.model88, '260114_allARMS/eData/true88_12m_ipdw.csv')

##### 0.2 harsh filter remove site #####
data.model.R <- data.model88 %>% filter(fullFrac=='yes')

## what might be corelated? 
plot(Shannon~TIN, data.model.R) # this look like something too 
plot(tpcIndex~chlA, data.model.R)
plot(Evenness~TIN, data.model.R)
plot(tpcIndex~TIN, data.model.R)


#### 1.0 build models with tin and chla #### 
##### 1.1 with one term #####
# Tin
glm1_TINlg <- glm(tpcIndex ~ log(TIN),data = data.model.R,family = Gamma(link = "log"))
glm1_TIN <- glm(tpcIndex ~ TIN,data = data.model.R,family = Gamma(link = "log"))
gam1_TINlg <- gam(tpcIndex ~ s(log(TIN)),data = data.model.R, family = Gamma(link = "log"),method = "REML") 
gam1_TIN <- gam(tpcIndex ~ s(TIN),data = data.model.R, family = Gamma(link = "log"),method = "REML") 


with(summary(glm1_TINlg), 1 - deviance/null.deviance) # deviance explained 0.07687278
with(summary(glm1_TIN), 1 - deviance/null.deviance) # deviance explained 0.1027999
summary(gam1_TINlg) # deviance explained 15.7%, significant correlation 
summary(gam1_TIN) # deviance explained 10.3%, significant correlation 
AIC(gam1_TINlg,gam1_TIN,glm1_TINlg,glm1_TIN) # best gam1_TINlg 5.161166 -9.209184
BIC(gam1_TINlg,gam1_TIN,glm1_TINlg,glm1_TIN) # best glm1_TIN   3.000000 -1.157748
# overall best model is gam1_TINlg, with low AIC/BIC and high deviance explained 

# chlA 
glm1_chlAlg <- glm(tpcIndex ~ log(chlA),data = data.model.R,family = Gamma(link = "log"))
glm1_chlAlgPlus <- glm(tpcIndex ~ log(chlA+1),data = data.model.R,family = Gamma(link = "log"))
glm1_chlA <- glm(tpcIndex ~ chlA,data = data.model.R,family = Gamma(link = "log"))
gam1_chlAlg <- gam(tpcIndex ~ s(log(chlA)),data = data.model.R, family = Gamma(link = "log"),method = "REML") 
gam1_chlAlgPlus <- gam(tpcIndex ~ s(log(chlA+1)),data = data.model.R, family = Gamma(link = "log"),method = "REML") 
gam1_chlA <- gam(tpcIndex ~ s(chlA),data = data.model.R, family = Gamma(link = "log"),method = "REML") 
gam1_chlAlg.d1 <- gam(tpcIndex ~ s(log(chlA)),data = data.model.R, family = Gamma(link = "log"),method = "REML") 
summary(gam1_chlAlg.d1)


with(summary(glm1_chlAlg), 1 - deviance/null.deviance) # deviance explained 0.2678505
with(summary(glm1_chlAlgPlus), 1 - deviance/null.deviance) # deviance explained 0.2673659
with(summary(glm1_chlA), 1 - deviance/null.deviance) # deviance explained 0.2498051
summary(gam1_chlAlg) # deviance explained  35.7%, significant correlation 
summary(gam1_chlA) # deviance explained 27.4%
AIC(gam1_chlAlg,gam1_chlA,glm1_chlAlg,glm1_chlA) # best gam1_chlAlg 6.947324 -28.19357
BIC(gam1_chlAlg,gam1_chlA,glm1_chlAlg,glm1_chlA) # best glm1_chlAlg 3.000000 -18.15502
# overall best model glm1_chlAlg because simple terms 



##### 1.2 with both terms #####
# tin and chlA 
glm2 <- glm(tpcIndex~ log(TIN)+log(chlA),data = data.model.R,family = Gamma(link = "log"))
glm2.d1 <- glm(tpcIndex~ TIN+log(chlA),data = data.model.R,family = Gamma(link = "log"))
gam2 <- gam(tpcIndex ~ s(log(TIN))+s(log(chlA)),data = data.model.R, family = Gamma(link = "log"),method = "REML") 
gam2.d1 <- gam(tpcIndex ~ s(TIN)+s(log(chlA)),data = data.model.R, family = Gamma(link = "log"),method = "REML") 
gam2.d2 <- gam(tpcIndex ~ s(log(TIN))+s(chlA),data = data.model.R, family = Gamma(link = "log"),method = "REML")
gam2.d3 <- gam(tpcIndex ~ s(TIN)+s(chlA),data = data.model.R, family = Gamma(link = "log"),method = "REML")

summary(gam2)    # Deviance explained = 36.6%, chlA sig
summary(gam2.d1) # Deviance explained = 38.6%, chlA sig
summary(gam2.d2) # Deviance explained = 27.5%, chlA sig
summary(gam2.d3) # Deviance explained = 31.2%, chlA sig
summary(glm2)    # Deviance explained = 0.293904, chlA sig
summary(glm2.d1) # Deviance explained = 0.2987193, chlA sig
with(summary(glm2), 1 - deviance/null.deviance) # 0.2746118
with(summary(glm2.d1), 1 - deviance/null.deviance) # 0.276129

AIC.temp <- AIC(gam2,glm2,gam2.d1,gam2.d2,gam2.d3)

AIC.temp%>% arrange(AIC)
summary(gam2.d1)
# remark: should just go back with chla 


#### 3.0 visulazation #### 
# Create prediction grid
data.plot <- data.model.R
tin_seq <- seq(min(data.plot$TIN), max(data.plot$TIN), length = 50) # log TIN
chla_seq <- seq(min(log(data.plot$chlA)), max(log(data.plot$chlA)), length = 50)
grid <- expand.grid(TIN=tin_seq, chlA = exp(chla_seq))  # Back to original scale

data.plot


# Predict
model.best <- glm1_chlAlg
grid$tpcIndex_pred <- exp(predict(model.best, newdata = grid))

TIN_unique <- sort(unique(grid$TIN))
chlA_unique <- sort(unique(grid$chlA))
z_matrix <- matrix(grid$tpcIndex_pred, 
                   nrow = length(TIN_unique),
                   ncol = length(chlA_unique),
                   byrow = TRUE)

# Plot 3D surface
plot_ly(z = ~z_matrix, 
        x = ~TIN_unique, 
        y = ~chlA_unique,
        type = "surface", colors = viridis(100)) %>%
  layout(scene = list(
    xaxis = list(title = "TIN (mg/L)", type = "log"),
    yaxis = list(title = "Chlorophyll (μg/L)", type = "log"),
    zaxis = list(title = "tpcIndex (Species Richness)")
  ))

# plot a headmap 
plot.heat <- ggplot(grid, aes(x = TIN, y = log(chlA), fill = tpcIndex_pred)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis", name = "tpcIndex") +
  
  # Customize x-axis: show TIN values at log positions
  scale_x_continuous(
    name = "TIN (mg/L)",  # Your actual units
#    breaks = log(c(0.04, 0.6, 0.6, 0.1, 0.15,0.2, 0.3)),  # Log of the values you want to show
#    labels = c("0.04", "0.6", "0.6", "0.1", "0.15", "0.2", "0.3")  # Original units
  ) +
  
  # Customize y-axis: show chlA values at log positions
  scale_y_continuous(
    name = "Chl-a (μg/L)",  # Your actual units
    breaks = log(c(0.5, 1, 1.2, 2, 2.5, 3, 4)),  # Log of the values
    labels = c("0.5", "1", "1.2", "2", "2.5", "3", "4")  # Original units
  ) +
  
  
  labs(x = "log(TIN)", y = "log(chlA)", 
       title = "tpcIndex Response Surface",
       subtitle = "Heatmap of predicted values") +
  theme_minimal() +
  
  # Add contour line for tpcIndex = 1
  geom_contour(aes(z = tpcIndex_pred), 
               breaks = c(0.9,1,1.1), 
               color = "white", 
               linewidth = 0.5,
               linetype = "solid") +
  
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 1.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )


# saveRDS(model.best, file = "260114_allARMS/diversityModel/1termModel_glm.rds")
