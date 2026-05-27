#### 0.0 prep #### 

## set WD 
#mac 
setwd("~/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#8 modeling/1. data")
rm(list=ls())
par(mfrow = c(1, 1))

## load lib
library(dplyr)
library(ggplot2)
library(Biostrings)
library(factoextra)
library(vegan)
library(pairwiseAdonis)
library(ggnewscale)    # to plot the gradient on x axis  
library(gridExtra)     # to plot the gradient on x axis 
library(patchwork)     # put them together  
library(ecodist)       # multiple regression on distance matries 


## get data
mTable <- read.csv("260114_allARMS/edata/true88_12m_ipdw.csv",row.names = 1) # it also double as the meta data 
fTable <- read.csv("260114_allARMS/fTable/freqTableARMS.csv", row.names = 1) # 88 ARMS with 12169 OTUs
names(mTable) <- c(names(mTable)[-1],'geometry2')
mTable$sitePhase <- paste0(mTable$siteAct,mTable$phase)
mTable <- mTable%>% arrange(ARMSno) 
wDist <- read.csv("260114_allARMS/eData/waterDist.csv",row.names = 1)


mTable$region <- 'west'
mTable[mTable$siteAct%in%c('CI', 'CLP'),]$region <- 'dTolo'
mTable[mTable$siteAct%in%c('LM', 'SW','CDA','NP'),]$region <- 'south'
mTable[mTable$siteAct%in%c('TPC', 'BI','PI','SK'),]$region <- 'east'

#### 1.0 chl x bc distance #### 
# make distance matrix
distMBCsr <-as.dist(vegdist(t(sqrt(fTable)), method = "bray")) # sqrt 

chla.m <- as.data.frame(mTable[,42])
row.names(chla.m) <- row.names(mTable)
chla.matrix  <-as.dist(vegdist(chla.m, method = "euclidean"))

# do mantel 
mantel_result.chla <- vegan::mantel(distMBCsr, chla.matrix, method = "spearman", permutations = 999)

#### 2.0 wDistance x bc distance ####
sample_to_site <- data.frame(
  Sample = mTable$ARMSno,           # e.g., "HK1", "HK2", etc.
  Site = mTable$siteAct              # e.g., "BI", "CC", "TPC", etc.
)

# 2. Create an expanded distance matrix for ALL 88 samples
n_samples <- nrow(sample_to_site)
sample_dist_matrix <- matrix(0, nrow = n_samples, ncol = n_samples)
rownames(sample_dist_matrix) <- sample_to_site$Sample
colnames(sample_dist_matrix) <- sample_to_site$Sample


# 3. Fill the matrix: for each sample pair, assign the distance between their sites
wDist.d1 <- wDist[,c(2,1, 3:11)]
names(wDist.d1)[c(1,2)] <- c('site2', 'site1')

wDist.fi <- rbind(wDist,wDist.d1)

for(i in 1:n_samples) {
  for(j in 1:n_samples) {
    site_i <- sample_to_site$Site[i]
    site_j <- sample_to_site$Site[j]
    
    # Get distance from your site-level matrix
    dist.temp <- ifelse(site_i==site_j, 0, 
                        wDist.fi %>% filter(site1==site_i & site2== site_j) %>% .$water_distance_km)
    sample_dist_matrix[i, j] <- dist.temp
  }
}

# 4. Convert to dist object (if needed for vegan functions)
w_dist <- as.dist(sample_dist_matrix)

# 5. mantel 
mantel_result.wDistance <- vegan::mantel(distMBCsr, w_dist, method = "spearman", permutations = 999)
mrm_result <- MRM(distMBCsr ~ w_dist + chla.matrix, nperm = 999)
print(mrm_result)


# 6. corelation? 
cor(as.vector(w_dist), as.vector(chla.matrix), method = "spearman") # no 


# 7. fair comparision 

# Standardize to 0-1 range
w_dist_std <- w_dist / max(w_dist)
chla_std <- chla.matrix / max(chla.matrix)

# Run MRM with standardized predictors
mrm_std <- MRM(distMBCsr ~ w_dist_std + chla_std, nperm = 999)
print(mrm_std$coef)


#### 3.0 true standard coefficient ####
# Standardize your distance matrices to mean = 0, SD = 1
distMBCsr_std <- as.dist(scale(distMBCsr, center = TRUE, scale = TRUE))
w_dist_std <- as.dist(scale(w_dist, center = TRUE, scale = TRUE))
chla_std <- as.dist(scale(log(chla.matrix+1), center = TRUE, scale = TRUE))

# Run MRM on standardized variables
mrm_std.T <- MRM(distMBCsr_std ~ w_dist_std + chla_std, nperm = 999)
mrm_std.T1 <- MRM(distMBCsr_std ~ chla_std, nperm = 999) # 
mrm_std.T2 <- MRM(distMBCsr_std ~ w_dist_std , nperm = 999)


print(mrm_std.T$coef)
print(mrm_std.T1$coef)
print(mrm_std.T2$coef)


print(chla.dist$coef)


