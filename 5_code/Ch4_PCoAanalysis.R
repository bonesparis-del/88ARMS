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


## get data
mTable <- read.csv("260114_allARMS/edata/true88_12m_ipdw.csv",row.names = 1) # it also double as the meta data 
fTable <- read.csv("260114_allARMS/fTable/freqTableARMS.csv", row.names = 1) # 88 ARMS with 12169 OTUs
names(mTable) <- c(names(mTable)[-1],'geometry2')
mTable$sitePhase <- paste0(mTable$siteAct,mTable$phase)
mTable <- mTable%>% arrange(ARMSno) 

#### 1.0 make data set #### 
# make distance matrix
distMeuc <- as.matrix(vegdist(t(fTable), method = "euclidean")) 
distMBC <-as.matrix(vegdist(t(fTable), method = "bray")) 
distMBCsr <-as.matrix(vegdist(t(sqrt(fTable)), method = "bray")) # sqrt 
distMJC <-as.matrix(vegdist(t(fTable), method = "jaccard")) 


mTable$region <- 'west'
mTable[mTable$siteAct%in%c('CI', 'CLP'),]$region <- 'dTolo'
mTable[mTable$siteAct%in%c('LM', 'SW','CDA','NP'),]$region <- 'south'
mTable[mTable$siteAct%in%c('TPC', 'BI','PI','SK'),]$region <- 'east'
  
###### plot all
  
  dist_matrix <- prcomp(distMBCsr) # bray crutit with squart root which is the best to visualization 
  # dist_matrix <- prcomp(distMJC[plot4ARMS,plot4ARMS]) #jaccard dist
  # dist_matrix <- prcomp(distMeuc[plot4ARMS,plot4ARMS]) # euclidean dist 
  # dist_matrix <- prcomp(distMBC[plot4ARMS,plot4ARMS])  # bray crutit 
  
  
  evalue <- fviz_eig(dist_matrix, addlabels = TRUE) # 27.6%, 16.2%
  dataPlot <- data.frame(matrix(nrow=nrow(mTable), ncol=5))
  row.names(dataPlot) <- mTable$ARMSno
  colnames(dataPlot) <- c("PC1", "PC2", "siteAct", "phase","sitePhase")
  
  dataPlot$PC1 <- dist_matrix$x[,1]
  dataPlot$PC2 <- dist_matrix$x[,2]
  dataPlot$siteAct <- mTable[,]$siteAct
  dataPlot$phase <- mTable[,]$phase
  dataPlot$sitePhase <- mTable[,]$sitePhase
  dataPlot$region <- mTable[,]$region
  dataPlot$chlA <- mTable[,]$chlA

  
  
#### 2.0 plot pca #### 
  
# try to do it with shleby's idea 
  # Create the revised plot with solid and hollow symbols
plotALL  <- ggplot(dataPlot, aes(x = PC1, y = PC2, color = region, shape = siteAct)) +
    geom_point(aes(shape = siteAct,
                   size = 80,
                   stroke=2)) +
    scale_shape_manual(values = c(1, 2, 0, 16, 17, 15, 3, 4, 8, 5, 18, 11, 12))+
    scale_color_manual(values = c(
    "west" = "#0072B2",   
    "dTolo" = "#E69F00",   
    "south" = "#D55E00", 
    "east" = "darkgreen")) +
    xlab("PC1 (27.6%)") +
    ylab("PC2 (16.2%)") +
    theme_classic()+
    theme(text = element_text(size = 17)) +
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_blank())
#          legend.position = "none")


plot(mTable$tpcIndex,mTable$Shannon)
plot(mTable$tpcIndex,mTable$InvSimpson)
plot(mTable$tpcIndex,mTable$Evenness)
plot(mTable$Richness,mTable$Shannon)
plot(mTable$Richness,mTable$InvSimpson)
plot(mTable$Richness,mTable$Evenness)

#### 3.0 Permanova #### 
row.names(mTable) <- mTable$ARMSno
set.seed(123)
model1 <- adonis2(distMBCsr~region, 
                  data=mTable,	by = "terms")
model1

set.seed(123)
model1.ph <- pairwise.adonis2(distMBCsr~region, 
                                      data=mTable,	by = "terms")
model1.ph


set.seed(123)
model2 <- adonis2(distMBCsr~region * chlA, 
                  data=mTable,	by = "terms")
model2

set.seed(123)
model2.ph <- pairwise.adonis2(distMBCsr~region*chlA, 
                              data=mTable,	by = "terms")
model2.ph


#### 4.0 lm the x axis with chla #### 
model3 <- lm(chlA~PC1,dataPlot)
summary(model3) # sig, x-axis is edata 

plot(PC1~chlA,dataPlot)

model4 <- lm(PC2~chlA,dataPlot)
summary(model4) # not sig, y-axis is region 

plot(PC2~chlA,dataPlot)


#### (VOID?) 5.0 plot with the chla gradient ####
# Create main plot
main_plot <- ggplot(dataPlot, aes(x = PC1, y = PC2, color = region, shape = siteAct)) +
  geom_point(aes(shape = siteAct), size = 4, stroke = 1) +
  scale_shape_manual(values = c(1, 2, 0, 16, 17, 15, 3, 4, 8, 5, 18, 11, 12)) +
  scale_color_manual(values = c(
    "west" = "#0072B2",   
    "dTolo" = "#E69F00",   
    "south" = "#D55E00", 
    "east" = "darkgreen"
  )) +
  xlab("PC1 (27.6%)") +
  ylab("PC2 (16.2%)") +
  theme_classic() +
  theme(
    text = element_text(size = 17),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(10, 10, 30, 10)  # Extra bottom margin for gradient bar
  )

# Create gradient bar plot
summary(dataPlot$PC1)
numbers <- seq(from = min(dataPlot$PC1), to = max(dataPlot$PC1), length.out = 100)
df.plot <- data.frame(numbers)
df.plot$chlA <- predict(model3, newdata = data.frame(PC1 = df.plot$numbers))


gradient_plot <- ggplot(df.plot, aes(x = numbers, y = 1, color = chlA)) +
  geom_point(size = 10, alpha = 1, shape=15) +
  scale_color_gradient(low = "yellow", high = "red", 
                       name = "Chl-a (μg/L)") +
  theme_void() +
  theme(
    plot.margin = margin(0, 10, 10, 10),
    legend.position = "bottom"
  ) +
  xlim(range(dataPlot$PC1)) +
  ylim(0.9, 1.1) +
  labs(color = "Chl-a")

# Combine using patchwork

combined_plot <- main_plot / gradient_plot + 
  plot_layout(heights = c(5, 1)) +
  plot_annotation(tag_levels = 'A')

print(combined_plot)
 

#### 6.0 a new plot ####
# 1. Create a matrix of your ordination scores (the axes you want to use)
ord_scores <- as.matrix(dataPlot[, c("PC1", "PC2")])

# 2. Select your environmental variable(s) as a data frame
#    Include 'chlA' and any other numeric variables you want to test.
env_data <- dataPlot[, c("chlA"), drop = FALSE]

# 3. Run envfit
set.seed(123) # for reproducibility of p-values
ef <- envfit(ord_scores ~ chlA, data = env_data, permutations = 999)
# If you have multiple variables, you can use:
# ef <- envfit(ord_scores ~ ., data = env_data, permutations = 999)

# 4. Check the results
ef

# Extract the scores (coordinates) of the fitted vectors
# The 'display = "vectors"' argument retrieves the arrow head coordinates [citation:8]
vector_coords <- as.data.frame(scores(ef, display = "vectors"))

# Add variable names as a column
vector_coords$variable <- rownames(vector_coords)

# The arrows from envfit are unit vectors. To make them visually useful,
# we need to scale them. A common approach is to multiply them by the
# square root of their R-squared (r2), which makes stronger predictors
# have longer arrows [citation:6][citation:8].
r2_data <- ef$vectors$r  # Extract r-squared values
scaling_factor <- sqrt(r2_data)

# Apply the scaling
vector_coords$PC1_scaled <- vector_coords$PC1 * scaling_factor
vector_coords$PC2_scaled <- vector_coords$PC2 * scaling_factor


# plot 
# Define a multiplier to make arrows a good length for the plot.
# You may need to adjust this number (e.g., 0.5, 1, 2) based on your data range.
arrow_multiplier <- 1

plotALL.d1  <- ggplot()+
  #dataPlot, aes(x = PC1, y = PC2, color = region, shape = siteAct)) +
  geom_point(data = dataPlot,
             aes(x = PC1, y = PC2, color = region,shape = siteAct),
             size = 3, alpha = 0.8) +
  scale_shape_manual(values = c(1, 2, 0, 16, 17, 15, 3, 4, 8, 5, 18, 11, 12))+
  scale_color_manual(values = c(
    "west" = "#0072B2",   
    "dTolo" = "#E69F00",   
    "south" = "#D55E00", 
    "east" = "darkgreen")) +
  xlab("PC1 (27.6%)") +
  ylab("PC2 (16.2%)") +
  # make the arrow 
  geom_segment(data = vector_coords,
               aes(x = 0, y = 0,
                   xend = PC1_scaled * arrow_multiplier,
                   yend = PC2_scaled * arrow_multiplier),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 1) +
  # Add labels for the variables at the arrowheads
  geom_text(data = vector_coords,
            aes(x = PC1_scaled * arrow_multiplier * 1.1,
                y = PC2_scaled * arrow_multiplier * 1.1,
                label = variable),
            color = "black", size = 4, fontface = "bold") +
  theme(text = element_text(size = 17)) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank())+
  theme_classic() 
#          legend.position = "none")



ggplot() +
  # Plot the community points
  geom_point(data = dataPlot,
             aes(x = PC1, y = PC2, color = region),
             size = 3, alpha = 0.8) +
  # Define colors for regions (customize as needed)
  scale_color_manual(values = c("east" = "darkgreen", "south" = "#D55E00",
                                "west" = "#0072B2", "dTolo" = "#E69F00")) +
  
  # Add the fitted environmental vectors as arrows
  geom_segment(data = vector_coords,
               aes(x = 0, y = 0,
                   xend = PC1_scaled * arrow_multiplier,
                   yend = PC2_scaled * arrow_multiplier),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 1) +
  
  # Add labels for the variables at the arrowheads
  geom_text(data = vector_coords,
            aes(x = PC1_scaled * arrow_multiplier * 1.1,
                y = PC2_scaled * arrow_multiplier * 1.1,
                label = variable),
            color = "black", size = 4, fontface = "bold") +
  
  # Labels with variance explained (from your earlier calculations)
  labs(x = "PC1 (27.6%)", y = "PC2 (16.2%)",
       color = "Region",
       title = "PCoA with Environmental Vectors (chl-a)") +
  theme_classic() +
  theme(text = element_text(size = 12),
        legend.position = "right")
