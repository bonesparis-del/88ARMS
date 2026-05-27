#### 0.0 prep ####
## set WD 
#mac 
setwd("~/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#8 modeling/1. data")
rm(list=ls())
par(mfrow = c(1, 1))

## load lib
library(dplyr)
library(tidyr)
library(ggplot2)
library(Biostrings)
library(factoextra)
library(vegan)
library(stats)

## get data
mTable <- read.csv("260114_allARMS/beta/basicTable.plus.ultra.MAX.csv",row.names = 1) # it also double as the meta data 
fTable <- read.csv("260114_allARMS/fTable/freqTableARMS.csv", row.names = 1) # 88 ARMS with 12169 OTUs
ABtax <- read.csv("260114_allARMS/beta/ABtaxTable.csv", row.names = 1) 
mTable$sitePhase <- paste0(mTable$siteAct,mTable$phase)
SOM <- read.csv("260114_allARMS/hpc/taxAss/TaxAsn_shelbyOmidori.csv", row.names = 1)


#### 1.0 check model before plot #### 
summary(lm(perArthropoda~Richness.x, mTable)) # Adjusted R-squared:  -0.01127 / p-value: 0.8621
summary(lm(perAnnelida~Richness.x, mTable)) # *** Adjusted R-squared:  0.6012 / p-value: < 2.2e-16
plot(perAnnelida~Richness.x,data.plot1)
summary(lm(perBacillariophyta~Richness.x, mTable)) # ** Adjusted R-squared:  0.1509 / p-value: 0.0001088
plot(perBacillariophyta~Richness.x, data.plot1)
summary(lm(perRhodophyta~Richness.x, mTable)) # ** Adjusted R-squared:  0.1632 / p-value: 5.62e-05
plot(perRhodophyta~Richness.x, data.plot1)
summary(lm(perMollusca~Richness.x, mTable)) # Adjusted R-squared:  0.02241 / p-value: 0.08712
summary(lm(perPorifera~Richness.x, mTable)) # Adjusted R-squared:  -0.01162  / p-value: 0.9804
summary(lm(Evenness~log(Richness.x), mTable)) # *** Adjusted R-squared:  0.7155  / p-value: 2.2e-16
plot(Evenness~Richness.x, mTable)
summary(lm(InvSimpson~Richness.x, mTable)) # *** Adjusted R-squared:  0.9028  / p-value: 2.2e-16
plot(InvSimpson~Richness.x, mTable)
summary(lm(perAnnelida~Strength, mTable)) # *** Adjusted R-squared:  0.9028  / p-value: 2.2e-16
plot(perAnnelida~Strength, mTable)



summary(lm((perBacillariophyta+perRhodophyta)~Strength, mTable)) # *** Adjusted R-squared:  0.9028  / p-value: 2.2e-16
plot(perBacillariophyta+perRhodophyta~Strength, mTable)

model_ar <- lm(Arthropoda~Richness.x, mTable) # slope = 0.132255
model_an <- lm(Annelida~Richness.x, mTable) # slope = 0.062780
model_dia <- lm(Bacillariophyta~Richness.x, mTable) # slope = 0.074286
model_red <- lm(Rhodophyta~Richness.x, mTable) # slope = 0.069235

# the math was right 
mTable$Arthropoda / mTable$Assigned /mTable$perArthropoda
mTable$Annelida / mTable$Assigned / mTable$perAnnelida
mTable$Bacillariophyta / mTable$Assigned / mTable$perBacillariophyta
mTable$Rhodophyta / mTable$Assigned / mTable$perRhodophyta


model_perAr <- lm(perArthropoda~Richness.x, mTable) # slope = -2.531e-06, NOT SIG
model_perAn <- lm(perAnnelida~Richness.x, mTable) # slope = -1.249e-04, SIG
model_perDia <- lm(perBacillariophyta~Richness.x, mTable) # slope = 5.536e-05, SIG
model_perRed <- lm(perRhodophyta~Richness.x, mTable) # slope = 5.426e-05, SIG
# because arthoropoda starts high so even when it has higher increase the propostion 
# increase might not be so much sinigicant 
# like a 2 meter guy grow to 2.4 meter but less significant than 0.3 meter guy
# grow to 0.6 meter 

model_AA <- lm(Assigned~Richness.x, mTable) # it's almost 1:1 
summary(model_AA)

summary(mTable$perArthropoda)
# 0.3510/0.1413, increased by 2.484076 times
summary(mTable$Arthropoda)
# 220.00/27.00, increased by 8.148148 times

summary(mTable$perBacillariophyta)
# 0.21808/0.01676, increased by 13.01193 times 
summary(mTable$Bacillariophyta)
# 167.00/3.00, increased by 55.6667 times 

# see how the diatom increase enven with a flater slop but still have a 
# higher percentage change than the dominatiing Arthoropoda 

  
#### 2.0 plot richness ~ different percentage 
mTable$region <- 'west'
mTable[mTable$siteAct%in%c('CI', 'CLP'),]$region <- 'dTolo'
mTable[mTable$siteAct%in%c('LM', 'SW','CDA','NP'),]$region <- 'south'
mTable[mTable$siteAct%in%c('TPC', 'BI','PI','SK'),]$region <- 'east'

data.plot1 <- mTable %>% select(Community, Richness.x, perArthropoda, perAnnelida,perBacillariophyta,
                                perRhodophyta, perMollusca, perCnidaria,perPorifera, Evenness, region)
data.plot1$pp <- data.plot1$perBacillariophyta + data.plot1$perRhodophyta

data.plot2 <- pivot_longer(data.plot1, cols = c(3:9,12), names_to = "phyla")

plot1 <- ggplot(data.plot2%>% filter(phyla %in% c('perAnnelida','perRhodophyta','perBacillariophyta')), 
                aes(x=Richness.x, y=value, color=phyla, shape=phyla)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  theme_classic()



plot2 <- ggplot(data.plot2%>% filter(phyla %in% c('perBacillariophyta')), 
                aes(x=Richness.x, y=Evenness)) +
  geom_point() + 
  geom_smooth(method=lm,formula = y ~ log(x),  se=F, fullrange=TRUE)+
  theme_classic()

plot2.region <- ggplot(data.plot2%>% filter(phyla %in% c('perBacillariophyta')), 
                aes(x=Richness.x, y=Evenness)) +
  geom_point(aes(color = region)) +  # Color only for points
  geom_smooth(method=lm,formula = y ~ log(x),  se=F, fullrange=F)+
  theme_classic()


plot3 <- ggplot(data.plot2%>% filter(phyla %in% c('perAnnelida','pp')), 
                aes(x=Richness.x, y=value, color=phyla, shape=phyla)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  theme_classic()

plot4 <- ggplot(data.plot2%>% filter(phyla %in% c('perArthropoda','perAnnelida')), 
                aes(x=Richness.x, y=value, color=phyla, shape=phyla)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  theme_classic()



# some simple model 
model2.region <- aov(Richness.x~region,data.plot2%>% filter(phyla %in% c('perBacillariophyta')))
summary(model2.region) 
TukeyHSD(model2.region) 

model2.regionplus <- lm(Evenness~Richness.x*region,data.plot2%>% filter(phyla %in% c('perBacillariophyta')))
summary(model2.regionplus) 
anova(model2.regionplus)

# the same regional difference like degree/strength 
# east/south have higher richness 
  
#### 3.0 plot richness ~ different percentage with AB data #### 
# write AB data into mTable
row.names(ABtax) <- paste0('AB_',ABtax$phylum)
ABtax.t <- ABtax %>% select(-phylum) %>% t() %>% as.data.frame()
ABtax.t$AB_total <- rowSums(ABtax.t)
ABtax.t$Community <- row.names(ABtax.t)
data.plotAB <- merge(mTable, ABtax.t, by='Community')

# write.csv(data.plotAB,'260114_allARMS/beta/basicTable.plus.ultra.MAX.AB.csv' )

plotAB1 <- ggplot(data.plotAB, 
                aes(x=Richness.x, y=AB_Annelida/AB_total)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  theme_classic()


plotAB2 <- ggplot(data.plotAB, 
                aes(x=Richness.x, y=AB_Bacillariophyta/AB_total)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  theme_classic()

plotAB3 <- ggplot(data.plotAB, 
                  aes(x=Richness.x, y=AB_Rhodophyta/AB_total)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  theme_classic()


# try to plot three together 
data.plotAB$perAB_Annelida <- data.plotAB$AB_Annelida/data.plotAB$AB_total
data.plotAB$perAB_Bacillariophyta <- data.plotAB$AB_Bacillariophyta/data.plotAB$AB_total
data.plotAB$perAB_Rhodophyta <- data.plotAB$AB_Rhodophyta/data.plotAB$AB_total
data.plotAB$perAB_Arthropoda <- data.plotAB$AB_Arthropoda/data.plotAB$AB_total

data.plotAB_long <- data.plotAB %>%
  select(Richness.x, 
         perAB_Annelida,
         perAB_Bacillariophyta,
         perAB_Rhodophyta,
         perAB_Arthropoda) %>%
  pivot_longer(cols = c(perAB_Annelida, perAB_Bacillariophyta, perAB_Rhodophyta,perAB_Arthropoda ),
               names_to = "Phylum",
               values_to = "Relative_Abundance")

# Create faceted plot
plot_all <- ggplot(data.plotAB_long, 
                   aes(x = Richness.x, y = Relative_Abundance, color=Phylum)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  labs(x = "Richness", y = "Relative Abundance", 
       title = "Relative Abundance vs Richness by Phylum") +
  theme_classic() +
  theme(strip.text = element_text(face = "bold", size = 12))
#    facet_wrap(~ Phylum, ncol = 3, scales = "fixed") +  # fixed scales for comparison

# modeling richness
target <- names(data.plotAB)[23:29]

dataStat <- as.data.frame(matrix(0, nrow=length(target),ncol=3))
names(dataStat) <- c('pvalue', 'r2','estimate')
rownames(dataStat) <- target

for (i in 1:length(target)) {
  formula_str <- as.formula(paste(target[i], "~Richness.x"))
  modelTemp <- lm(formula_str, data.plotAB)
  dataStat[i,1] <- summary(modelTemp)$coefficients[2,4] 
  dataStat[i,2] <- summary(modelTemp)$adj.r.squared
  dataStat[i,3] <- summary(modelTemp)$coefficients[2,1] 
}

dataStat$adjPvalueBH <- p.adjust(dataStat$pvalue, method = "BH") 

dataStat


# modeling AB
target <- names(data.plotAB)[57:68]

dataStat.ab <- as.data.frame(matrix(0, nrow=length(target),ncol=3))
names(dataStat.ab) <- c('pvalue', 'r2','estimate')
rownames(dataStat.ab) <- target

for (i in 1:length(target)) {
formula_str <- as.formula(paste(target[i], "~Richness.x"))
modelTemp <- lm(formula_str, data.plotAB)
dataStat.ab[i,1] <- summary(modelTemp)$coefficients[2,4] 
dataStat.ab[i,2] <- summary(modelTemp)$adj.r.squared
dataStat.ab[i,3] <- summary(modelTemp)$coefficients[2,1] 
}

dataStat.ab$adjPvalueBonferroni <- p.adjust(dataStat.ab$pvalue, method = "bonferroni") 

# modeling AB
target <- names(data.plotAB)[57:68]

dataStat.PerAb <- as.data.frame(matrix(0, nrow=length(target),ncol=3))
names(dataStat.PerAb) <- c('pvalue', 'r2','estimate')
rownames(dataStat.PerAb) <- target

for (i in 1:length(target)) {
  formula_str <- as.formula(paste(target[i], "/~Richness.x"))
  modelTemp <- lm(formula_str, data.plotAB)
  dataStat.PerAb[i,1] <- summary(modelTemp)$coefficients[2,4] 
  dataStat.PerAb[i,2] <- summary(modelTemp)$adj.r.squared
  dataStat.PerAb[i,3] <- summary(modelTemp)$coefficients[2,1] 
}

dataStat.PerAb$adjPvalueBonferroni <- p.adjust(dataStat.PerAb$pvalue, method = "bonferroni") 









#### 4.0 high low diversity group ####
data.plot1$diversity <- 'high'
data.plot1[data.plot1$Richness.x<
             quantile(data.plot1$Richness.x, probs = 0.66),]$diversity <- 'mid'
data.plot1[data.plot1$Richness.x<
             quantile(data.plot1$Richness.x, probs = 0.33),]$diversity <- 'low'

colSums(data.plotAB[data.plot1$diversity=='high',57:114]) / sum(data.plot1$diversity=='high')
# AB_Annelida/AB_Arthropoda/AB_Porifera/AB_Echinodermata/AB_Bryozoa/AB_Mollusca
# AB_Rhodophyta: 7.564567e+03, AB_Bacillariophyta: 1.902533e+03 , AB_Annelida: 1.077194e+05 
colSums(data.plotAB[data.plot1$diversity=='mid',57:114]) / sum(data.plot1$diversity=='mid')
# AB_Annelida/AB_Arthropoda/AB_Porifera/AB_Echinodermata/AB_Bryozoa/AB_Mollusca
# AB_Rhodophyta: 2.804655e+03, AB_Bacillariophyta: 5.497241e+02, AB_Annelida: 7.700400e+04
colSums(data.plotAB[data.plot1$diversity=='low',57:114]) / sum(data.plot1$diversity=='low')
# AB_Annelida/AB_Arthropoda/AB_Porifera/AB_Echinodermata/AB_Bryozoa/AB_Mollusca
# AB_Rhodophyta: 1.581345e+03, AB_Bacillariophyta: 2.095862e+02, AB_Annelida: 8.311562e+04

model.An.Rich <- aov(perAnnelida~diversity, data.plot1)
model.D.Rich <- aov(perBacillariophyta~diversity, data.plot1)
model.R.Rich <- aov(perAnnelida~diversity, data.plot1)
model.Ar.Rich <- aov(perArthropoda~diversity, data.plot1)

summary(model.An.Rich)
summary(model.D.Rich)
summary(model.R.Rich)
summary(model.Ar.Rich)

TukeyHSD(model.An.Rich)
TukeyHSD(model.D.Rich)
TukeyHSD(model.R.Rich)



#### 5.0 and regional difference? #### 
model.R.region <- aov(perRhodophyta~region, mTable)
model.D.region <- aov(perBacillariophyta~region, mTable)
model.An.region <- aov(perAnnelida~region, mTable)
model.PP.region <- aov((perRhodophyta+perAnnelida)~region, mTable)


summary(model.R.region)
summary(model.D.region)
summary(model.An.region)
summary(model.PP.region)


TukeyHSD(model.R.region)
TukeyHSD(model.D.region)
TukeyHSD(model.An.region)
TukeyHSD(model.PP.region)


#### 6.0 richness/abundnace ~ richness #### 
data.plotAB_long.d1 <- data.plotAB %>%
  select(Richness.x, 
         AB_Annelida,
         AB_Bacillariophyta,
         AB_Rhodophyta,
         AB_Arthropoda) %>%
  pivot_longer(cols = c(AB_Annelida, AB_Bacillariophyta, AB_Rhodophyta,AB_Arthropoda ),
               names_to = "Phylum",
               values_to = "Abundance")

plot_all.d1 <- ggplot(data.plotAB_long.d1, 
                   aes(x = Richness.x, y = Abundance, color=Phylum)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  labs(x = "Richness", y = "Relative Abundance", 
       title = "Relative Abundance vs Richness by Phylum") +
  theme_classic() +
  theme(strip.text = element_text(face = "bold", size = 12))
#    facet_wrap(~ Phylum, ncol = 3, scales = "fixed") +  # fixed scales for comparison

model_all.d1 <- lm(Abundance~Richness.x+Phylum, data.plotAB_long.d1)
summary(model_all.d1)

anova(model_all.d1)


data.plotAB_long.d2 <- data.plotAB %>%
  select(Richness.x, 
         Annelida,
         Bacillariophyta,
         Rhodophyta,
         Arthropoda) %>%
  pivot_longer(cols = c(Annelida, Bacillariophyta, Rhodophyta,Arthropoda ),
               names_to = "Phylum",
               values_to = "Richness")

plot_all.d2 <- ggplot(data.plotAB_long.d2, 
                      aes(x = Richness.x, y = Richness, color=Phylum)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  labs(x = "Richness", y = "OTU", 
       title = "OTU vs Richness by Phylum") +
  theme_classic() +
  theme(strip.text = element_text(face = "bold", size = 12))
#    facet_wrap(~ Phylum, ncol = 3, scales = "fixed") +  # fixed scales for comparison

model_all.d1 <- lm(Abundance~Richness.x+Phylum, data.plotAB_long.d1)
summary(model_all.d1)

anova(model_all.d1)

#### 7.0 anneli is polychaet #### 
SOM %>% filter(phylum == 'Annelida') %>% 
  group_by(class) %>% summarize(obs=n()) %>% arrange(-obs)
# Polychaeta/Clitellata/Sipuncula
# Sipuncula was once classified as its own phylum but no more 

SOM %>% filter(phylum == 'Sipuncula') 

# okay dive into polychaeta then 
SOM %>% filter(class == 'Polychaeta') %>% 
  group_by(order) %>% summarize(obs=n()) %>% arrange(-obs)
# Polychaeta/Clitellata/Sipuncula
# Sipuncula was once classified as its own phylum but no more 

SOM %>% filter(class == 'Polychaeta' & identity<100) %>% 
  group_by(order) 


#### 8.0 final figure #### 
plotFIN <- ggplot(data.plot2%>% 
              filter(phyla %in% c('perAnnelida','pp','perArthropoda', 'perMollusca','perCnidaria')), 
                aes(x=Richness.x, y=value, color=phyla, shape=phyla)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  theme_classic()

summary(lm(perAnnelida~Richness.x, data.plot1)) # r2 = 0.6012, p < 2.2e-16
summary(lm(pp~Richness.x, data.plot1)) # r2 = 0.2107, p < 4.082e-06

summary(lm(perArthropoda~Richness.x, data.plot1)) 
summary(lm(perMollusca~Richness.x, data.plot1)) 
summary(lm(perCnidaria~Richness.x, data.plot1)) 


sum(data.plotAB$AB_Annelida) # 7875050
sum(data.plotAB$AB_Mollusca) # 672675

#### 9.0 regional diversity #### 
sTable <- read.csv("260114_allARMS/beta/FINAL.basicTable.csv",row.names = 1) 
sTable$region <- 'west'
sTable[sTable$siteAct%in%c('CI', 'CLP'),]$region <- 'dTolo'
sTable[sTable$siteAct%in%c('LM', 'SW','CDA','NP'),]$region <- 'south'
sTable[sTable$siteAct%in%c('TPC', 'BI','PI','SK'),]$region <- 'east'


model.aov <- aov(Richness.x ~ region, sTable%>% filter(fullFrac=='yes'))
TukeyHSD(model.aov)

model.str <- aov(Strength~ region, sTable)
TukeyHSD(model.str)

model.deg <- aov( Degree~ region, sTable)
TukeyHSD(model.deg)

model.eve <- aov(Evenness~ region, mTable)
TukeyHSD(model.eve)



# abundance also follow the OTU richness pattern so all good 

