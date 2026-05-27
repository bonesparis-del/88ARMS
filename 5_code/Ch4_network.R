#### 0.0 prep ####

# load package 
library(ggplot2)   # plot 
library(dplyr)     # explore data
library(tidyr)     # explore data
library(vegan)     # community analysis 
library(igraph)    # network analysis
library(qgraph)    # network visualization and analysis
library(tidygraph) # more plot 
library(pheatmap)  # more plot 
library(ggraph)    # more plot 

# work D
setwd("~/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#8 modeling/1. data")
rm(list=ls())
par(mfrow = c(1, 1))

# load data & model 
dataAll <- read.csv("260114_allARMS/eData/allArms.mix.Gmean.b.idp0.8.csv",row.names = 1)[,-c(3,4,5)]
diversity <- read.csv("260114_allARMS/beta/basicTable.plus.ultra.csv",row.names = 1)
fTable <- read.csv("260114_allARMS/fTable/freqTableARMS.csv",row.names = 1)

#### 1.0 Community Similarity/Dissimilarity Matrix ####
# Calculate beta diversity/dissimilarity
bray_dist <- vegdist(t(fTable), method = "bray")  # 88 × 88 distance matrix
jaccard_dist <- vegdist(t(fTable), method = "jaccard", binary = TRUE)

# Or use correlation as similarity
community_cor <- cor(fTable, method = "spearman")  # 88 × 88 correlation matrix


#### 2.0 Network Construction & Centrality Analysis ####

# Convert dissimilarity to similarity (for network edges)
bray_sim <- 1 - as.matrix(bray_dist)

# Threshold to create network (keep top 30% connections)
threshold <- quantile(bray_sim, 0.70)  # Keep values above 70th percentile
adj_matrix <- ifelse(bray_sim > threshold, bray_sim, 0)
# if two community has similarity score over 0.7 percentile then we consider they are connected

# Create network
g <- graph_from_adjacency_matrix(adj_matrix, 
                                 mode = "undirected", # Connections go BOTH ways
                                 weighted = TRUE, # account for bray-curtis score 
                                 diag = FALSE) # community dont connect to itself

str(g) 
# 88 communities 
E(g)$weight 
# 1117 connections each has it's own weight

# Calculate centrality measures
centrality <- data.frame(
  Community = V(g)$name,
  Degree = degree(g),           # Number of connections: how many communities with over 0.7 percentile similarity 
  Strength = strength(g),       # Sum of connection weights: sum all the similarity value of the connected community 
  Betweenness = betweenness(g, normalized = TRUE),  
  # Bridge communities: for any community to connect, how many shortest route passes through this one and normalize it.  
  Closeness = closeness(g, normalized = TRUE),      
  # inverse of the sum of distances to all the other vertices, higher-> not close to others.  
  Eigenvector = eigen_centrality(g)$vector,         
  # vertices with high eigenvector centralities are those which are connected to many other vertices which are
  # in turn, connected to many others
  PageRank = page_rank(g)$vector                    
  # Importance
)


# Identify source communities (top 20% strong)
source_threshold <- quantile(centrality$Strength, 0.8)
source_comms <- centrality[centrality$Strength >= source_threshold, ]
sourceARMS <- source_comms$Community
diversity %>% filter(ARMSno %in% sourceARMS)

diversity.d1 <- merge(centrality,diversity, by.x="Community", by.y="ARMSno")

#### 3.0 Spatial-Temporal Analysis  ####
# Add coordinates to centrality
centrality.plot <- merge(centrality, diversity, by.x = "Community", by.y="ARMSno")

# Plot spatial pattern of connectivity
# need to learn how to plot it on the map !!!
# Degree: how many direct friends you have 
# Strength: how close these friendships are 
plot.site <- ggplot(centrality.plot, aes(x = longitude, y = latitude, 
                       size = Strength, color = Degree)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(2, 10)) +
  theme_minimal() +
  labs(title = "Spatial Distribution of Community Connectivity")

# Temporal pattern (not what i want to look into)
centrality.plot$sampledDate <- as.Date(centrality.plot$sampledDate, "%m/%d/%y" )
plot.time <- ggplot(centrality.plot, aes(x = sampledDate, y = Strength, color = siteAct, group = siteAct)) +
  geom_line() +
  geom_point(aes(size = Betweenness)) +
  theme_minimal()

#### 4.0 Source-Sink Analysis  ####
# what is a source 
# high species richness/manu unique species/exports species to others
# Calculate potential source strength
# Communities that share many species with others
fTable.pa <- fTable
fTable.pa[fTable.pa>0] <- 1
source_strength <- data.frame(
  Community = colnames(fTable),
  Richness = colSums(fTable > 0),
  TotalAbundance = colSums(fTable),
  # Number of OTUs unique to this community
  UniqueOTUs = colSums(fTable.pa[rowSums(fTable.pa)==1,]),
  # Proportion of OTUs shared with at least 50% of communities
  SharingProp.50 = colSums(fTable.pa[rowSums(fTable.pa)>=ncol(fTable.pa)/2,])/diversity$Richness
)

# Potential source index
source_strength$SourceIndex <- 
  (scale(source_strength$SharingProp.50) + 
     scale(source_strength$Richness) + 
     scale(centrality$Strength)) / 3

diversity.d2 <- merge(source_strength,diversity.d1,by="Community")




#### 6.0 Network Modularity (Detect Sub-communities) ####
# Detect modules/clusters
set.seed(123)
modules <- cluster_louvain(g)
V(g)$module <- modules$membership
# we have four groups of cluster here 

# Module-level analysis
module_summary <- data.frame(
  # Column 1: Module number, one of the four cluster 
  Module = 1:length(modules),
  # Column 2: How many communities in this module
  N_Communities = sizes(modules),
  # Column 3: How important/central are communities in this module
  Mean_Centrality = tapply(centrality$Strength, modules$membership, mean),
  # Column 4: How spread out is this module geographically
  Spatial_Extent = tapply(diversity$siteAct, modules$membership, 
                          function(x) length(unique(x)))
)

# Identify which modules are sources
module_Betweenness <- aggregate(centrality$Betweenness, 
                               by = list(Module = modules$membership), 
                               mean)

diversity.d2$cluster <- modules$membership
diversity.d2%>% filter(cluster==1) %>% group_by(siteAct) %>% 
  summarise(n_obs = n()) %>% arrange(-n_obs)
# the good northeast ones Mean_Centrality = 10.003226 / betweenness = 0.010738416
diversity.d2%>% filter(cluster==2) %>% group_by(siteAct) %>% 
         summarise(n_obs = n()) %>% arrange(-n_obs)
# all CI and CLP, deepTolo ARMS Mean_Centrality = 5.084595 / betweenness = 0.008324741
diversity.d2%>% filter(cluster==3) %>% group_by(siteAct) %>% 
  summarise(n_obs = n()) %>% arrange(-n_obs)
# also good ones but more CDA and southeast Mean_Centrality = 10.503040 / betweenness = 0.012487112
diversity.d2%>% filter(cluster==4) %>% group_by(siteAct) %>% 
  summarise(n_obs = n()) %>% arrange(-n_obs)
# the bad ones, PC/SK/SSW and one SK Mean_Centrality = 3.788284 / betweenness = 0.008506691

# write it out 
# write.csv(diversity.d2,"260114_allARMS/beta/basicTable.plus.ultra.MAX.csv")

#### 7.0 Complete Analysis Pipeline ####
analyze_community_connectivity <- function(fTable, diversity.d2) {
  
  # 1. Calculate dissimilarity
  bray_dist <- vegdist(t(fTable), "bray")
  
  # 2. Create network
  sim_matrix <- 1 - as.matrix(bray_dist)
  adj_matrix <- ifelse(sim_matrix > quantile(sim_matrix, 0.70), sim_matrix, 0)
  g <- graph_from_adjacency_matrix(adj_matrix, weighted = TRUE)
  
  # 3. Calculate centrality
  centrality_metrics <- data.frame(
    Community = V(g)$name,  # sample/ARMS numbers 
    Degree = degree(g),     # how many neighbors (is neighbors when top 30% similar)
    Strength = strength(g), # How similar you are to all others combined 
    Betweenness = betweenness(g, normalized = TRUE), # Habitat corridor significance
    Closeness = closeness(g, normalized = TRUE),     # How quickly species can reach everyone
    Eigenvector = eigen_centrality(g)$vector         # Being friends with popular reefs
  )
                                
  # 4. Add metadata
  centrality_metrics <- merge(centrality_metrics, diversity, by.x = "Community",by.y="ARMSno")
  
  
  # 5. Calculate source score
  richness <- colSums(fTable > 0)
  sharing <- colSums(fTable.pa[rowSums(fTable.pa)>=ncol(fTable.pa)/2,])/diversity$Richness
  
  source_score <- data.frame(
    Community = colnames(fTable),
    Richness = richness,
    SharingProp = sharing,
    SourceScore = scale(richness) + scale(sharing) + scale(centrality_metrics$Strength)
  ) # i need to read and think about if this is sensible 
  
  # 6. Identify source communities
  centrality_metrics$IsSource <- source_score$SourceScore >= 
    quantile(source_score$SourceScore, 0.90)
  
  # 7. Return results
  return(list(
    Network = g,
    Centrality = centrality_metrics,
    SourceScore = source_score,
    Dissimilarity = bray_dist
  ))
}

# Run analysis
results <- analyze_community_connectivity(fTable, diversity)

# Identify top source communities
top_sources <- results$Centrality %>%
  filter(IsSource == TRUE) %>%
  arrange(desc(Strength))



#### 8.0 Visualization ####

# include region 
diversity.d2$region <- 'west'
diversity.d2[diversity.d2$siteAct%in%c('CI', 'CLP'),]$region <- 'dTolo'
diversity.d2[diversity.d2$siteAct%in%c('LM', 'SW','CDA','NP'),]$region <- 'south'
diversity.d2[diversity.d2$siteAct%in%c('TPC', 'BI','PI','SK'),]$region <- 'east'



# Network visualization
network_plot <- as_tbl_graph(results$Network) %>%
  # Join your centrality and metadata to nodes
  activate(nodes) %>%
  # Use left_join to add data (nodes are in same order as results$Centrality)
  mutate(
    Community = results$Centrality$Community,
    Site = results$Centrality$Site,
    Strength = results$Centrality$Strength,
    Betweenness = results$Centrality$Betweenness,
    Degree = centrality_degree(),
    siteAct = results$Centrality$siteAct,
    cluster = diversity.d2$cluster,
    region = diversity.d2$region,
    sourceIndex =results$SourceScore$SourceScore
  )





# Now plot
set.seed(38)  
p <- ggraph(network_plot, layout = "fr") +
  geom_edge_link(aes(alpha = weight), width = 0.5, color = "gray70") +
  geom_node_point(aes(size = Strength, color = as.factor(cluster), shape=region), 
                  alpha = 0.8) +
  scale_shape_manual(values = c(16, 17, 15, 18))+
  geom_node_label(
#  aes(label = ifelse(sourceIndex > quantile(sourceIndex, 0.9), paste0('Source_',Community,"_",siteAct),
#                     ifelse(sourceIndex < quantile(sourceIndex, 0.1), paste0('Sink_',Community,"_",siteAct),""))),
 aes(label = paste0(Community,"_",siteAct)),
    
    repel = TRUE,
    size = 3,
    fill = alpha("white", 0.7)
  ) +
  scale_size_continuous(range = c(2, 10)) +
  labs(title = "Community Network",
       subtitle = "Node size = Strength, Color = Site",
       color = "cluster") +
  theme_void()
p

