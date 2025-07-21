# Date: 02/20/2025 
# Goal: Process gcor results to identify traits for MVMR
# Using Gcor results of PheWAS-Ttest snp-trait outliers
# Author: KJS

#Set working direcotry
setwd("/pl/active/colelab/users/kjames/refinedMR/")

#install.packages("ggraph")
#install.packages("tidygraph")  # Needed for network visualization
#install.packages("igraph") # for louvain cluster
library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(ggplot2)

# May use this in the future for programatically completing these tasks
# File names
#files <- c("interim_data/ldsc/PheWAS_Ttest_LDSC_input/oilyfish_GCST90239664_TG_Graham_descending_TTest_filenames.txt",
#"interim_data/ldsc/PheWAS_Ttest_LDSC_input/oilyfish_GCST90132314_CVD_Aragam_descending_TTest_filenames.txt",
#"interim_data/ldsc/PheWAS_Ttest_LDSC_input/oilyfish_GCST006901_Height_Yengo_2018_descending_TTest_filenames.txt",
#"interim_data/ldsc/PheWAS_Ttest_LDSC_input/whiteXwholegrain_GCST90239658_LDL_Graham_descending_TTest_filenames.txt",
#"interim_data/ldsc/PheWAS_Ttest_LDSC_input/whiteXwholegrain_GCST90132314_CVD_Aragam_descending_TTest_filenames.txt",
#"interim_data/ldsc/PheWAS_Ttest_LDSC_input/whiteXwholegrain_GCST006901_Height_Yengo_2018_descending_TTest_filenames.txt",
#"interim_data/ldsc/PheWAS_Ttest_LDSC_input/alcohol_GCST90013405_ALT_Pazoki_descending_TTest_filenames.txt",
#"interim_data/ldsc/PheWAS_Ttest_LDSC_input/alcohol_GCST90319877_cir_Ghouse_descending_TTest_filenames.txt",
#"interim_data/ldsc/PheWAS_Ttest_LDSC_input/alcohol_GCST006901_Height_Yengo_2018_descending_TTest_filenames.txt") 

# File with gcors from traits identified in PheWAS-Ttest method
df <- read.table("/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/processed_results/gcor_PheWASttest_results_all.txt")

# Header
gcor_header <- c("p1", "p2", "rg", "se", "z", "p", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se", "gcov_int", "gcov_int_se")

# Assign header to df
colnames(df) <- gcor_header

# How many times does "1448_4" occur in the p1 column?
count_1448_4 <- df %>% filter(p1 == "1448_4.ldsc.imputed_v3.both_sexes.tsv.bgz") %>% tally() #121
count_1448_4_df <- df %>% filter(p1 == "1448_4.ldsc.imputed_v3.both_sexes.tsv.bgz")
count_1687 <- df %>% filter(p1=="1687.ldsc.imputed_v3.both_sexes.tsv.bgz") %>% tally() #120, expected

# Ensure all pairs are in the upper triangle (P1 < P2)
df_upper <- df %>%
  mutate(P1_new = pmin(p1, p2),  # Assign the smaller value to P1
         P2_new = pmax(p1, p2)) %>%  # Assign the larger value to P2
  select(P1_new, P2_new, rg) %>%
  rename(p1 = P1_new, p2 = P2_new)

# Is it whole?
ggplot(df_upper, aes(x = p1, y = p2, fill = rg)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(df_upper$rg)) +
  labs(title = "Upper Triangle Heatmap", x = "Trait 1", y = "Trait 2") +
  theme_minimal()
# THANK GOD! Complete upper triangle.

# Replace file names with trait descriptions
key <- read.csv("/pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/metadata/UKBB GWAS Imputed v3 - File Manifest Release 20180731 - Description Lookup.csv",header = FALSE)

# Rename dfupper df to keep code simpler
df <- df_upper

# Isolate the ID portion from file name
df$p1id <- sub("\\..*","", df$p1)
df$p2id <- sub("\\..*","", df$p2)

# Merge
dfm <- merge(df, key, by.x = "p1id", by.y = "V1", all.x = TRUE, all.y = FALSE)
dfm <- merge(dfm, key, by.x = "p2id", by.y = "V1", all.x = TRUE, all.y = FALSE)

# Keep unique
dfmu <- unique(dfm)

# Double check that everything is accounted for
# Make sure top triangle is full
ggplot(dfmu, aes(x = p1, y = p2, fill = rg)) +
  geom_tile() +  # Heatmap tiles
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) 

# Rename 
names(dfmu)[which(names(dfmu)=="V2.x")] <- "p1_description"
names(dfmu)[which(names(dfmu)=="V2.y")] <- "p2_description"

# Move descriptions to first two positions
dfmu2 <- dfmu %>% select(p1_description, p2_description, everything())

# Ensure all pairs are in the upper triangle (P1 < P2)
# Have to redo with the description columns
dfmu2 <- dfmu %>%
  mutate(P1_desc_new = pmin(p1_description, p2_description),  # Assign the smaller value to P1
         P2_desc_new = pmax(p1_description, p2_description)) %>%  # Assign the larger value to P2
  select(P1_desc_new, P2_desc_new, p2id, p1id, p1, p2, rg) %>%
  rename(p1_description = P1_desc_new, p2_description = P2_desc_new)

ggplot(dfmu2, aes(x = p1_description, y = p2_description, fill = rg)) +
  geom_tile() +  # Heatmap tiles
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Color scale
  theme_minimal() +
  labs(title = "Genetic Correlations",
       x = "P1",
       y = "P2",
       fill = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Develop weights for Louvain clusters
# Test weights
weight <- read.csv("/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_FilePaths_oilyfish_GCST90239664_TG_Graham.csv")
weight_count <- weight %>%
  distinct(rmd_SNPnames, rmd_traitnames) %>%  # Remove duplicate SNP-Trait pairs
  count(rmd_traitnames, name = "unique_SNPS_per_trait")  # Count unique SNPs per trait

# How many SNPs are unique?
unique_snps <- length(unique(weight$rmd_SNPnames))

# For the weights, (# unique snps per trait)/(# SNPs removed)
weight_count$trait_weight <- weight_count$unique_SNPS_per_trait/unique_snps

# Isolate the ID portion from file name
weight_count$rmd_trait_id <- sub("\\..*","", weight_count$rmd_traitnames)

#############
# Get traits specific to a diet-outcome pair

# File with traits specific to an MR analysis
traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/oilyfish_GCST90239664_TG_Graham_descending_TTest_filenames.txt", header = FALSE)
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/oilyfish_GCST90132314_CVD_Aragam_descending_TTest_filenames.txt", header = FALSE)
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/oilyfish_GCST006901_Height_Yengo_2018_descending_TTest_filenames.txt", header = FALSE)
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/whiteXwholegrain_GCST90239658_LDL_Graham_descending_TTest_filenames.txt", header = FALSE)
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/whiteXwholegrain_GCST90132314_CVD_Aragam_descending_TTest_filenames.txt", header = FALSE)
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/whiteXwholegrain_GCST006901_Height_Yengo_2018_descending_TTest_filenames.txt", header = FALSE)
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/alcohol_GCST90013405_ALT_Pazoki_descending_TTest_filenames.txt", header = FALSE)
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/alcohol_GCST90319877_cir_Ghouse_descending_TTest_filenames.txt", header = FALSE)
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/alcohol_GCST006901_Height_Yengo_2018_descending_TTest_filenames.txt", header = FALSE)

# Get unique traits
traitsu <- unique(traits)
traitsu$V1 <- gsub("/pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/LDSC/", "", traitsu$V1)

# Filter df to traitsu
df2 <- dfmu2[dfmu2$p1 %in% traitsu$V1 & dfmu2$p2 %in% traitsu$V1,] #use AND not OR
df <- df2

# Merge weights with df
#df <- merge(df, weight_count, by.x = "p1id", by.y = "rmd_trait_id")

# Move descriptions to first two positions [order matters for clustering!]
#df <- df %>% select(p1_description, p2_description, everything())


# Ensure all pairs are in the upper triangle (P1 < P2)
# Have to redo with the description columns
df <- df %>%
  mutate(P1_desc_new = pmin(p1_description, p2_description),  # Assign the smaller value to P1
         P2_desc_new = pmax(p1_description, p2_description)) %>%  # Assign the larger value to P2
  select(P1_desc_new, P2_desc_new, p2id, p1id, p1, p2, rg) %>%
  rename(p1_description = P1_desc_new, p2_description = P2_desc_new)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot to get the vibe
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Heat map
#df_pfiltered <- df[df$p<0.05,]
ggplot(df, aes(x = p1_description, y = p2_description, fill = rg)) +
  geom_tile() +  # Heatmap tiles
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Color scale
  theme_minimal() +
  labs(title = "Genetic Correlations of Significant Pairs",
       x = "P1",
       y = "P2",
       fill = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Try louvain algorithm to match like things together
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter for strong associations (rg value > 0.6)
df_rgfiltered <- subset(df, rg > 0.4) # high so that we can find things that are the "same" to avoid colinearity problems
#df_rgfiltered <- subset(dfmu, rg > 0.4)
summary(df_rgfiltered) #note, by subsetting to rg>0.6 all p-values are higly sig

# Create a graph object
g <- igraph::graph_from_data_frame(df_rgfiltered, directed = FALSE)

# Apply a community detection algorithm (Louvain method)
clusters_louvain <- igraph::cluster_louvain(g, weights = E(g)$trait_weight)
#clusters_louvain <- igraph::cluster_louvain(g) # weights = E(g)$trait_weight)
#clusters_leiden <- igraph::cluster_leiden(g)

# NOTE
# Being lazy... assign which model you want to play with
clusters <- clusters_louvain
#clusters <- clusters_leiden


# Get cluster memberships df
membership_df <- data.frame(Node = names(membership(clusters)), Cluster = as.vector(membership(clusters)))

# Print groups
print(membership_df)

# Plot, simple
plot(clusters,g)

# Plot, pretty
# Convert igraph object to a tidygraph for ggraph
g_tidy <- as_tbl_graph(g)

# Add cluster membership as a node attribute
g_tidy <- g_tidy %>%
  mutate(cluster = as.factor(membership(clusters)))#, degree = degree(g)) 
  #mutate(degree = degree(g)) # degree centrality as a node attribute

# Plot with ggraph
ggraph(g_tidy, layout = "fr") +  # Force-directed layout
  #geom_edge_link(aes(width = trait_weight), alpha = 0.5, color = "gray50") +  # Edge thickness based on weight
  geom_edge_link(aes(alpha = 0.5), color = "gray50") +  # Edges (relationships) # no weight
  #geom_node_point(aes(size = degree, color = degree), show.legend = TRUE) + # Size will be based on degree centrality
  geom_node_point(aes(color = cluster), size = 5) +  # Nodes colored by cluster
  geom_node_text(aes(label = name), repel = TRUE, size = 5) +  # Labels with repulsion
  scale_color_viridis_d(option = "plasma") +  # Color clusters with Viridis
  theme_void() +  # Clean background
  labs(title = "Network Graph of Significant Associations", color = "Cluster")


# What gets removed? Not correlated with anything?
df_rg_no <- subset(df, rg <= 0.4)

# Combine df2$p1 and df2$p2 to create a reference list
unq_traits_df <- unique(c(df_rgfiltered$p1_description, df_rgfiltered$p2_description))

# Find traits in df$p1 that are NOT in df2$p1 or df2$p2
missing_traits <- setdiff(df_rg_no$p1_description, unq_traits_df)
missing_traits

############
# Measures the number of direct connections a node has.
# High-degree nodes are often hubs in the network.
# degree centrality
degree_centrality <- degree(g)
# Find the top most connected nodes
top_degree_nodes <- sort(degree_centrality, decreasing = TRUE)[1:10]
print(top_degree_nodes)

# Measures how often a node lies on the shortest path between other nodes.
# High betweenness nodes are important for communication flow.
# Compute betweenness centrality (think bridges)
betweenness_centrality <- betweenness(g, normalized = TRUE)
# Get the top nodes with highest betweenness
top_betweenness_nodes <- sort(betweenness_centrality, decreasing = TRUE)[1:10]
print(top_betweenness_nodes)

# Nodes are important if they connect to other important nodes.
# For finding nodes that influence other highly connected nodes
# Compute eigenvector centrality
eigenvector_centrality <- eigen_centrality(g)$vector
# Find the most influential nodes
top_eigenvector_nodes <- sort(eigenvector_centrality, decreasing = TRUE)[1:10]
print(top_eigenvector_nodes)

# Most "authoratative/influential" nodes
# Compute PageRank
pagerank_scores <- page_rank(g)$vector
# Find the highest-ranked nodes
top_pagerank_nodes <- sort(pagerank_scores, decreasing = TRUE)[1:10]
print(top_pagerank_nodes)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Try hierarchical clustering
library(tidyr)

# Convert long format to wide format
#dfsub <- select(df, c(p1id, p2id, rg))
df_wide <- tidyr::pivot_wider(dfsub, names_from = p2id, values_from = rg, values_fill = 0)

dfsub <- select(df, c(p1, p2, rg))
df_wide <- tidyr::pivot_wider(dfsub, names_from = p2, values_from = rg)
df_wide <- as.data.frame(df_wide)

#dfsub <- select(df, c(p1_description, p2_description, rg))
#df_wide <- tidyr::pivot_wider(dfsub, names_from = p2_description, values_from = rg, values_fill = 0)
#df_wide <- as.data.frame(df_wide)
# Set row names as SNP1 and remove SNP1 column
#rownames(df_wide) <- df_wide$p1_description
#df_wide$p1_description <- NULL

# Set row names as SNP1 and remove SNP1 column
rownames(df_wide) <- df_wide$p1id
df_wide$p1id <- NULL

# Add usefull names
df_wide2 <- merge(df_wide, key, by.x = 0, by.y = "V1")
rownames(df_wide2) <- df_wide2$V2

# Convert correlation to a distance metric
dist_matrix <- as.dist(1 - df_wide)  # Distance = 1 - correlation

# Print distance matrix
print(dist_matrix)

# Perform hierarchical clustering
hc <- hclust(dist_matrix, method = "ward.D")  # Ward's method minimizes variance

# Plot dendrogram
plot(hc, main = "Hierarchical Clustering of SNPs", xlab = "", sub = "")

# Cut into 3 clusters
clusters <- cutree(hc, k = 3)

# Print SNP cluster assignments
print(clusters)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Try k-means clustering
library(cluster)   # For silhouette score
library(factoextra) # For visualization
library(tidyr)

# Convert to wide format
wide_df <- reshape2::dcast(df, p1 ~ p2, value.var = "rg")

# Ensure row names and column names match
rownames(wide_df) <- wide_df$p1
wide_df$p1 <- NULL # get rid of p1

# Find row names that are not in column names
missing_columns <- setdiff(rownames(wide_df), colnames(wide_df))
missing_rows <- setdiff(colnames(wide_df), rownames(wide_df))

# Add the missing column "1448..."
if (!("1448_4.ldsc.imputed_v3.both_sexes.tsv.bgz" %in% colnames(wide_df))) {
  wide_df$"1448_4.ldsc.imputed_v3.both_sexes.tsv.bgz" <- NA
}

new_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(wide_df)))
rownames(new_row) <- "78_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz"
colnames(new_row) <- colnames(wide_df)
wide_df <- rbind(wide_df, new_row)

# Add the missing row "78..."
if (!("78_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz" %in% rownames(wide_df))) {
  new_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(wide_df)))
  colnames(new_row) <- colnames(wide_df)
  rownames(new_row) <- "78_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz"
  wide_df <- rbind(wide_df, new_row)
}

# Reorder columns to match row order
wide_df2 <- wide_df[, intersect(rownames(wide_df), colnames(wide_df))]

# Back up
wide_matrix <- wide_df2

# Fill lower triangle with upper triangle values
wide_matrix[lower.tri(wide_matrix)] <- t(wide_matrix)[lower.tri(wide_matrix)]

# Replace diagonal with 1 (self-correlation)
diag(wide_matrix) <- 1

# Print result
print(wide_matrix)

# Convert correlation matrix to distance metric (1 - correlation)
# euclidean by default
distance_matrix <- as.dist(1 - wide_matrix)

# wss=within cluster sums of squares
fviz_nbclust(as.matrix(distance_matrix), kmeans, method = "wss") +
  ggtitle("Elbow Method for Optimal K")

# silhouette=average sihouette width; high is good 
fviz_nbclust(as.matrix(distance_matrix), kmeans, method = "silhouette") +
  ggtitle("Silhouette Score for Optimal K")

# Gap statistic
#gap_stat <- clusGap(as.matrix(distance_matrix), FUN = kmeans, nstart = 25,
#                    K.max = 30, B = 50)
#fviz_gap_stat(gap_stat)

#~~~~~~~~~~~~~~~~~~~~~
# Source: https://github.com/LizaDarrous/PheWAS-cluster/blob/main/R/Clustering.R

# AIC 
# https://stackoverflow.com/questions/15839774/how-to-calculate-bic-for-k-means-clustering-in-r/25557162
kmeansAIC = function(fit){
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(data.frame(AIC = D + 2*m*k,
                    BIC = D + log(n)*m*k))
  
}

fit <- kmeans(x = distance_matrix, centers = 6)
kmeansAIC(fit)

cluster_list = list()
IC_df = data.frame(matrix(data=NA, nrow = 49, ncol=1))
colnames(IC_df) = c("AIC")
IC_df$nCluster = 2:50

BIC_df = data.frame(matrix(data=NA, nrow = 49, ncol=1))
colnames(BIC_df) = c("BIC")
BIC_df$nCluster = 2:50

for(i in 2:50){
  set.seed(240) # setting seed
  kmeans.re <- kmeans(distance_matrix, centers = i, nstart = 50, iter.max = 300)
  cluster_list[[length(cluster_list)+1]] = kmeans.re
  IC = kmeansAIC(kmeans.re)
  IC_df[(i-1),1] = IC$AIC
  BIC_df[(i-1),1] = IC$BIC
}


# AIC-cluster number identification for each observation
kmeans.minAIC = cluster_list[[which(IC_df$AIC==min(IC_df$AIC))]]
nClust.AIC = max(kmeans.minAIC$cluster)

# plot AIC
#pdf(file=paste0(res_dir,"/AIC-ncluster_",EXP_pheno,".pdf"), width=5, height = 7)
plot(IC_df$nCluster,IC_df$AIC, xlab="nCluster", ylab="AIC")
abline(v=nClust.AIC,col="red")
#dev.off()

# BIC-cluster number identification for each observation
kmeans.minBIC = cluster_list[[which(BIC_df$BIC==min(BIC_df$BIC))]]
nClust.BIC = max(kmeans.minBIC$cluster)

# plot AIC
#pdf(file=paste0(res_dir,"/AIC-ncluster_",EXP_pheno,".pdf"), width=5, height = 7)
plot(BIC_df$nCluster, BIC_df$AIC, xlab="nCluster", ylab="BIC")
abline(v=nClust.BIC,col="red")
#dev.off()

# assigning SNPs to clusters
AICclusters_rsid = list()
for(i in 1:nClust.AIC){
  AICclusters_rsid[[i]] = names(kmeans.minAIC$cluster)[which(kmeans.minAIC$cluster==i)]
}
print(paste0("Number of SNPs in each AIC grouped clusters: ", paste0(lengths(AICclusters_rsid), collapse = ", ")))
AICclusters_rsid_df = t(plyr::ldply(AICclusters_rsid, rbind))
#write.csv(AICclusters_rsid_df, paste0(res_dir,"/AICclusters_rsid_",EXP_pheno,".csv"), row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~
set.seed(123) # set a fixed seed for reproducibility, this will apply downwards

# Look at different numbers of clusters
k2 <- kmeans(as.matrix(distance_matrix), centers = 2, nstart = 25)
k3 <- kmeans(as.matrix(distance_matrix), centers = 3, nstart = 25)
k4 <- kmeans(as.matrix(distance_matrix), centers = 4, nstart = 25)
k5 <- kmeans(as.matrix(distance_matrix), centers = 5, nstart = 25)
k6 <- kmeans(as.matrix(distance_matrix), centers = 6, nstart = 25)
k7 <- kmeans(as.matrix(distance_matrix), centers = 7, nstart = 25)
k8 <- kmeans(as.matrix(distance_matrix), centers = 8, nstart = 25)
k9 <- kmeans(as.matrix(distance_matrix), centers = 9, nstart = 25)
k10 <- kmeans(as.matrix(distance_matrix), centers = 10, nstart = 25)

# plots to compare
p2 <- fviz_cluster(k2, geom = "point", data = as.matrix(distance_matrix)) + ggtitle("k = 2")
p3 <- fviz_cluster(k3, geom = "point", data = as.matrix(distance_matrix)) + ggtitle("k = 3")
p4 <- fviz_cluster(k4, geom = "point", data = as.matrix(distance_matrix)) + ggtitle("k = 4")
p5 <- fviz_cluster(k5, geom = "point", data = as.matrix(distance_matrix)) + ggtitle("k = 5")
p6 <- fviz_cluster(k6, geom = "point",  data = as.matrix(distance_matrix)) + ggtitle("k = 6")
p7 <- fviz_cluster(k7, geom = "point",  data = as.matrix(distance_matrix)) + ggtitle("k = 7")
p8 <- fviz_cluster(k8, geom = "point",  data = as.matrix(distance_matrix)) + ggtitle("k = 8")
p9 <- fviz_cluster(k9, geom = "point",  data = as.matrix(distance_matrix)) + ggtitle("k = 9")
p10 <- fviz_cluster(k10, geom = "point",  data = as.matrix(distance_matrix)) + ggtitle("k = 10")

library(gridExtra)
grid.arrange(p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 3)
#~~~~~~~~~~~~~~~~~~~~~

# Choose the optimal K (from the previous step, e.g., K = 2)
# for oily fish on TG, 8 has the most distinct clusters
optimal_k <- 2 # WHAT IS THIS NUMBER....

# Run K-means clustering
kmeans_result <- kmeans(as.matrix(distance_matrix), centers = optimal_k, nstart = 10)

# Add cluster assignments to a dataframe
cluster_df <- data.frame(Distance = as.vector(distance_matrix), Cluster = kmeans_result$cluster)

# Print first few rows
head(cluster_df)

table(kmeans_result$cluster)

# Create a dataframe with traits and their cluster assignments
clustered_traits <- data.frame(Trait = rownames(wide_matrix), Cluster = kmeans_result$cluster)

# View first few rows
head(clustered_traits)

# Count how many traits are in each cluster
table(clustered_traits$Cluster)

# Merge with descriptions
keyu <- unique(key)
clustered_traits$prefix <- sub("\\..*", "", clustered_traits$Trait)
dat <- merge(clustered_traits, keyu, by.x = "prefix", by.y = "V1", all.x = TRUE, all.y = FALSE)

# Save dat

#~~~~~~~~~~~~~~~~~~~~~~~~
# With dat, select one trait per cluster and begin formatting for MVMR
# We want to grab the trait with the highest heritability per cluster, so reload df 

# File with gcors from traits identified in PheWAS-Ttest method
gcor <- read.table("/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/processed_results/gcor_PheWASttest_results_all.txt")

# Header
gcor_header <- c("p1", "p2", "rg", "se", "z", "p", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se", "gcov_int", "gcov_int_se")

# Assign header to df
colnames(gcor) <- gcor_header

# Subset
gcor <- gcor %>% select(p2, h2_obs, h2_obs_se)

# Get trait ID
gcor$prefix <- sub("\\..*", "", gcor$p2)

# Get unique
gcor <- unique(gcor)

# Manually add heritability for 1448 trait (1st trait that gets left out in heritability table)
# See file: ldsc_results_starting_with_0.txt.log
# Heritability of phenotype 1
# ---------------------------
# Total Observed scale h2: 0.0047 (0.0015)
manual_row <- data.frame("p2" = "1448_4.ldsc.imputed_v3.both_sexes.tsv.bgz",
                         "h2_obs" = 0.0047,
                         "h2_obs_se"= 0.0015,
                         "prefix"= "1448_4")
gcor <- rbind(manual_row, gcor)

# Merge with dat
dat2 <- merge(dat, gcor, by = "prefix")

# Find the trait with the highest h2_obs within each Cluster
dat3 <- dat2 %>%
  group_by(Cluster) %>%
  slice_max(h2_obs, with_ties = FALSE)

print(dat3)

# Update p2 to be gwas
dat3$p2 <- gsub("ldsc", "gwas", dat3$p2)

# Save dat3
write.csv(dat3, "interim_data/ldsc/MVMR_NL_traits/oilyfish_GCST90239664_TG_Graham_TTest_MVMR_filenames.txt")
