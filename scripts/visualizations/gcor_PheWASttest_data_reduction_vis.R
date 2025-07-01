# Date: 02/07/2025 (HBD Alaina)
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
count_1687 <- df %>% filter(p1=="1687.ldsc.imputed_v3.both_sexes.tsv.bgz") %>% tally()

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

# THANK GOD!

#######################################################
#######################################################
# Kept for historical sake 
# Helped me trouble shoot file sizes etc.
# Doesn't do anything now.
# Can delete this chunk in the future.

# Do you have all pairwise comparisons?
# Step 1: Get all unique items
#all_items <- unique(c(df$p1, df$p2))

# Step 2: Generate all possible pairwise comparisons (upper triangle only)
#all_pairs <- expand.grid(p1 = all_items, p2 = all_items) # this contains self-comparisons, remove next
#all_pairs$p1 <- as.character(all_pairs$p1)
#all_pairs$p2 <- as.character(all_pairs$p2)
#all_pairs_noself <- all_pairs %>% filter(p1 < p2) #7381

# Now, subset df
#filtered_df <- inner_join(df, all_pairs_noself, by = c("p1", "p2"))

# Count occurrences in both columns
#count_table <- df %>%
#  tidyr::pivot_longer(cols = c(p1, p2), values_to = "Trait") %>%  # Stack P1 and P2 into one column
#  count(Trait, name = "Count") %>%  # Count occurrences
#  arrange(desc(Count))  # Sort by frequency
# Well, they aren't the same... so clearly missing stuff

#count_table2 <- filtered_df %>%
#  tidyr::pivot_longer(cols = c(p1, p2), values_to = "Trait") %>%  # Stack P1 and P2 into one column
#  count(Trait, name = "Count") %>%  # Count occurrences
#  arrange(desc(Count))  # Sort by frequency
# Sigh, something is not right.

#Ensure p1 and p2 are characters so filer < will work below
#str(all_pairs$p1)
#all_pairs$p1 <- as.character(all_pairs$p1)
#all_pairs$p2 <- as.character(all_pairs$p2)
#all_pairs_top <- all_pairs %>% filter(p1 > p2)  # Keep only unique comparisons (no duplicates)

# Step 3: Find missing pairs
#missing_pairs <- anti_join(all_pairs_top, df, by = c("p1", "p2"))

# Print missing pairs
#print(missing_pairs)
#######################################################
#######################################################

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


# Read files in a loop
#for(i in seq_along(files)){
#  # Grep trait names from file name
#  exposure <- 
#  outcome <- 
#}
#############
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
df2 <- dfmu[dfmu$p1 %in% traitsu$V1,]
df <- df2

# Merge weights with df
df <- merge(df, weight_count, by.x = "p1id", by.y = "rmd_trait_id")

# Move descriptions to first two positions [order matters for clustering!]
df <- df %>% select(p1_description, p2_description, everything())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot to get the vibe
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Heat map
df_pfiltered <- df[df$p<0.05,]
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
  geom_edge_link(aes(width = trait_weight), alpha = 0.5, color = "gray50") +  # Edge thickness based on weight
  #geom_edge_link(aes(alpha = 0.5), color = "gray50") +  # Edges (relationships)
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






