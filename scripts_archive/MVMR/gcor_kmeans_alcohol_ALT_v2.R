# Date: 03/28/2025 
# Updated: 04/21/2025
# Goal: Process gcor results to identify traits for MVMR
# Using Gcor results of PheWAS-Ttest snp-trait outliers. We tried using k-means
# but now we are going to lean into hierarchical clustering. 
#
# Author: KJS

#Set working direcotry
setwd("/pl/active/colelab/users/kjames/refinedMR/")

# Set plot directory for saving
res_dir <- "/pl/active/colelab/users/kjames/refinedMR/interim_data/MVMR/oilyfish_tg/plots"

#install.packages("ggraph")
#install.packages("tidygraph")  # Needed for network visualization
#install.packages("igraph") # for louvain cluster
library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(ggplot2)
library(cluster)   # For silhouette score
library(factoextra) # For visualization
library(tidyr)
library(RColorBrewer)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# File with gcors from traits identified in PheWAS-Ttest method
df <- read.table("/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/processed_results/gcor_PheWASttest_results_all.txt")

# Header
gcor_header <- c("p1", "p2", "rg", "se", "z", "p", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se", "gcov_int", "gcov_int_se")

# Assign header to df
colnames(df) <- gcor_header

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

# Ensure all pairs are in the upper triangle (P1 < P2) #4/21/25-think this causes errors between p1, p1desription, and p1id; try a fix below
# Have to redo with the description columns
#dfmu2 <- dfmu %>%
#  mutate(P1_desc_new = pmin(p1_description, p2_description),  # Assign the smaller value to P1
#         P2_desc_new = pmax(p1_description, p2_description)) %>%  # Assign the larger value to P2
#  select(P1_desc_new, P2_desc_new, p2id, p1id, p1, p2, rg) %>%
#  rename(p1_description = P1_desc_new, p2_description = P2_desc_new)

# Fix that maintains the relationship of all 3 types of id columns
dfmu2 <- dfmu %>%
  mutate(swap = p1_description > p2_description) %>%
  mutate(
    p1_final = if_else(swap, p2, p1),
    p2_final = if_else(swap, p1, p2),
    p1id_final = if_else(swap, p2id, p1id),
    p2id_final = if_else(swap, p1id, p2id),
    p1_description_final = if_else(swap, p2_description, p1_description),
    p2_description_final = if_else(swap, p1_description, p2_description)
  ) %>%
  select(p1_final, p2_final, p1id_final, p2id_final, p1_description_final, p2_description_final, rg) %>%
  rename(
    p1 = p1_final, p2 = p2_final,
    p1id = p1id_final, p2id = p2id_final,
    p1_description = p1_description_final,
    p2_description = p2_description_final
  )

ggplot(dfmu2, aes(x = p1_description, y = p2_description, fill = rg)) +
  geom_tile() +  # Heatmap tiles
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Color scale
  theme_minimal() +
  labs(title = "Genetic Correlations",
       x = "P1",
       y = "P2",
       fill = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define exposure and outcome of analysis
EXP_pheno <- "alcohol"
OUT_pheno <- "ALT"

# File with traits specific to an MR analysis
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/oilyfish_GCST90239664_TG_Graham_descending_TTest_filenames.txt", header = FALSE)
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/oilyfish_GCST90132314_CVD_Aragam_descending_TTest_filenames.txt", header = FALSE)
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/oilyfish_GCST006901_Height_Yengo_2018_descending_TTest_filenames.txt", header = FALSE)
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/whiteXwholegrain_GCST90239658_LDL_Graham_descending_TTest_filenames.txt", header = FALSE)
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/whiteXwholegrain_GCST90132314_CVD_Aragam_descending_TTest_filenames.txt", header = FALSE)
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/whiteXwholegrain_GCST006901_Height_Yengo_2018_descending_TTest_filenames.txt", header = FALSE)
traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/alcohol_GCST90013405_ALT_Pazoki_descending_TTest_filenames.txt", header = FALSE)
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/alcohol_GCST90319877_cir_Ghouse_descending_TTest_filenames.txt", header = FALSE)
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/alcohol_GCST006901_Height_Yengo_2018_descending_TTest_filenames.txt", header = FALSE)

# Get unique traits
traitsu <- unique(traits)
traitsu$V1 <- gsub("/pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/LDSC/", "", traitsu$V1)

# Filter df to traitsu
df2 <- dfmu2[dfmu2$p1 %in% traitsu$V1 & dfmu2$p2 %in% traitsu$V1,] #use AND not OR
df <- df2


# Fix that maintains the relationship of all 3 types of id columns
df <- df %>%
  mutate(swap = p1_description > p2_description) %>%
  mutate(
    p1_final = if_else(swap, p2, p1),
    p2_final = if_else(swap, p1, p2),
    p1id_final = if_else(swap, p2id, p1id),
    p2id_final = if_else(swap, p1id, p2id),
    p1_description_final = if_else(swap, p2_description, p1_description),
    p2_description_final = if_else(swap, p1_description, p2_description)
  ) %>%
  select(p1_final, p2_final, p1id_final, p2id_final, p1_description_final, p2_description_final, rg) %>%
  rename(
    p1 = p1_final, p2 = p2_final,
    p1id = p1id_final, p2id = p2id_final,
    p1_description = p1_description_final,
    p2_description = p2_description_final
  )

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

# We can do better than this graph - we want the tiles to be organized by 
# cluster and not by alphabetical order

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make wide data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4/21/25 - Having trouble with the chunk below making a complete matrix. 
# Try with new code
# Make sure p1-p2 and p2-p1 exists (dcast has issues if not)
df_full <- df %>%
  mutate(p1_alt = p2, p2_alt = p1) %>%
  rename(rg_alt = rg) %>%
  select(p1 = p1_alt, p2 = p2_alt, rg = rg_alt) %>%
  bind_rows(df)  # combine original and flipped
# Now cast
wide_df_test <- reshape2::dcast(df_full, p1 ~ p2, value.var = "rg")

# Set rownames
rownames(wide_df_test) <- wide_df_test[,"p1"]
wide_df_test <- wide_df_test %>% select(-p1)

# Reorder rows and columns
#rownames_matched <- wide_df_test$p1
#wide_mat <- wide_df_test %>% column_to_rownames("p1") %>% as.matrix()
#wide_mat <- wide_mat[rownames_matched, rownames_matched]  # enforce symmetry

wide_matrix <- wide_df_test
str(wide_matrix)

# Add diagonal
diag(wide_matrix) <- 1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Hierarchical clustering plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(789678)

# Side adventure 
# Plot based on rg value instead of alphabetical order
hc = hclust(as.dist(1 - wide_matrix), method = "complete") # hclust(dist(wide_matrix)) - dist() is for calculating distances but I just need to convert my correlations (similarity) into distances
ordered_traits = hc$labels[hc$order]
# Add clustering-based factor levels
dfplot <- df %>% mutate(p1 = factor(p1, levels = ordered_traits),
                        p2 = factor(p2, levels = ordered_traits))

# Get unique mapping of trait -> label
#label_map <- df %>% select(p1, p1_description) %>% distinct() #OLD - fixed below 
#label_map <- setNames(label_map$p1_description, label_map$p1) #OLD
#dim(label_map) # this doesn't match wide_matrix. why? #OLD
#setdiff(rownames(wide_df_test), names(label_map)) #OLD
# This is missing... 23102_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz - whole body water mass

# To fix, build label map from p1 and p2! (this one trait is only represented in p2)
label_map <- df %>%
  select(p1, p1_description) %>%
  rename(trait = p1, label = p1_description) %>%
  bind_rows(
    df %>% select(p2, p2_description) %>%
      rename(trait = p2, label = p2_description)
  ) %>%
  distinct(trait, label) %>%
  tibble::deframe()  # Makes it a named vector


# This is ordered by clustering but has results in both triangles
ggplot(dfplot, aes(x = p1, y = p2, fill = rg)) +
  geom_tile() +
  scale_x_discrete(labels = label_map) +
  scale_y_discrete(labels = label_map) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Prep to get all values in one triangle
# Reorder matrix based on clustering
wide_ordered <- wide_matrix[ordered_traits, ordered_traits]

# Keep upper triangle and convert to long format
wide_ordered[lower.tri(wide_ordered)] <- NA

# For plotting, turn back into long data
# Turn into matrix
str(wide_ordered)
wide_ordered <- wide_ordered %>% as.matrix()
# Turn long
df_upper <- as.data.frame(as.table(wide_ordered)) %>%
  filter(!is.na(Freq)) %>%
  rename(p1 = Var1, p2 = Var2, rg = Freq)

# Plot
ggplot(df_upper, aes(x = p1, y = p2, fill = rg)) +
  geom_tile() +
  scale_x_discrete(labels = label_map) +
  scale_y_discrete(labels = label_map) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/pl/active/colelab/users/kjames/refinedMR/interim_data/MVMR/alcohol_ALT/plots/alcohol_alt_hc_complete_heatmap.png")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Hierarchical clustering
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert correlation matrix to distance metric (1 - correlation)
# euclidean by default
distance_matrix <- as.dist(1 - wide_matrix)
hc <- hclust(distance_matrix, method = "complete") # ward.D2, complete, etc.

# Preserve a copy so that you can keep the file names, which are needed for the cluster step
hc_fn <- hc
clusters_filename <- cutree(hc_fn, k=5)
cluster_df_filename <- data.frame(
  trait_id = names(clusters_filename),
  cluster = clusters_filename
)
write.csv(cluster_df_filename, paste0("interim_data/MVMR/",EXP_pheno,"_",OUT_pheno,"/hclust_complete_cut5.csv"), row.names = FALSE)

# This for readability
# Add labels - helps for plotting
hc$labels <- label_map[hc$labels] #label_map built above

# Plot tree
plot(hc, main = "Trait Clustering Based on Genetic Correlation")
abline(h=1.15, col="red", lty=2)

# Cut tree
clusters <- cutree(hc, k = 5)  # Choose number of groups

# Combine with labels
cluster_df <- data.frame(
  trait_id = names(clusters),
  cluster = clusters
)

# Plot with colors
dend <- as.dendrogram(hc)
dend <- dendextend::color_branches(dend, k = 5)
plot(dend, cex = 0.5, main = "")

# Plot with pheatmap
png("/pl/active/colelab/users/kjames/refinedMR/interim_data/MVMR/alcohol_ALT/plots/alcohol_alt_hc_complete_pheatmap.png", width = 3000, height = 2000, res = 200)
pheatmap::pheatmap(wide_matrix,
         labels_row = label_map[rownames(wide_matrix)],
         labels_col = label_map[colnames(wide_matrix)],
         clustering_distance_rows = as.dist(1 - abs(wide_matrix)),
         clustering_distance_cols = as.dist(1 - abs(wide_matrix)),
         clustering_method = "complete",
         #col = brewer.pal(10,"Spectral"),
         main = "")
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract traits for MVMR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# One idea for deciding what trait to extract is... 
# Look at the number of significant SNPs that overlap with the GI for diet and select the trait
# with the most













#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# k-means
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert correlation matrix to distance metric (1 - correlation)
# euclidean by default
distance_matrix <- as.dist(1 - wide_matrix)

# plot directory for saving
res_dir <- "/pl/active/colelab/users/kjames/refinedMR/interim_data/MVMR/oilyfish_tg/plots"

# How many clusters?
# wss=within cluster sums of squares
# plot elbow method
pdf(file=paste0(res_dir,"/fvis_wss_elbow_ncluster_",EXP_pheno,"_",OUT_pheno,".pdf"), width=7, height = 5)
fviz_nbclust(as.matrix(distance_matrix), kmeans, method = "wss") +
  ggtitle("Elbow Method for Optimal K")
dev.off()

# silhouette=average sihouette width; high is good 
pdf(file=paste0(res_dir,"/fvis_silhouette_ncluster_",EXP_pheno,"_",OUT_pheno,".pdf"), width=7, height = 5)
fviz_nbclust(as.matrix(distance_matrix), kmeans, method = "silhouette") +
  ggtitle("Silhouette Score for Optimal K")
dev.off()

# AIC
# Source: https://github.com/LizaDarrous/PheWAS-cluster/blob/main/R/Clustering.R
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
nClust.AIC

# plot directory for saving
res_dir <- "/pl/active/colelab/users/kjames/refinedMR/interim_data/MVMR/oilyfish_tg/plots"

# plot AIC
pdf(file=paste0(res_dir,"/AIC-ncluster_",EXP_pheno,"_",OUT_pheno,".pdf"), width=5, height = 7)
plot(IC_df$nCluster,IC_df$AIC, xlab="nCluster", ylab="AIC")
abline(v=nClust.AIC,col="red")
dev.off()

# BIC-cluster number identification for each observation
kmeans.minBIC = cluster_list[[which(BIC_df$BIC==min(BIC_df$BIC))]]
nClust.BIC = max(kmeans.minBIC$cluster)

# plot AIC
pdf(file=paste0(res_dir,"/BIC-ncluster_",EXP_pheno,"_",OUT_pheno,".pdf"), width=5, height = 7)
plot(BIC_df$nCluster, BIC_df$AIC, xlab="nCluster", ylab="BIC")
abline(v=nClust.BIC,col="red")
dev.off()

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
p2 <- fviz_cluster(k2, geom = "point", data = as.matrix(distance_matrix), label = rownames(wide_matrix)) + ggtitle("k = 2")
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

# Save
png("/pl/active/colelab/users/kjames/refinedMR/interim_data/MVMR/oilyfish_tg/plots/oilyfish_tg_kmeans_kclusters.png", width = 3000, height = 2000, res = 200)
grid.arrange(p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 3)
dev.off()
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
write.csv(dat3, "interim_data/ldsc/MVMR_NL_traits/alcohol_ALT_TTest_kmeans2_MVMR_filenames.txt")



