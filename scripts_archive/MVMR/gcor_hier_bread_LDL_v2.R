# Date: 03/28/2025 
# Updated: 04/21/2025
# Updated: 05/10/2025
# Goal: Process gcor results to identify traits for MVMR
# Using Gcor results of PheWAS-Ttest snp-trait outliers. 
# Approach is based on hierarchical clustering. 
#
# Author: KJS
#
# Traits: bread type on LDL-C

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET UP

#Set working directory
setwd("/pl/active/colelab/users/kjames/refinedMR/")

# Set plot directory for saving
res_dir <- "/pl/active/colelab/users/kjames/refinedMR/interim_data/MVMR/bread_LDLc/plots"

# Load libraries
library(dplyr)
library(ggplot2)
library(factoextra) # For visualization
library(tidyr)
library(RColorBrewer)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GENETIC CORRELATION PROCESSING

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
# Yes! Complete upper triangle.

# Replace file names with trait descriptions
key <- read.csv("/pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/metadata/UKBB GWAS Imputed v3 - File Manifest Release 20180731 - Description Lookup.csv",header = FALSE)

# Rename dfupper df to keep code simpler
df <- df_upper

# Isolate the ID portion from file name
df$p1id <- sub("\\..*","", df$p1)
df$p2id <- sub("\\..*","", df$p2)

# Merge to get human readable names
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
# Solution that maintains the relationship of all 3 types of id columns
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

# Another visual check
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
# DIET-HEALTH OUTCOME WRANGLING AND CLUSTERING

# Define exposure and outcome of analysis
EXP_pheno <- "whiteXwholegrain"
OUT_pheno <- "LDL"

# File with traits specific to an MR analysis
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/oilyfish_GCST90239664_TG_Graham_descending_TTest_filenames.txt", header = FALSE)
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/oilyfish_GCST90132314_CVD_Aragam_descending_TTest_filenames.txt", header = FALSE)
#traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/oilyfish_GCST006901_Height_Yengo_2018_descending_TTest_filenames.txt", header = FALSE)
traits <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/whiteXwholegrain_GCST90239658_LDL_Graham_descending_TTest_filenames.txt", header = FALSE)
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
df <- df2 # simplify name

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
# Make wide data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make a complete matrix
# Make sure p1-p2 and p2-p1 exists (dcast has issues if not)
df_full <- df %>%
  mutate(p1_alt = p2, p2_alt = p1) %>%
  rename(rg_alt = rg) %>%
  select(p1 = p1_alt, p2 = p2_alt, rg = rg_alt) %>%
  bind_rows(df)  # combine original and flipped

# Now cast
wide_matrix <- reshape2::dcast(df_full, p1 ~ p2, value.var = "rg")

# Set rownames
rownames(wide_matrix) <- wide_matrix[,"p1"]
wide_matrix <- wide_matrix %>% select(-p1)

# Check structure
str(wide_matrix)

# Add diagonal
diag(wide_matrix) <- 1

# Save
write.csv(wide_matrix, paste0(res_dir,"/",EXP_pheno,"_",OUT_pheno,"_gcor_wide_matrix.csv"), row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HIERARCHICAL CLUSTERING VISUALS
# Consider how many clusters 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Good habit
set.seed(789678)

# Plot based on rg value instead of alphabetical order
hc = hclust(as.dist(1 - wide_matrix), method = "complete") # note: hclust(dist(wide_matrix)) - dist() is for calculating distances but I just need to convert my correlations (similarity) into distances
ordered_traits = hc$labels[hc$order]

# How many clusters?
plot(rev(hc$height), type = "b") #reverse of height
# Silhouette method
factoextra::fviz_nbclust(matrix(as.dist(1 - wide_matrix)), FUN = hcut, method = "silhouette")

# Add clustering-based factor levels
dfplot <- df %>% mutate(p1 = factor(p1, levels = ordered_traits),
                        p2 = factor(p2, levels = ordered_traits))

# Get unique mapping of trait -> label
# Build label map from p1 and p2! (this one trait is only represented in p2)
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

# Save
ggsave(paste0(res_dir,"/",EXP_pheno,"_", OUT_pheno, "_alt_hc_complete_heatmap.png"),
       width = 20,
       height = 16,
       units = "in",
       dpi = 300,
       limitsize = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HIERARCHICAL CLUSTERING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert correlation matrix to distance metric (1 - correlation)
# euclidean by default
distance_matrix <- as.dist(1 - wide_matrix)
hc <- hclust(distance_matrix, method = "complete") 

# Preserve a copy so that you can keep the file names, which are needed for the cluster step
hc_fn <- hc

# From plots above, define how many clusters
k=3 # define how many clusters

# Cut tree
clusters_filename <- cutree(hc_fn, k=k)

# Turn into data frame for results processing
cluster_df_filename <- data.frame(
  trait_id = names(clusters_filename),
  cluster = clusters_filename
)

# Save
write.csv(cluster_df_filename, paste0(res_dir,"/",EXP_pheno,"_",OUT_pheno,"_hclust_complete_cut_",k,".csv"), row.names = FALSE)

# This for readability
# Add labels - helps for plotting
hc$labels <- label_map[hc$labels] #label_map built above

# Plot tree
plot(hc, main = "Trait Clustering Based on Genetic Correlation")
abline(h=1.15, col="red", lty=2)

# Cut tree
clusters <- cutree(hc, k = k)  # Choose number of groups

# Combine with labels
cluster_df <- data.frame(
  trait_id = names(clusters),
  cluster = clusters
)

# Plot with colors
dend <- as.dendrogram(hc)
dend <- dendextend::color_branches(dend, k = k)
plot(dend, cex = 0.5, main = "")

# Plot with pheatmap
png(paste0(res_dir,"/",EXP_pheno,"_",OUT_pheno,"_hclust_complete_cut_",k,"_pheatmap.png"),
    width = 3000,
    height = 2000,
    res = 200)
pheatmap::pheatmap(wide_matrix,
         labels_row = label_map[rownames(wide_matrix)],
         labels_col = label_map[colnames(wide_matrix)],
         clustering_distance_rows = as.dist(1 - abs(wide_matrix)),
         clustering_distance_cols = as.dist(1 - abs(wide_matrix)),
         clustering_method = "complete",
         #col = brewer.pal(10,"Spectral"),
         main = "")
dev.off()



