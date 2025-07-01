# Manhattan plot


# Libraries
library(ggplot2)
library(tidyr)


#file_list <- paste0("interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/4100_irnt/clumped_results_4100_irnt_chr", c(1:20, 22),".tsv.clumped") # selected in kmeans 2, ankle spacing width
#dfsig <- lapply(file_list, function(f) read.table(f, header = TRUE)) %>% bind_rows()

# Significant SNPs
file_list <- paste0("interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/50_irnt/clumped_results_50_irnt_chr", 1:22,".tsv.clumped") # selected in kmeans 2, standing height
dfsig <- lapply(file_list, function(f) read.table(f, header = TRUE)) %>% bind_rows()
head(dfsig)
dim(dfsig)
# Get chr and bp columns
dfsig <- separate_wider_delim(dfsig, cols = SNP, delim = ":", names = c("chr", "bp", "ref", "alt"))
# Convert 'bp' to numeric for proper data handling
dfsig$bp <- as.numeric(dfsig$bp)

# 50_irnt
# Full summary stats
df_50_full <- read.table("/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/50_irnt.gwas.imputed_v3.both_sexes.tsv", header=TRUE)
# Extract SNPs from GWAS
df_50_subset <- df_50_full[df_50_full$variant %in% snp_list,]
# Save subsetted SNPs in case R runs out of memory and closes
saveRDS(df_50_subset, "interim_data/MVMR/oilyfish_tg/df_50_subset.rds")
# Remove full summary stats to clear space
rm(df_50_full)

df <- df_50_full

# Get chr and bp columns
df <- separate_wider_delim(df, cols = variant, delim = ":", names = c("chr", "bp", "ref", "alt"))

# Convert 'bp' to numeric for proper data handling
df$bp <- as.numeric(df$bp)

# Create a cumulative position for each SNP
df <- df %>%
  group_by(chr) %>%
  mutate(chr_len = max(bp)) %>%
  ungroup() %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  mutate(BPcum = bp + tot)

# Prepare X-axis labels
axisdf <- df %>%
  group_by(chr) %>%
  summarize(center = mean(BPcum))

# log10 pvalue
df$log10p <- -log10(df$pval)

# Subset - makes graphing manageable
dfs <- df[df$pval < 0.00005,]

# Plot
ggplot(dfs, aes(x = BPcum, y = log10p, color = as.factor(chr))) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = rep(c("gray10", "gray60"), 12)) +
  scale_x_continuous(label = axisdf$chr, breaks = axisdf$center) +
  labs(x = "Chromosome", y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

# Check that the spacing between the significant hits is indeed at least 50K
# Assuming your dataframe is named df with columns 'chr' and 'bp'
filtered_df <- dfsig %>%
  group_by(chr) %>%                      # Group by chromosome
  arrange(bp, .by_group = TRUE) %>%      # Arrange positions within each chromosome
  filter(row_number() == 1 | bp - lag(bp) >= 50000) %>%
  ungroup()


