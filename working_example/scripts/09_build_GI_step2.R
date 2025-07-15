# Step 9: 09_build_GI_step2.R
# Goal: Make the GI for MVMR using the clumped files from Plink.

# Overall Strategy:
# - Gather all significant SNPs from the clumped files (snp_list)
# - Subset each full GWAS summary stat file by those SNPs
# - Combine by columns to make the GI

# Specifics:
# - Load GWAS sum-stat file, subset, save subset, unload GWAS sum-stat file
# - Repeat with next file. Note, some files, like outcome files, need special care.

# Load libraries
library(dplyr)
library(stringr)
library(purrr)

# ===================== USER CONFIGURATION =====================
# Set working direcotry
setwd("~/working_example")

# Set traits
EXP_pheno <- "alcohol"
OUT_pheno <- "ALT"
diet_trait <- paste0(EXP_pheno, "_", OUT_pheno)

# Define where the updated plink clumped files are
file_list <- paste0("~/working_example/results/plink_clumped/alcohol_ALT/GI_clumped/GI_flipped_clumped_results_chr", 1:22,".tsv.clumped")

# Set output file directory
output_dir <- paste0("~/working_example/results/GI/", diet_trait)
# Check if directory exists and make it if not
if (!dir.exists(output_dir)) {
  # If it doesn't exist, create it recursively
  dir.create(output_dir, recursive = TRUE)
  print(paste0("Nested directory '", output_dir, "' created successfully."))
} else {
  print(paste0("Nested directory '", output_dir, "' already exists."))
}

# Set diet file directory
data_dir <- paste0("~/working_example/data/")

# Define phewas traits
conf_traits <- c("30080_irnt",
                "30270_irnt",
                "23098_irnt",
                "20015_irnt",
                "23106_irnt")

# Define diet trait
diet_trait <- c("alcohol")

# Define diet path
diet_path <- paste0(data_dir, "/formatted_alcohol_clumped_full_sumstats.rds")
# ==============================================================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1: Load clumped data and save a SNP list
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load data
clumped_data <- lapply(file_list, function(f) read.table(f, header = TRUE)) %>% bind_rows()
length(unique(clumped_data$SNP)) # all unique
# This is what we will use for MVMR!

# Save
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE) #create folders if needed
write.table(clumped_data, 
            file = file.path(output_dir, paste0(output_dir,"/clumped_flipped_", diet_trait, "_MVMR_SNPs_GI.txt")),
            row.names = FALSE, quote = FALSE)

# Load if needed
#clumped_data <- read.table(paste0(output_dir,"/clumped_flipped_", diet_trait, "_MVMR_SNPs_GI.txt"), header = TRUE)

# Get a vector
snp_list <- clumped_data$SNP

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2: Subset the GWAS data to SNPs in snp_list
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define helper functions
# This might need to be adjusted based on the structure of your exposure data
subset_and_save_rds <- function(trait, path) {
  df <- readRDS(path)
  df_sub <- df[df$variant %in% snp_list, ]
  saveRDS(df_sub, file = file.path(output_dir, paste0("/df_", trait, "_subset.rds")))
  rm(df, df_sub); gc()
}

# This might need to be adjusted based on the structure of your exposure data
subset_and_save_txt <- function(trait) {
  infile <- file.path(paste0(data_dir ,"/GWAS_NL/",trait, ".gwas.imputed_v3.both_sexes.tsv"))
  df <- read.table(infile, header = TRUE)
  df_sub <- df[df$variant %in% snp_list, ]
  saveRDS(df_sub, file = file.path(output_dir, paste0("df_", trait, "_subset.rds")))
  rm(df, df_sub); gc()
}

# This might need to be adjusted based on the structure of your outcome data
subset_ALT_outcome <- function() {
  df <- readRDS(paste0(data_dir, "/formatted_GCST90013405_ALT_Pazoki.rds"))
  df$chr_pos <- str_replace(df$variant, "^([^:]+:[0-9]+):.*", "\\1")
  df_sub <- df[df$chr_pos %in% snp_list_noA0A1, ]
  saveRDS(df_sub, file = file.path(output_dir, "/df_ALT_subset.rds"))
  rm(df, df_sub); gc()
}

# Apply the functions to subset the full GWAS summary stats to the SNP list and 
# save the output. Note, this will take a while. 
walk2(conf_traits, subset_and_save_txt)
walk2(diet_trait, diet_paths, subset_and_save_rds)
subset_ALT_outcome()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 3: Combine subsets to create GI
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_traits <- c(conf_traits, diet_trait, "ALT")

load_and_standardize <- function(trait) {
  df <- readRDS(file.path(output_dir, paste0("df_", trait, "_subset.rds")))
  colnames(df) <- paste0(colnames(df), ".", trait) #add trait extension to column name
  df <- df %>%
    rename(variant = paste0("variant.", trait)) %>%
    mutate(chr_pos = str_replace(variant, "^([^:]+:[0-9]+):.*", "\\1"))
  return(df)
}

dfs <- lapply(all_traits, load_and_standardize)

mvmr_dat <- reduce(dfs, left_join, by = "chr_pos") %>%
  na.omit()

cat("Final dimensions:", dim(mvmr_dat), "\n")
