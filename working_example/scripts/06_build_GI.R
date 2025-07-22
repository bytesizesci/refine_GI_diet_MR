# Step 6: 06_build_GI.R
# Goal: Make the GI for MVMR using the clumped files from Plink.

# Load libraries
library(dplyr)

# ===================== USER CONFIGURATION =====================
# Set working direcotry
#setwd("/pl/active/colelab/users/kjames/refinedMR/")
setwd("~/working_example")

# Set traits
EXP_pheno <- "alcohol"
OUT_pheno <- "ALT"
diet_trait <- paste0(EXP_pheno, "_", OUT_pheno)

# Define PheWAS traits
# Identified in 04_snp_overlap.R
phewas_traits <- c("20015_irnt",
                   "23098_irnt",
                   "23106_irnt",
                   "30080_irnt",
                   "30270_irnt")

# Set clumped file directory
clumped_dir <- paste0("~/working_example/results/plink_clumped/",diet_trait)

# Set output file directory
output_dir <- paste0("~/working_example/results/GI_temp/", diet_trait) #this is created below

# Set diet file directory
data_dir <- paste0("~/working_example/data")

# ==============================================================

# Function to load clumped data for one trait
load_clumped_trait <- function(trait, diet_trait) {
  files <- paste0(clumped_dir, "/", trait,
                  "/clumped_results_", trait, "_chr", 1:22, ".tsv.clumped")
  dat <- lapply(files, function(f) read.table(f, header = TRUE)) %>%
    bind_rows()
  dat$pheno <- trait
  return(dat)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1: Load clumped data and combine
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load all candidate PheWAS traits (conf=confounder)
conf_snp <- do.call(rbind, lapply(phewas_traits, load_clumped_trait, diet_trait = diet_trait))
cat("Total SNPs loaded:", nrow(conf_snp), "\n")

# Save confounder SNPs
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE) #create folders if needed
write.table(conf_snp,
            file = file.path(output_dir, paste0(output_dir,"/confounder_clumped_snps_", diet_trait, "_MVMR_GI.txt")),
            row.names = FALSE, quote = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2: Load clumped diet data and combine with confounders
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load top hits from the diet exposure (e.g., alcohol)
diet_sig <- readRDS(data_dir, "/formatted_alcohol_clumped_noLD_nobiallelic.rds") %>%
  select(variant, pval.exposure)

# Rename columns
colnames(diet_sig) <- c("SNP", "P")

# Combine with confounder SNPs
conf_snp_sub <- conf_snp %>% select(SNP, P)
mvmr_snps <- bind_rows(conf_snp_sub, diet_sig)

# Save
write.table(mvmr_snps,
            file = file.path(output_dir, paste0(output_dir, "/all_clumped_snps_", diet_trait, "_MVMR_GI.txt")),
            row.names = FALSE, quote = FALSE)


