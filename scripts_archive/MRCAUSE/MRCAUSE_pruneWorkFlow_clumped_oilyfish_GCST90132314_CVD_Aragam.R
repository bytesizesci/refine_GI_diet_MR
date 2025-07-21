# Step 3 of MR-Cause: LD Pruning
# Parallelize 1 job per chromosome
# Sources: 
# LD Estimations from 1K Genomes CEU Population - https://zenodo.org/records/1464357#.W8a-fxROmV4
# MR Cause LD-pruning tutorial - https://jean997.github.io/cause/ldl_cad.html#LD_pruning_using_built_in_function

# Load required library
library(parallel)
library(cause)
library(dplyr)

# Set traits for naming files 
exposure <- "oilyfish"
outcome <- "GCST90132314_CVD_Aragam" 

# Set wd
wdir <- "/pl/active/colelab/users/kjames/enviroMR/"
setwd(wdir)

# Be in enviroMR 
dest <- "/pl/active/colelab/users/kjames/enviroMR/LDEstimated_1KG_CEU" # general
res <- paste0("/pl/active/colelab/users/kjames/enviroMR/LD_1KG_pruned/", exposure, "_", outcome) # trait specific
#inputdat <- "/pl/active/colelab/users/kjames/enviroMR/MR_CAUSE"

# Get ld_prune
source("scripts/ld_prune.R")

# Load X
X <- readRDS(paste0("interim_data/MRCAUSE/", exposure, "_", outcome, "_nonLDpruned.rds"))
variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))

# BUild a loop
chrN <- c(1:22)
for(chrN in 1:length(chrN)){
  # Download to Alpine
  #download.file(paste0("https://zenodo.org/record/1464357/files/chr",chrN,"_AF0.05_0.1.RDS?download=1"), destfile = paste0(dest, "/chr",chrN,"_AF0.05_0.1.RDS"))
  #download.file(paste0("https://zenodo.org/record/1464357/files/chr",chrN,"_AF0.05_snpdata.RDS?download=1"), destfile= paste0(dest, "/chr",chrN,"_AF0.05_snpdata.RDS"))
  
  # If already downloaded, comment out above and just load file
  # Load into R
  ld <- readRDS(paste0(dest,"/chr", chrN,"_AF0.05_0.1.RDS"))
  snp_info <- readRDS(paste0(dest,"/chr", chrN, "_AF0.05_snpdata.RDS"))
  
  # Prune with MR-cause function
  pruned <- ld_prune(variants = variants, 
                     ld = ld, 
                     total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))
  
  # How many independent snps?
  n_pruned_SNPs <- length(pruned)
  
  # Save
  write.csv(pruned, file = paste0(res, "/",exposure ,"_", outcome,"_nSNP",n_pruned_SNPs,"_chr", chrN, "_AF0.05.RDS"))
  
  # Optionally, you can delete the downloaded file if it's no longer needed
  #file.remove(file)
}

# Define the number of URLs to process
num_urls <- 22
#num_urls <- 2 #for getting the last 2 chr to run after having a memory problem

# Process URLs in parallel
cl <- makeCluster(detectCores())
clusterExport(cl, c("dest", "res", "chrN", "exposure", "outcome", "variants"))#deleted "inputdat" bc it was commented out above
clusterApply(cl, 1:num_urls, process_data)
stopCluster(cl)
