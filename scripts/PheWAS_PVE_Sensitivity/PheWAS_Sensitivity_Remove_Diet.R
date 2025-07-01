# MR Analysis: PheWAS-PVE - Get diet IDs for sensitivity analysis 
# Author: KJS
# Date: 12/03/24 

# Goal: Get a list of dietary intake traits in the PheWAS to remove them 

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

#~~~~~~~~~~~~~~~~~~~~~
# 1. Load PheWAS data
#~~~~~~~~~~~~~~~~~~~~~
pwas <- read.csv("interim_data/PheWAS_Grep_SNPs/PheWAS_Oilyfish_WhiteXWhole_Alcochol_UpdatedProxies.csv")
pwas <- pwas %>% select(-c(X, V1))

# Make file name prettier. Get rid of extra junk.
pwas$trait <- sub("\\..*$", "", pwas$file_name)

# PWAS includes some duplicated traits that are in both raw and 
# inverse rank normalized transformed states. 
# In these cases, we want to keep the IRNT versions. 
# Remove the traits from pwas that end in _raw
pwas <- pwas[!grepl("_raw$", pwas$trait),]

#~~~~~~~~~~~~~~~~~~~~~
# 2. Filter PheWAS traits to those with an effective sample size of 50K
# For binary traits:
# Neff = (4 x Ncases x Ncontrols)/(Ncases + Ncontrols)
# For quantitative traits: 
# Neff = samplesize [... so this step is geared towards binary traits]
#~~~~~~~~~~~~~~~~~~~~~
md <- read.table(gzfile("../../../common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/metadata/phenotypes.both_sexes.v2.tsv.bgz"),header=TRUE, sep="\t", fill=TRUE, quote = "")
# quote="" important to get a full file

#~~~~~~~~~~~~~~~~~~~~~
# 2.1. Filter based on general sample size and ordinal variable sample size
#~~~~~~~~~~~~~~~~~~~~~
# from, https://github.com/LizaDarrous/PheWAS-cluster/blob/main/R/getSNPs_PheWAS.R
n_thresh <- 50000
md$n_eff = md$n_non_missing
ncont_ind = which(!is.na(md$n_cases)) #which rows have values in n_cases 
md$n_eff[ncont_ind] = (4*md$n_cases[ncont_ind]*md$n_controls[ncont_ind]) / (md$n_cases[ncont_ind]+md$n_controls[ncont_ind]) # (4*(control*case))/totN
small_n = which(md$n_eff < n_thresh)
md = md[-small_n,]

# Only keep traits in PheWAS that have sample size >50K
pwas2 <- pwas %>% filter((trait %in% md$phenotype))

#~~~~~~~~~~~~~~~~~~~~~
# 2.2 Distill to unique traits
#~~~~~~~~~~~~~~~~~~~~~
pwas_unique <- pwas2[!duplicated(pwas2$trait),]

#~~~~~~~~~~~~~~~~~~~~~
# 2.3. Merge pwas data with other meta data
#~~~~~~~~~~~~~~~~~~~~~
# Get the manifest key for PheWAS
key <- read.csv("/pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/metadata/UKBB GWAS Imputed v3 - File Manifest Release 20180731 - Description Lookup.csv",header = FALSE)
key <- key[!duplicated(key),]

# Merge to make a description column. Also, name the data frame traits to match LD's workflow.
trait_info <- merge(pwas_unique, key, by.x = "trait", by.y = "V1", all.x = TRUE, all.y = FALSE)

# new_name = old_name
trait_info = trait_info %>% dplyr::rename("description" = "V2")

# Only keep trait and description
trait_info = trait_info %>% select(trait, description)

#~~~~~~~~~~~~~~~~~~~~~
# 3. Save list, download, manually curate and upload to alpine
#~~~~~~~~~~~~~~~~~~~~~
write.csv(trait_info, "interim_data/PheWAS_PVE/PheWAS_traitinfo_manualcuration_sensitivity_rmdiet.csv")

