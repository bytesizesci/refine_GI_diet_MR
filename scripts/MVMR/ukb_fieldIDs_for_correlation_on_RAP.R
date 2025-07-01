# Date: 06/16/25
# Goal: Get field IDs for table exporter on UKB Rap. For Supplemental Table of manuscript
#
# Author: KJS

# Notes: The pwas2 trait IDs include medical ids not in UKB's field ID pattern. We need to 
# remove these and format everything in pXXXXX_i0 format. The i0 is for instance 0, aka
# the baseline time point

# Get traits (file generated in PheWAS_Ttest_clumped_alcohol_GCST90319877_cir_Ghouse.R)
pwas2_trait_IDs <- read.csv("/pl/active/colelab/users/kjames/refinedMR/interim_data/PheWAS_Grep_SNPs/pwas2_trait_ids.csv")

