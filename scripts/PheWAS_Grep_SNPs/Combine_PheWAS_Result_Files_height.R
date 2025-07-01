# Goal: Combine PheWAS grepped files into one master file for height
# Author: KJS
# Date: November 20, 2024

# Run locally
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

#Merge the files together
library(data.table)

#Read files with full path
#files <- list.files("/pl/active/colelab/users/kjames/enviroMR/results_array", full.names = TRUE)
#files <- list.files("/pl/active/colelab/users/kjames/enviroMR/results_PheWAS/freshfruit_clumped_GCST90310295_DBP", full.names = TRUE)
#files <- list.files("/pl/active/colelab/users/kjames/enviroMR/results_PheWAS/RefiningGI_DietMR_Alcochol_WhiteXWhole_OilyFish", full.names = TRUE)
#files <- list.files("/pl/active/colelab/users/kjames/refinedMR/interim_data/PheWAS_Grep_SNPs/PheWAS_SNPs_Oilyfish_WhiteXWhole_Alcochol_UpdatedProxies/", full.names = TRUE)
files <- list.files("/pl/active/colelab/users/kjames/refinedMR/interim_data/PheWAS_Grep_SNPs/PheWAS_SNPs_Height_Proxies/", full.names = TRUE)

#Combine all data
data <- rbindlist(lapply(files, fread), fill = TRUE)

#Save
#write.csv(data, "/pl/active/colelab/users/kjames/enviroMR/results/SNP_PVE_freshfruit_clumped_GCST90310295_DBP_CombinedArrays_240821.csv")
#write.csv(data, "/pl/active/colelab/users/kjames/refinedMR/interim_data/PheWAS_Grep_SNPs/PheWAS_Oilyfish_WhiteXWhole_Alcochol_UpdatedProxies.csv")
write.csv(data, "/pl/active/colelab/users/kjames/refinedMR/interim_data/PheWAS_Grep_SNPs/PheWAS_Height_Outcome.csv")
