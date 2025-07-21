# Updated: 11/12/24 - with new proxy method, need to regrep SNPs; switches wd to refinedMR
# Updated: 9/17/24
# Update: Generalizing workflow to run for MR Diet Refinement paper with Julie
# Added code that specifies exposure
# Author: KJS

# R Script that loads the .bgz Neale results
# Finds desired SNPs, calculates PVE, and saves the output.
# Output is processed in a subsequent workflow.
# When we have a final file, we can sort by SNP, see which pheno is the 
# most explanatory and decide whether to keep it or not

# Used in an array. Runs 1 job.

#Check that the packages have been installed 
#In the terminal
#acompile --> module load R/4.2.2 --> R --> check/install
#need dplyr and R.utils

#library(R.utils) #specify instead of load globally
library(dplyr)

#Set working direcotry
setwd("/pl/active/colelab/users/kjames/refinedMR/")

# Set traits for naming files 
exposure <- "PheWAS_SNPs_Oilyfish_WhiteXWhole_Alcochol_UpdatedProxies"

#Load SNPs that you want to query (gv=genetic variants)
gv <- read.csv("/pl/active/colelab/users/kjames/refinedMR/interim_data/merged_GI/PheWAS_variants_oilyfish_whiteXwholegrain_alcohol_updated_exposure_proxies.csv", row.names = 1)
head(gv)
names(gv) <- "variant.exposure"

#Get TASK_ID - refers to which file to take
taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
str(taskid)
taskid_num <- as.numeric(taskid)
print(paste0("the number processed was... ", taskid))

#Get list of files in GWAS both_sexes directory
#4586 files... many are duplicates (i.e., run raw and inverse rank normalized)

# 11/13/24 - The files are already unzipped and stored in the following path in gwas_files2
#gwas_files <- list.files("/pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/GWAS2", pattern = "*both_sexes*", full.names = TRUE)
#head(gwas_files)
#str(gwas_files)
gwas_files <- list.files("/pl/active/colelab/users/kjames/enviroMR/temp", full.names = TRUE)
head(gwas_files)
str(gwas_files)

#Specify file
gwas_file <- gwas_files[taskid_num]
str(gwas_file)
print(paste0("the file processed was... ", gwas_file))

# Dev check in .out file
print("Heading to the calculations...")

#Rename
file_i <- gwas_file

#Make a name for the file destination
file_dest <- gsub(".bgz", "", file_i)
file_dest <- gsub("/pl/active/colelab/users/kjames/enviroMR/temp/", "", file_dest) # Update - adjusted here to show new path, not common/published_gwas path
file_dest #should just be the filename.tsv now

#Step 1: Unzip the file
#Unzip file using gunzip
#df_unzip <- R.utils::gunzip(file_i, destname = paste0("temp/",file_dest), skip = TRUE) #produces a warning but its okay

# 11/13/24 - files are already unzipped in enviroMR location; use them for this rerun rather than restoring/taking more storage
df_unzip_enviroMR <- paste0("/pl/active/colelab/users/kjames/enviroMR/temp/", file_dest)

#Step 2: Load the tsv using read.table
#df_unzipped <- read.table(df_unzip, sep = "\t")
df_unzipped <- read.table(df_unzip_enviroMR, sep = "\t") #use this if you set the wd to something other than enviroMR; also see line 93-94

#Set colnames
colnames(df_unzipped) <- df_unzipped[1,]
df_unzipped <- df_unzipped[-1,]

# Added this variant focused approach because the chr:pos:AX:AX notation was causing errors. 
# Will take caution when performing the filtering/MR but can proceed with grepping SNPs 
#df_unzipped <- df_unzipped %>% rowwise() %>% mutate(chr_pos = paste0(str_split(variant, '\\:')[[1]][1], ':', str_split(variant, '\\:')[[1]][2])) %>% ungroup() #this solution is brutally slow
#gv <- gv %>% rowwise() %>% mutate(chr_pos = paste0(str_split(variant.exposure, '\\:')[[1]][1], ':', str_split(variant.exposure, '\\:')[[1]][2])) %>% ungroup()
df_unzipped$chr_pos <- regmatches(df_unzipped$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", df_unzipped$variant, ignore.case = TRUE)) 
gv$chr_pos <- regmatches(gv$variant.exposure, regexpr("[a-z0-9]+\\:[a-z0-9]+", gv$variant.exposure, ignore.case = TRUE)) 

# From a local test wher I loaded in a .bgz file and looked at the expected outcome from the grep
#df_sub <- df_unzipped %>% filter(chr_pos %in% gv$chr_pos) #205 variants
#missing <- gv %>% anti_join(df_sub, by="chr_pos") # 8 missing
#length(unique(gv$chr_pos)) #235

#Subset to desired SNPs
df_sub <- df_unzipped %>% filter(chr_pos %in% gv$chr_pos) #228 variants; missing 7
#df_sub <- df_unzipped %>% filter(variant %in% gv$variant.exposure) #gv.V2 if using older read.table inputs

#Add col for filename
#This will be used to figure out which pheno the row corresponds to
#file_name <- gsub("/pl/active/colelab/users/kjames/enviroMR/temp/", "", df_unzip_enviroMR)
#file_name <- gsub("temp/", "", df_unzip)
df_sub$file_name <- file_dest

#Make data numeric (not all cols, just the ones supposed to be numeric)
df_sub$minor_AF <- as.numeric(df_sub$minor_AF)
df_sub$n_complete_samples <- as.numeric(df_sub$n_complete_samples)
df_sub$AC <- as.numeric(df_sub$AC)
df_sub$ytx <- as.numeric(df_sub$ytx)
df_sub$beta <- as.numeric(df_sub$beta)
df_sub$se <- as.numeric(df_sub$se)
df_sub$tstat <- as.numeric(df_sub$tstat)
df_sub$pval <- as.numeric(df_sub$pval)

#Calculate percent variance explained
#Proportion of variance in phenotype explained by a given 
#SNP (PVE) = [2*(beta^2)*MAF*(1-MAF)]/[2*(beta^2)*MAF(1-MAF)+((se(beta))^2)*2*N*MAF*(1-MAF)]
#
#Source: Shim, H., Chasman, D.I., Smith, J.D., Mora, S., Ridker, P.M., Nickerson, D.A., Krauss, R.M.,
#and Stephens, M. (2015). A multivariate genome-wide association analysis of 10 LDL subfractions,
#and their response to statin treatment, in 1868 Caucasians. PLoS One 10, e0120758.
#
#Define variables to make equation simpler to type
MAF <- df_sub$minor_AF
N  <- df_sub$n_complete_samples
se <- df_sub$se
beta <- df_sub$beta
#Calculate PVE variable
numerator <- 2*(beta^2)*MAF*(1-MAF)
denominator <- (2*(beta^2)*MAF*(1-MAF)+((se^2)*2*N*MAF*(1-MAF))) #fixed since _v2
df_sub$PVE <- numerator/denominator

#Remove giant unzipped file and start the loop again
#rm(df_unzip)
#rm(df_unzipped)
#Garbage collection
#https://stackoverflow.com/questions/1467201/forcing-garbage-collection-to-run-in-r-with-the-gc-command
#gc() 

#Save
write.csv(df_sub, paste0("/pl/active/colelab/users/kjames/refinedMR/interim_data/PheWAS_Grep_SNPs/", exposure, "/", exposure, "_", taskid,".csv"))