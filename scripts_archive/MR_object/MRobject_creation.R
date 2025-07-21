# Goal: Create an MRobject for the refining GI of MR project - updated with new proxy exposure SNPs and neg controls 
# Author: KJS
# Date: November 13, 2024

# Create list object for input data
df1 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_oilyfish_GCST90239664_TG_Graham_GRCh37.rds")
df2 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_oilyfish_GCST90132314_CVD_Aragam_GRCh37.rds")
df3 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_oilyfish_GCST006901_Height_Yengo_2018.rds")

df4 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_whiteXwholegrain_GCST90239658_LDL_Graham.rds")
df5 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_whiteXwholegrain_GCST90132314_CVD_Aragam.rds")
df6 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_whiteXwholegrain_GCST006901_Height_Yengo_2018.rds")

df7 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_alcohol_GCST90013405_ALT_Pazoki.rds")
df8 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_alcohol_GCST90319877_cir_Ghouse.rds")
df9 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_alcohol_GCST006901_Height_Yengo_2018.rds")

# MRobj
mrobj <- list()
mrobj[[1]] <- df1
mrobj[[2]] <- df2
mrobj[[3]] <- df3
mrobj[[4]] <- df4
mrobj[[5]] <- df5
mrobj[[6]] <- df6
mrobj[[7]] <- df7
mrobj[[8]] <- df8
mrobj[[9]] <- df9

# Name slots
names(mrobj) <- c("oilyfish_TG", "oilyfish_CVD", "oilyfish_height",
                  "whiteXwholegrain_LDL", "whiteXwholegrain_CVD", "whiteXwholegrain_height",
                  "alcohol_ALT", "alcohol_cir", "alcohol_height") 

# Save
saveRDS(mrobj, "interim_data/merged_GI/MRobject_241113.rds")
