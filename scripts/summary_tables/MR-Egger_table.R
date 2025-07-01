# Goal: Grab MR-Egger Pleiotropy results
# Date: 06/12/25
# Author: KJS

# We previously did all the calculations but didn't include MR-Egger results in the 
# final tables. We now need to add MR-Egger. 

#Set working directory
setwd("/pl/active/colelab/users/kjames/refinedMR/")

combined_table_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/Nonfiltered_GI/tables"

# Load results generated from scripts/Nonfiltered_GI/MR_Nonfiltered_GI.R (not from the sensitivity analysis) 
df1 <- read.csv("/pl/active/colelab/users/kjames/refinedMR/results/Nonfiltered_GI/tables/MR_NonFiltered_and_Filtered_Results_Formattedoilyfish_GCST90239664_TG_Graham.csv")
df2 <- read.csv("/pl/active/colelab/users/kjames/refinedMR/results/Nonfiltered_GI/tables/MR_NonFiltered_and_Filtered_Results_Formattedoilyfish_GCST90132314_CVD_Aragam.csv")
df3 <- read.csv("/pl/active/colelab/users/kjames/refinedMR/results/Nonfiltered_GI/tables/MR_NonFiltered_and_Filtered_Results_Formattedoilyfish_GCST006901_Height_Yengo_2018.csv")

df4 <- read.csv("/pl/active/colelab/users/kjames/refinedMR/results/Nonfiltered_GI/tables/MR_NonFiltered_and_Filtered_Results_FormattedwhiteXwholegrain_GCST90239658_LDL_Graham.csv")
df5 <- read.csv("/pl/active/colelab/users/kjames/refinedMR/results/Nonfiltered_GI/tables/MR_NonFiltered_and_Filtered_Results_FormattedwhiteXwholegrain_GCST90132314_CVD_Aragam.csv")
df6 <- read.csv("/pl/active/colelab/users/kjames/refinedMR/results/Nonfiltered_GI/tables/MR_NonFiltered_and_Filtered_Results_FormattedwhiteXwholegrain_GCST006901_Height_Yengo_2018.csv")

df7 <- read.csv("/pl/active/colelab/users/kjames/refinedMR/results/Nonfiltered_GI/tables/MR_NonFiltered_and_Filtered_Results_Formattedalcohol_GCST90013405_ALT_Pazoki.csv")
df8 <- read.csv("/pl/active/colelab/users/kjames/refinedMR/results/Nonfiltered_GI/tables/MR_NonFiltered_and_Filtered_Results_Formattedalcohol_GCST90319877_cir_Ghouse.csv")
df9 <- read.csv("/pl/active/colelab/users/kjames/refinedMR/results/Nonfiltered_GI/tables/MR_NonFiltered_and_Filtered_Results_Formattedalcohol_GCST006901_Height_Yengo_2018.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# df1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Grab MR-Egger and MR-Egger intercept rows
# Manually look at rows to select
# Rows 36 to 45
df1_mre_top <- df1[36:40,]
df1_mre_bot <- df1[41:45,]
df1_mre_all <- cbind(df1_mre_top, df1_mre_bot)
df1_mre_all <- df1_mre_all[,c(2:9, 12:16)]
# Define column names
newnames <- c("Method","Estimate","Std.Error","X95.CI.low","X95.CI.high","P.value","GI","N_snps",
              "Intercept.Estimate","Intercept.Std.Error","Intercept.X95.CI.low","Intercept.X95.CI.high","Intercept.P.value")
colnames(df1_mre_all) <- newnames
# Make a pheno column to play it safe
df1_mre_all$diet_out <- "oilyfish_TG"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# df2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Grab MR-Egger and MR-Egger intercept rows
# Manually look at rows to select
# Rows 36 to 45
df2_mre_top <- df2[36:40,]
df2_mre_bot <- df2[41:45,]
df2_mre_all <- cbind(df2_mre_top, df2_mre_bot)
df2_mre_all <- df2_mre_all[,c(2:9, 12:16)]
# Define column names
newnames <- c("Method","Estimate","Std.Error","X95.CI.low","X95.CI.high","P.value","GI","N_snps",
              "Intercept.Estimate","Intercept.Std.Error","Intercept.X95.CI.low","Intercept.X95.CI.high","Intercept.P.value")
colnames(df2_mre_all) <- newnames
# Make a pheno column to play it safe
df2_mre_all$diet_out <- "oilyfish_CVD"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# df3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Grab MR-Egger and MR-Egger intercept rows
# Manually look at rows to select
# Rows 36 to 45
df3_mre_top <- df3[36:40,]
df3_mre_bot <- df3[41:45,]
df3_mre_all <- cbind(df3_mre_top, df3_mre_bot)
df3_mre_all <- df3_mre_all[,c(2:9, 12:16)]
# Define column names
newnames <- c("Method","Estimate","Std.Error","X95.CI.low","X95.CI.high","P.value","GI","N_snps",
              "Intercept.Estimate","Intercept.Std.Error","Intercept.X95.CI.low","Intercept.X95.CI.high","Intercept.P.value")
colnames(df3_mre_all) <- newnames
# Make a pheno column to play it safe
df3_mre_all$diet_out <- "oilyfish_height"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# df4
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Grab MR-Egger and MR-Egger intercept rows
# Manually look at rows to select
# Rows 36 to 45
df4_mre_top <- df4[36:40,]
df4_mre_bot <- df4[41:45,]
df4_mre_all <- cbind(df4_mre_top, df4_mre_bot)
df4_mre_all <- df4_mre_all[,c(2:9, 12:16)]
# Define column names
newnames <- c("Method","Estimate","Std.Error","X95.CI.low","X95.CI.high","P.value","GI","N_snps",
              "Intercept.Estimate","Intercept.Std.Error","Intercept.X95.CI.low","Intercept.X95.CI.high","Intercept.P.value")
colnames(df4_mre_all) <- newnames
# Make a pheno column to play it safe
df4_mre_all$diet_out <- "whiteXwholegrain_LDL"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# df5
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Grab MR-Egger and MR-Egger intercept rows
# Manually look at rows to select
# Rows 36 to 45
df5_mre_top <- df5[36:40,]
df5_mre_bot <- df5[41:45,]
df5_mre_all <- cbind(df5_mre_top, df5_mre_bot)
df5_mre_all <- df5_mre_all[,c(2:9, 12:16)]
# Define column names
newnames <- c("Method","Estimate","Std.Error","X95.CI.low","X95.CI.high","P.value","GI","N_snps",
              "Intercept.Estimate","Intercept.Std.Error","Intercept.X95.CI.low","Intercept.X95.CI.high","Intercept.P.value")
colnames(df5_mre_all) <- newnames
# Make a pheno column to play it safe
df5_mre_all$diet_out <- "whiteXwholegrain_CVD"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# df6
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Grab MR-Egger and MR-Egger intercept rows
# Manually look at rows to select
# Rows 36 to 45
df6_mre_top <- df6[36:40,]
df6_mre_bot <- df6[41:45,]
df6_mre_all <- cbind(df6_mre_top, df6_mre_bot)
df6_mre_all <- df6_mre_all[,c(2:9, 12:16)]
# Define column names
newnames <- c("Method","Estimate","Std.Error","X95.CI.low","X95.CI.high","P.value","GI","N_snps",
              "Intercept.Estimate","Intercept.Std.Error","Intercept.X95.CI.low","Intercept.X95.CI.high","Intercept.P.value")
colnames(df6_mre_all) <- newnames
# Make a pheno column to play it safe
df6_mre_all$diet_out <- "whiteXwholegrain_height"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# df7
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Grab MR-Egger and MR-Egger intercept rows
# Manually look at rows to select
# Rows 36 to 45
df7_mre_top <- df7[36:40,]
df7_mre_bot <- df7[41:45,]
df7_mre_all <- cbind(df7_mre_top, df7_mre_bot)
df7_mre_all <- df7_mre_all[,c(2:9, 12:16)]
# Define column names
newnames <- c("Method","Estimate","Std.Error","X95.CI.low","X95.CI.high","P.value","GI","N_snps",
              "Intercept.Estimate","Intercept.Std.Error","Intercept.X95.CI.low","Intercept.X95.CI.high","Intercept.P.value")
colnames(df7_mre_all) <- newnames
# Make a pheno column to play it safe
df7_mre_all$diet_out <- "alcohol_ALT"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# df8
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Grab MR-Egger and MR-Egger intercept rows
# Manually look at rows to select
# Rows 36 to 45
df8_mre_top <- df8[36:40,]
df8_mre_bot <- df8[41:45,]
df8_mre_all <- cbind(df8_mre_top, df8_mre_bot)
df8_mre_all <- df8_mre_all[,c(2:9, 12:16)]
# Define column names
newnames <- c("Method","Estimate","Std.Error","X95.CI.low","X95.CI.high","P.value","GI","N_snps",
              "Intercept.Estimate","Intercept.Std.Error","Intercept.X95.CI.low","Intercept.X95.CI.high","Intercept.P.value")
colnames(df8_mre_all) <- newnames
# Make a pheno column to play it safe
df8_mre_all$diet_out <- "alcohol_cirrhosis"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# df9
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Grab MR-Egger and MR-Egger intercept rows
# Manually look at rows to select
# Rows 36 to 45
df9_mre_top <- df9[36:40,]
df9_mre_bot <- df9[41:45,]
df9_mre_all <- cbind(df9_mre_top, df9_mre_bot)
df9_mre_all <- df9_mre_all[,c(2:9, 12:16)]
# Define column names
newnames <- c("Method","Estimate","Std.Error","X95.CI.low","X95.CI.high","P.value","GI","N_snps",
              "Intercept.Estimate","Intercept.Std.Error","Intercept.X95.CI.low","Intercept.X95.CI.high","Intercept.P.value")
colnames(df9_mre_all) <- newnames
# Make a pheno column to play it safe
df9_mre_all$diet_out <- "alcohol_height"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combine in one large table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
egger_res <- rbind(df1_mre_all, 
                   df2_mre_all,
                   df3_mre_all,
                   df4_mre_all,
                   df5_mre_all,
                   df6_mre_all,
                   df7_mre_all,
                   df8_mre_all,
                   df9_mre_all)

# Save
write.csv(egger_res, "results/general/tables/MR_NonFiltered_and_Filtered_MREgger_Alltraits.csv")
