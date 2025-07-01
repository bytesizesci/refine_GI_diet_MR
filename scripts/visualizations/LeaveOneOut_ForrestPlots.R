# Goal: Perform leave one out analyses
# Author: KJS
# Date: 12/18/24

# Libraries 
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(MendelianRandomization)

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

# Load harmonized and standardized genetic instruments
mrobj <- readRDS("interim_data/merged_GI/MRobject_steiger_PheWAStt_PheWASPVE_241120.rds")

# Create a data frame with naming conventions for loop
df_names <- data.frame(EXP_pheno = c("oilyfish", "oilyfish", "oilyfish", "whiteXwholegrain", "whiteXwholegrain",  "whiteXwholegrain", "alcohol", "alcohol","alcohol"),
                       OUT_pheno = c("GCST90239664_TG_Graham", "GCST90132314_CVD_Aragam", "GCST006901_Height_Yengo_2018", "GCST90239658_LDL_Graham","GCST90132314_CVD_Aragam", "GCST006901_Height_Yengo_2018", "GCST90013405_ALT_Pazoki", "GCST90319877_cir_Ghouse", "GCST006901_Height_Yengo_2018"),
                       EXP_names = c("Oily Fish", "Oily fish", "Oily fish", "White vs Whole Grain Bread", "White vs Whole Grain Bread", "White vs Whole Grain Bread","Alcohol", "Alcohol", "Alcohol"),
                       OUT_names = c("Triglycerides", "Cardiovascular Disease", "Height", "LDL Cholesterol", "Cardiovascular Disease", "Height", "Alanine Aminotransferase", "Liver Cirrhosis", "Height"))

# Perform MR 
for (i in seq_along(mrobj)){
  # Get data frame
  mrdat <- mrobj[[i]]
  # Get trait names
  # Get names
  EXP_pheno <- df_names[i,1]
  OUT_pheno <- df_names[i,2]
  EXP_name <- df_names[i,3]
  OUT_name <- df_names[i,4]
  
  # For no-filtering instrument, use the whole data frame. Comment out filtering
  #mrdat_stgfilt <- mrdat[mrdat$steiger_dir == "TRUE" | mrdat$steiger_pval >= 0.05,]
  # Perform MR
  mrInput <- MendelianRandomization:: mr_input(bx = mrdat$beta.exposure,
                                               bxse = mrdat$se.exposure,
                                               by = mrdat$beta.outcome,
                                               byse = mrdat$se.outcome)
  mr_LOO_res <- MendelianRandomization::mr_loo(mrInput)
  mr_LOO_res
  
  # Custom plot for external facing meetings
  custom_LOO_forest <-  mr_LOO_res + #ggplot(mrdat, aes(x = beta.exposure, y = beta.outcome)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      #legend.position="none", 
      axis.title.y=ggplot2::element_text(size=14),
      axis.text.y=ggplot2::element_text(size=14), 
      axis.ticks.y=ggplot2::element_line(size=0),
      axis.title.x=ggplot2::element_text(size=14),
      axis.text.x=ggplot2::element_text(size=14),
      axis.ticks.x=ggplot2::element_line(size=0),
      panel.border = element_rect(color = "black", 
                                  fill = NA, 
                                  size = 3)) +
    #labs(x = paste0("Genetic Effect on ", EXP_name),
    #     y = paste0( "Genetic Effect on ", OUT_name)) +
    ggtitle(paste0(EXP_name, " on ",  OUT_name)) 
    #scale_color_manual(values = c("#3E4A89FF", "#6DCD59FF"))
  
    custom_LOO_forest
  
    # Save the plots
    plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots"
  
    ggsave(paste0(plots_dir,"/NonfilteredGI_MR_forrestplot_",EXP_pheno, "_", OUT_pheno,".png"),
         plot = custom_LOO_forest,
         device = NULL,
         path = NULL,
         scale = 1,
         width = 6,
         height = 12,
         units = c("in"),
         dpi = 300,
         limitsize = TRUE)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Custom forest plots with variant names
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Try again but we want to retain variant names...

# Perform MR 
for (i in seq_along(mrobj)){
  # Get data frame
  mrdat <- mrobj[[i]]
  # Get trait names
  # Get names
  EXP_pheno <- df_names[i,1]
  OUT_pheno <- df_names[i,2]
  EXP_name <- df_names[i,3]
  OUT_name <- df_names[i,4]

  # Uncomment for testing
  #mrdat <- mrobj[[4]]

  # Get mr object
  mrInput <- MendelianRandomization:: mr_input(bx = mrdat$beta.exposure,
                                               bxse = mrdat$se.exposure,
                                               by = mrdat$beta.outcome,
                                               byse = mrdat$se.outcome)
  # Run LOO
  mr_LOO_res <- MendelianRandomization::mr_loo(mrInput)
  # Get LOO data
  mr_LOO_plotdata <- mr_LOO_res$data # mr_LOO_plotdata$snps contains a row for every SNP and the last row is the ivw estimate
  
  # We want to rename all the rows except the last one
  # Note, some slots have proxies and some don't
  # And, some slots have SNP.x, SNP, and/or SNP_proxy 
  # Use a custom function to account for this behavior
  SNP_proxy <- "SNP_proxy" # slot 4
  SNP.x <- "SNP.x" # slot 3,6,8,9
  
  if (SNP_proxy %in% colnames(mrdat)) {
    # Case A: SNP_proxy exists
    print(paste("Slot", i, "- SNP_proxy column exists"))
    # Combine SNP and SNP_proxy (complementary)
    mrdat2 <- mrdat %>% mutate(combined_variant = coalesce(SNP, SNP_proxy))
    variant_names <- c(mrdat2$combined_variant, "IVW estimate")
    
  } else if ("SNP.x" %in% colnames(mrdat)) {
    # Case B: SNP.x exists
    print(paste("Slot", i, "- SNP.x column exists"))
    # Combine SNP and SNP.x (complementary)
    mrdat2 <- mrdat %>% mutate(combined_variant = coalesce(SNP.x, SNP))
    variant_names <- c(mrdat2$combined_variant, "IVW estimate")
  
  } else {
    # Case C: Neither SNP_proxy nor SNP.x exists
    print(paste("Slot", i, "- Neither SNP_proxy nor SNP.x column exists"))
    # Use SNP alone
    variant_names <- c(mrdat$SNP, "IVW estimate")
  }
  
  #mrdat2 <- mrdat %>% mutate(combined_variant = coalesce(SNP, SNP_proxy))
  #variant_names <- c(mrdat2$combined_variant, "IVW estimate") #add ivw estimate
  
  # Now set names
  mr_LOO_plotdata$snps <- variant_names
  
  # Make ivw estimate the last row of the group
  # Specify the group you want to appear last 
  mr_LOO_plotdata$snps <- forcats::fct_relevel(mr_LOO_plotdata$snps, "IVW estimate", after = Inf)
  
  # Check the levels
  levels(mr_LOO_plotdata$snps)
  
  # Make custom plot
  custom_LOO_forest2 <- ggplot(mr_LOO_plotdata, aes(x = estimates, y = snps, xmin=CI_lower, xmax = CI_upper)) +
    geom_point() +
    geom_linerange() +
    scale_y_discrete(limits = rev(mr_LOO_plotdata$snps)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Variant", y = "Effect Size", title = paste0(EXP_name, " on ",  OUT_name)) +
    theme_bw() +
    theme(
      #legend.position="none", 
      axis.title.y=ggplot2::element_text(size=14),
      axis.text.y=ggplot2::element_text(size=14), 
      axis.ticks.y=ggplot2::element_line(size=0),
      axis.title.x=ggplot2::element_text(size=14),
      axis.text.x=ggplot2::element_text(size=14),
      axis.ticks.x=ggplot2::element_line(size=0),
      panel.border = element_rect(color = "black", 
                                  fill = NA, 
                                  size = 3))

  custom_LOO_forest2
  
  # Save the plots
  plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots"
  
  ggsave(paste0(plots_dir,"/forestplots/NonfilteredGI_MR_custom_forrestplot_",EXP_pheno, "_", OUT_pheno,".png"),
         plot = custom_LOO_forest2,
         device = NULL,
         path = NULL,
         scale = 1,
         width = 10,
         height = 12,
         units = c("in"),
         dpi = 300,
         limitsize = TRUE)
}
 