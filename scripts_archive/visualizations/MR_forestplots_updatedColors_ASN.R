# Goal: Forest plot code for MR-diet project - Updated for ASN presentation
# Author: KJS
# Date: 11/27/2024
# Date: 05/20/2025... sigh

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

# Load libraries
library(dplyr)
library(ggplot2)
library(MendelianRandomization)
library(RColorBrewer)
library(cowplot)

# Get updated colors for plots
paired_colors <- brewer.pal(12, "Paired")
print(paired_colors)
colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")

# Load IVW results
data <- read.csv("results/general/tables/MR_NonFiltered_and_Filtered_IVW_RandomEffects_Results_AllTraitPairs.csv", row.names = 1)

# Load WM results
#data_wm <- read.csv("results/general/tables/MR_NonFiltered_and_Filtered_WeightedMedian_AllTraitPairs_ProcessedManuscriptReady.csv", row.names = 1)

# Make label for traits
data$exp_out <- paste0(data$exposure, " \u2192 ", data$outcome) #\u2192 prints an arrow

# Factor exp_out 
desired_order <- c("Alcohol → Height",
                   "Alcohol → Liver Cirrhosis",
                   "Alcohol → Alanine Aminotransferase",
                   "White vs Whole Grain Bread → Height",
                   "White vs Whole Grain Bread → Cardiovascular Disease",
                   "White vs Whole Grain Bread → LDL Cholesterol",
                   "Oily fish → Height",
                   "Oily fish → Cardiovascular Disease",
                   "Oily Fish → Triglycerides")
data$exp_out <- factor(data$exp_out, levels = desired_order)

# Check
levels(data$exp_out)

# Now shorten them for plotting purposes
desired_order_short <- c("Alcohol → Height",
                   "Alcohol → Liver Cirrhosis",
                   "Alcohol → ALT",
                   "White vs Whole Grain Bread → Height",
                   "White vs Whole Grain Bread → CVD",
                   "White vs Whole Grain Bread → LDL-C",
                   "Oily fish → Height",
                   "Oily fish → CVD",
                   "Oily Fish → TG")
levels(data$exp_out) <- desired_order_short

# Edit for plot
# Rename test column with manuscript labels
data$test[data$test == "mr_core_ivw_randomeffects"] <- "Core"
data$test[data$test == "mr_ivw_randomeffects"] <- "All"
data$test[data$test == "mr_pve_ivw_randomeffects"] <- "PheWAS-PVE"
data$test[data$test == "mr_stg_ivw_randomeffects"] <- "Steiger"
data$test[data$test == "mr_ttest_ivw_randomeffects"] <- "PheWAS-TT"

# Add colors
#my_colors <- c("#d73027","#fc8d59", "#fee090", "#e0f3f8", "#91bfdb", "#4575b4")
my_colors <- c("#33A02C","#FB9A99", "#E31A1C", "#A6CEE3", "#1F78B4")

# Create the forest plot
forestplot <- ggplot(data, aes(x = estimate, y = exp_out, color = test)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3) + # Error bars
  scale_color_manual(values = my_colors) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "Category",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "bottom"
  )

forestplot

# Save
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
#ggsave(paste0(plots_dir,"/MR_IVW_forestplot_alltraits_allmethods_updatedColors.png"),
       plot = forestplot,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 10,
       height = 4,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Updated way to show the same data - Still IVW
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pos_df <- data %>% filter(outcome != "Height")
neg_df <- data %>% filter(outcome == "Height")

# Create the forest plot for positive controls
forestplot_pos <- ggplot(pos_df, aes(x = estimate, y = exp_out, color = test)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3) + # Error bars
  scale_color_manual(values = my_colors) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    legend.position = "bottom",
    legend.justification = "left"
  )

forestplot_pos

# Create the forest plot for negative controls
forestplot_neg <- ggplot(neg_df, aes(x = estimate, y = exp_out, color = test)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3) + # Error bars
  scale_color_manual(values = my_colors) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    legend.position = "bottom",
    legend.justification = "left"
  )

forestplot_neg

# Arrange with cowplot
# Arrange plots
combined_plot <- cowplot::plot_grid(
  forestplot_pos, forestplot_neg,
  ncol = 1,
  align = "v", # horizontal alignment (v=stacked)
  axis = "l", # align the left (y-axis) of the plots
  labels = c("A", "B"),  # Auto-label plots
  label_size = 16
)

combined_plot

# Save
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
ggsave(paste0(plots_dir,"/MR_IVW_forestplot_posnegFacet_allmethods_updatedColors.png"),
       plot = combined_plot,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 10,
       height = 10,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Separate by filtered and non-filtered GI - IVW
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pos_df <- data %>% filter(outcome != "Height")
pos_df_nf <- pos_df %>% filter(test == "All") #nf=no filter
pos_df_nf$color <- c("biomarker", "disease", "biomarker", "disease", "biomarker", "disease") # for illustrating biomarkers/disease 

neg_df <- data %>% filter(outcome == "Height")
neg_df_nf <- neg_df %>% filter(test == "All")
neg_df_nf$color <- c("disease", "disease","disease") # for illustrating biomarkers/disease 

# This will help illustrate the point about biomarkers + diseases
my_colors_nf <- c("#33A02C","#1F78B4")

# Create the forest plot for positive controls
forestplot_pos_nf <- ggplot(pos_df_nf, aes(x = estimate, y = exp_out, color = color)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3) + # Error bars
  scale_color_manual(values = my_colors_nf) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "none", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_pos_nf

# Create the forest plot for negative controls
forestplot_neg_nf <- ggplot(neg_df_nf, aes(x = estimate, y = exp_out, color = color)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3) + # Error bars
  scale_color_manual(values = my_colors_nf[2]) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "none", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_neg_nf

# Arrange with cowplot
# Arrange plots
combined_plot <- cowplot::plot_grid(
  forestplot_pos_nf, forestplot_neg_nf,
  ncol = 1,
  align = "v", # horizontal alignment (v=stacked)
  axis = "l", # align the left (y-axis) of the plots
  labels = c("A", "B"),  # Auto-label plots
  label_size = 16
)

combined_plot

# Save
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
ggsave(paste0(plots_dir,"/MR_IVW_forestplot_posNegFacet_nonfilteredGI_updatedColors.png"),
       plot = combined_plot,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 10,
       height = 10,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Separate by filtered and non-filtered GI - IVW - Adjustments for ASN Talk
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pos_df <- data %>% filter(outcome != "Height")
pos_df_nf <- pos_df %>% filter(test == "All") #nf=no filter
pos_df_nf$color <- c("biomarker", "disease", "biomarker", "disease", "biomarker", "disease") # for illustrating biomarkers/disease 

neg_df <- data %>% filter(outcome == "Height")
neg_df_nf <- neg_df %>% filter(test == "All")
neg_df_nf$color <- c("disease", "disease","disease") # for illustrating biomarkers/disease 

# Make the colors the same so I don't have to describe the difference between a biomarker 
# and a disease in my limited time.
my_colors_nf <- c("#1F78B4","#1F78B4") #"#1F78B4")

# Create the forest plot for positive controls
forestplot_pos_nf <- ggplot(pos_df_nf, aes(x = estimate, y = exp_out, color = color)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.4, linewidth = 1) + # Error bars
  scale_color_manual(values = my_colors_nf) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1.6) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "none", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_pos_nf

# Create the forest plot for negative controls
forestplot_neg_nf <- ggplot(neg_df_nf, aes(x = estimate, y = exp_out, color = color)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.4, linewidth = 1) + # Error bars
  scale_color_manual(values = my_colors_nf[2]) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1.6) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "none", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_neg_nf

# Arrange with cowplot
# Arrange plots
combined_plot <- cowplot::plot_grid(
  forestplot_pos_nf, forestplot_neg_nf,
  ncol = 1,
  align = "v", # horizontal alignment (v=stacked)
  axis = "l", # align the left (y-axis) of the plots
  labels = c("A", "B"),  # Auto-label plots
  label_size = 16
)

combined_plot

# Save
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
ggsave(paste0(plots_dir,"/MR_IVW_forestplot_posNegFacet_nonfilteredGI_updated4ASN.png"),
       plot = combined_plot,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 10,
       height = 10,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For the second half of the ASN talk, compare PheWAS-ttest to the unfiltered GI
# IVW
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pos_df <- data %>% filter(outcome != "Height")
pos_df_ttest_nf <- pos_df[pos_df$test == "PheWAS-TT" | pos_df$test == "All",] #nf=no filter

neg_df <- data %>% filter(outcome == "Height")
neg_df_ttest_nf <- neg_df[neg_df$test == "PheWAS-TT" | neg_df$test == "All",]

# Make the colors the same so I don't have to describe the difference between a biomarker 
# and a disease in my limited time.
my_colors_nf <- c("#1F78B4","#FF7F00") #"#1F78B4")

# Create the forest plot for positive controls
forestplot_pos_all_tt <- ggplot(pos_df_ttest_nf, aes(x = estimate, y = exp_out, color = test)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.4, linewidth = 1) + # Error bars
  scale_color_manual(values = my_colors_nf) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1.6) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "none", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_pos_all_tt

# Create the forest plot for negative controls
forestplot_neg_all_tt <- ggplot(neg_df_ttest_nf, aes(x = estimate, y = exp_out, color = test)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.4, linewidth = 1) + # Error bars
  scale_color_manual(values = my_colors_nf) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1.6) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "none", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_neg_all_tt

# Arrange with cowplot
# Arrange plots
combined_plot <- cowplot::plot_grid(
  forestplot_pos_all_tt, forestplot_neg_all_tt,
  ncol = 1,
  align = "v", # horizontal alignment (v=stacked)
  axis = "l", # align the left (y-axis) of the plots
  labels = c("A", "B"),  # Auto-label plots
  label_size = 16
)

combined_plot

# Save
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
ggsave(paste0(plots_dir,"/MR_IVW_forestplot_posNegFacet_nonfiltered_and_filtered_GI_ASN.png"),
       plot = combined_plot,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 10,
       height = 10,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Separate by filtered and non-filtered GI - IVW
# Filtered data here
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pos_df <- data %>% filter(outcome != "Height")
pos_df_f <- pos_df %>% filter(test != "All") #nf=no filter
pos_df_f$color <- c("biomarker", "disease", "biomarker", "disease", "biomarker", "disease") # for illustrating biomarkers/disease 

neg_df <- data %>% filter(outcome == "Height")
neg_df_f <- neg_df %>% filter(test != "All")
neg_df_f$color <- c("disease", "disease","disease") # for illustrating biomarkers/disease 

# Colors for the methods
# Don't use the colors for the IVW/WM plots above, "#33A02C" and "#1F78B4"
my_colors_f <- c("#6A3D9A", "#A6CEE3","#FB9A99","#FF7F00","#CAB2D6", "#FDBF6F", "#E31A1C")

# Create the forest plot for positive controls
forestplot_pos_f <- ggplot(pos_df_f, aes(x = estimate, y = exp_out, color = test)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3, linewidth = 1) + # Error bars
  scale_color_manual(values = my_colors_f[1:4]) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "bottom", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_pos_f

# Create the forest plot for negative controls
forestplot_neg_f <- ggplot(neg_df_f, aes(x = estimate, y = exp_out, color = test)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3, linewidth = 1) + # Error bars
  scale_color_manual(values = my_colors_f) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "bottom", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_neg_f

# Arrange with cowplot
# Arrange plots
combined_plot <- cowplot::plot_grid(
  forestplot_pos_f, forestplot_neg_f,
  ncol = 1,
  align = "v", # horizontal alignment (v=stacked)
  axis = "l", # align the left (y-axis) of the plots
  labels = c("A", "B"),  # Auto-label plots
  label_size = 16
)

combined_plot

# Save
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
ggsave(paste0(plots_dir,"/MR_IVW_forestplot_posNegFacet_filteredGI_updatedColors.png"),
       plot = combined_plot,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 10,
       height = 10,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Separate by filtered and non-filtered GI - IVW
# Filtered data here
# Use a subset of methods!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pos_df <- data %>% filter(outcome != "Height")
pos_df_f <- pos_df %>% filter(test != "All") #nf=no filter
pos_df_f$color <- c("biomarker", "disease", "biomarker", "disease", "biomarker", "disease") # for illustrating biomarkers/disease 
# Remove Steiger
pos_df_f2 <- pos_df_f[!(pos_df_f$test == "Steiger"),]

neg_df <- data %>% filter(outcome == "Height")
neg_df_f <- neg_df %>% filter(test != "All")
neg_df_f$color <- c("disease", "disease","disease") # for illustrating biomarkers/disease 
# Remove Steiger
neg_df_f2 <- neg_df_f[!(neg_df_f$test == "Steiger"),]

# Colors for the methods
# Don't use the colors for the IVW/WM plots above, "#33A02C" and "#1F78B4"
colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
my_colors_f <- c("#FB9A99", "#FF7F00", "#33A02C")
# Create the forest plot for positive controls
forestplot_pos_f2 <- ggplot(pos_df_f2, aes(x = estimate, y = exp_out, color = test)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3, linewidth = 1) + # Error bars
  scale_color_manual(values = my_colors_f[1:4]) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "bottom", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_pos_f2

# Create the forest plot for negative controls
forestplot_neg_f2 <- ggplot(neg_df_f2, aes(x = estimate, y = exp_out, color = test)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3, linewidth = 1) + # Error bars
  scale_color_manual(values = my_colors_f) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "bottom", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_neg_f2

# Arrange with cowplot
# Arrange plots
combined_plot <- cowplot::plot_grid(
  forestplot_pos_f2, forestplot_neg_f2,
  ncol = 1,
  align = "v", # horizontal alignment (v=stacked)
  axis = "l", # align the left (y-axis) of the plots
  labels = c("A", "B"),  # Auto-label plots
  label_size = 16
)

combined_plot

# Save
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
ggsave(paste0(plots_dir,"/MR_IVW_forestplot_posNegFacet_filteredGI_updatedColors_v2.png"),
       plot = combined_plot,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 10,
       height = 10,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Repeat for weighted median method
# WM!!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load WM results
data_wm <- read.csv("results/general/tables/MR_NonFiltered_and_Filtered_WeightedMedian_AllTraitPairs_ProcessedManuscriptReady.csv", row.names = 1)

# Make label for traits
data_wm$exp_out <- paste0(data_wm$exposure, " \u2192 ", data_wm$outcome) #\u2192 prints an arrow

# Factor exp_out 
desired_order <- c("Alcohol → Height",
                   "Alcohol → Liver Cirrhosis",
                   "Alcohol → Alanine Aminotransferase",
                   "White vs Whole Grain Bread → Height",
                   "White vs Whole Grain Bread → Cardiovascular Disease",
                   "White vs Whole Grain Bread → LDL Cholesterol",
                   "Oily fish → Height",
                   "Oily fish → Cardiovascular Disease",
                   "Oily Fish → Triglycerides")
data_wm$exp_out <- factor(data_wm$exp_out, levels = desired_order)

# Check
levels(data_wm$exp_out)

# Now shorten them for plotting purposes
desired_order_short <- c("Alcohol → Height",
                         "Alcohol → Liver Cirrhosis",
                         "Alcohol → ALT",
                         "White vs Whole Grain Bread → Height",
                         "White vs Whole Grain Bread → CVD",
                         "White vs Whole Grain Bread → LDL-C",
                         "Oily fish → Height",
                         "Oily fish → CVD",
                         "Oily Fish → TG")
levels(data_wm$exp_out) <- desired_order_short

# Edit for plot
# Rename test column with manuscript labels
data_wm$GI[data_wm$GI == "Core GI"] <- "Core"
data_wm$GI[data_wm$GI == "all"] <- "All"
data_wm$GI[data_wm$GI == "PheWAS PVE Filtered"] <- "PheWAS-PVE"
data_wm$GI[data_wm$GI == "Steiger Filtered"] <- "Steiger"
data_wm$GI[data_wm$GI == "PheWAS T-test Filtered"] <- "PheWAS-TT"

# Rename CI 
data_wm <- data_wm %>% rename(cilower = X95..CI.,
                              ciupper = X.,
                              estimate = Estimate)

# Add colors
#my_colors <- c("#d73027","#fc8d59", "#fee090", "#e0f3f8", "#91bfdb", "#4575b4")
my_colors <- c("#33A02C","#FB9A99", "#E31A1C", "#A6CEE3", "#1F78B4")

# Create the forest plot
forestplot <- ggplot(data_wm, aes(x = estimate, y = exp_out, color = GI)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3) + # Error bars
  scale_color_manual(values = my_colors) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "Category",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "bottom"
  )

forestplot

# Save
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
ggsave(paste0(plots_dir,"/MR_WM_forestplot_alltraits_allmethods_updatedColors.png"),
plot = forestplot,
device = NULL,
path = NULL,
scale = 1,
width = 10,
height = 4,
units = c("in"),
dpi = 300,
limitsize = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Separate by positive/negative controls - WM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pos_df <- data_wm %>% filter(outcome != "Height")
neg_df <- data_wm %>% filter(outcome == "Height")

# Create the forest plot for positive controls
forestplot_pos <- ggplot(pos_df, aes(x = estimate, y = exp_out, color = GI)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3) + # Error bars
  scale_color_manual(values = my_colors) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    legend.position = "bottom",
    legend.justification = "left"
  )

forestplot_pos

# Create the forest plot for negative controls
forestplot_neg <- ggplot(neg_df, aes(x = estimate, y = exp_out, color = GI)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3) + # Error bars
  scale_color_manual(values = my_colors) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    legend.position = "bottom",
    legend.justification = "left"
  )

forestplot_neg

# Arrange with cowplot
# Arrange plots
combined_plot <- cowplot::plot_grid(
  forestplot_pos, forestplot_neg,
  ncol = 1,
  align = "v", # horizontal alignment (v=stacked)
  axis = "l", # align the left (y-axis) of the plots
  labels = c("A", "B"),  # Auto-label plots
  label_size = 16
)

combined_plot

# Save
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
ggsave(paste0(plots_dir,"/MR_WM_forestplot_posnegFacet_allmethods_updatedColors.png"),
       plot = combined_plot,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 10,
       height = 10,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Separate by filtered and non-filtered GI - WM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pos_df <- data_wm %>% filter(outcome != "Height")
pos_df_nf <- pos_df %>% filter(GI == "All") #nf=no filter
pos_df_nf$color <- c("biomarker", "disease", "biomarker", "disease", "biomarker", "disease") # for illustrating biomarkers/disease 

neg_df <- data_wm %>% filter(outcome == "Height")
neg_df_nf <- neg_df %>% filter(GI == "All")
neg_df_nf$color <- c("disease", "disease","disease") # for illustrating biomarkers/disease 

# This will help illustrate the point about biomarkers + diseases
my_colors_nf <- c("#33A02C","#1F78B4")

# Create the forest plot for positive controls
forestplot_pos_nf <- ggplot(pos_df_nf, aes(x = estimate, y = exp_out, color = color)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3) + # Error bars
  scale_color_manual(values = my_colors_nf) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "none", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_pos_nf

# Create the forest plot for negative controls
forestplot_neg_nf <- ggplot(neg_df_nf, aes(x = estimate, y = exp_out, color = color)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3) + # Error bars
  scale_color_manual(values = my_colors_nf[2]) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "none", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_neg_nf

# Arrange with cowplot
# Arrange plots
combined_plot <- cowplot::plot_grid(
  forestplot_pos_nf, forestplot_neg_nf,
  ncol = 1,
  align = "v", # horizontal alignment (v=stacked)
  axis = "l", # align the left (y-axis) of the plots
  labels = c("C", "D"),  # Auto-label plots
  label_size = 16
)

combined_plot

# Save
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
ggsave(paste0(plots_dir,"/MR_WM_forestplot_posNegFacet_nonfilteredGI_updatedColors.png"),
       plot = combined_plot,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 10,
       height = 10,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Separate by filtered and non-filtered GI - WM - Updated for ASN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pos_df <- data_wm %>% filter(outcome != "Height")
pos_df_nf <- pos_df %>% filter(GI == "All") #nf=no filter
pos_df_nf$color <- c("biomarker", "disease", "biomarker", "disease", "biomarker", "disease") # for illustrating biomarkers/disease 

neg_df <- data_wm %>% filter(outcome == "Height")
neg_df_nf <- neg_df %>% filter(GI == "All")
neg_df_nf$color <- c("disease", "disease","disease") # for illustrating biomarkers/disease 

# Updated for ASN - simplify the colors/info the colors are conveying
my_colors_nf <- c("#1F78B4","#1F78B4")

# Create the forest plot for positive controls
forestplot_pos_nf <- ggplot(pos_df_nf, aes(x = estimate, y = exp_out, color = color)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.4, linewidth = 1) + # Error bars
  scale_color_manual(values = my_colors_nf) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1.6) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "none", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_pos_nf

# Create the forest plot for negative controls
forestplot_neg_nf <- ggplot(neg_df_nf, aes(x = estimate, y = exp_out, color = color)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.4, linewidth = 1) + # Error bars
  scale_color_manual(values = my_colors_nf[2]) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1.6) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "none", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_neg_nf

# Arrange with cowplot
# Arrange plots
combined_plot <- cowplot::plot_grid(
  forestplot_pos_nf, forestplot_neg_nf,
  ncol = 1,
  align = "v", # horizontal alignment (v=stacked)
  axis = "l", # align the left (y-axis) of the plots
  labels = c("C", "D"),  # Auto-label plots
  label_size = 16
)

combined_plot

# Save
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
ggsave(paste0(plots_dir,"/MR_WM_forestplot_posNegFacet_nonfilteredGI_updated4ASN.png"),
       plot = combined_plot,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 10,
       height = 10,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For the second half of the ASN talk, compare PheWAS-ttest to the unfiltered GI
# WM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pos_df <- data_wm %>% filter(outcome != "Height")
pos_df_ttest_nf <- pos_df[pos_df$GI == "PheWAS-TT" | pos_df$GI == "All",] #nf=no filter

neg_df <- data_wm %>% filter(outcome == "Height")
neg_df_ttest_nf <- neg_df[neg_df$GI == "PheWAS-TT" | neg_df$GI == "All",]

# Make the colors the same so I don't have to describe the difference between a biomarker 
# and a disease in my limited time.
my_colors_nf <- c("#1F78B4","#FF7F00") #"#1F78B4")

# Create the forest plot for positive controls
forestplot_pos_all_tt <- ggplot(pos_df_ttest_nf, aes(x = estimate, y = exp_out, color = GI)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.4, linewidth = 1) + # Error bars
  scale_color_manual(values = my_colors_nf) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1.6) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "none", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_pos_all_tt

# Create the forest plot for negative controls
forestplot_neg_all_tt <- ggplot(neg_df_ttest_nf, aes(x = estimate, y = exp_out, color = GI)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.4, linewidth = 1) + # Error bars
  scale_color_manual(values = my_colors_nf) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1.6) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "none", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_neg_all_tt

# Arrange with cowplot
# Arrange plots
combined_plot <- cowplot::plot_grid(
  forestplot_pos_all_tt, forestplot_neg_all_tt,
  ncol = 1,
  align = "v", # horizontal alignment (v=stacked)
  axis = "l", # align the left (y-axis) of the plots
  labels = c("C", "D"),  # Auto-label plots
  label_size = 16
)

combined_plot

# Save
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
ggsave(paste0(plots_dir,"/MR_WM_forestplot_posNegFacet_nonfiltered_and_filtered_GI_ASN.png"),
       plot = combined_plot,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 10,
       height = 10,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Separate by filtered and non-filtered GI - WM
# Filtered data here
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pos_df <- data_wm %>% filter(outcome != "Height")
pos_df_f <- pos_df %>% filter(GI != "All") #nf=no filter
pos_df_f$color <- c("biomarker", "disease", "biomarker", "disease", "biomarker", "disease") # for illustrating biomarkers/disease 

neg_df <- data_wm %>% filter(outcome == "Height")
neg_df_f <- neg_df %>% filter(GI != "All")
neg_df_f$color <- c("disease", "disease","disease") # for illustrating biomarkers/disease 

# Colors for the methods
# Don't use the colors for the IVW/WM plots above, "#33A02C" and "#1F78B4"
my_colors_f <- c("#6A3D9A", "#A6CEE3","#FB9A99","#FF7F00","#CAB2D6", "#FDBF6F", "#E31A1C")

# Create the forest plot for positive controls
forestplot_pos_f <- ggplot(pos_df_f, aes(x = estimate, y = exp_out, color = GI)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3, linewidth = 1) + # Error bars
  scale_color_manual(values = my_colors_f[1:4]) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "bottom", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_pos_f

# Create the forest plot for negative controls
forestplot_neg_f <- ggplot(neg_df_f, aes(x = estimate, y = exp_out, color = GI)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3, linewidth = 1) + # Error bars
  scale_color_manual(values = my_colors_f) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "bottom", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_neg_f

# Arrange with cowplot
# Arrange plots
combined_plot <- cowplot::plot_grid(
  forestplot_pos_f, forestplot_neg_f,
  ncol = 1,
  align = "v", # horizontal alignment (v=stacked)
  axis = "l", # align the left (y-axis) of the plots
  labels = c("C", "D"),  # Auto-label plots
  label_size = 16
)

combined_plot

# Save
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
ggsave(paste0(plots_dir,"/MR_WM_forestplot_posNegFacet_filteredGI_updatedColors.png"),
       plot = combined_plot,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 10,
       height = 10,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Separate by filtered and non-filtered GI - WM
# Filtered data here
# Use a subset of methods!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pos_df <- data_wm %>% filter(outcome != "Height")
pos_df_f <- pos_df %>% filter(GI != "All") #nf=no filter
pos_df_f$color <- c("biomarker", "disease", "biomarker", "disease", "biomarker", "disease") # for illustrating biomarkers/disease 
# Remove Steiger
pos_df_f2 <- pos_df_f[!(pos_df_f$GI == "Steiger"),]

neg_df <- data_wm %>% filter(outcome == "Height")
neg_df_f <- neg_df %>% filter(GI != "All")
neg_df_f$color <- c("disease", "disease","disease") # for illustrating biomarkers/disease 
# Remove Steiger
neg_df_f2 <- neg_df_f[!(neg_df_f$GI == "Steiger"),]

# Colors for the methods
# Don't use the colors for the IVW/WM plots above, "#33A02C" and "#1F78B4"
colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
my_colors_f <- c("#FB9A99", "#FF7F00", "#33A02C")
# Create the forest plot for positive controls
forestplot_pos_f2 <- ggplot(pos_df_f2, aes(x = estimate, y = exp_out, color = GI)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3, linewidth = 1) + # Error bars
  scale_color_manual(values = my_colors_f[1:4]) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "bottom", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_pos_f2

# Create the forest plot for negative controls
forestplot_neg_f2 <- ggplot(neg_df_f2, aes(x = estimate, y = exp_out, color = GI)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3, linewidth = 1) + # Error bars
  scale_color_manual(values = my_colors_f) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"), #14 looks good in print
    legend.text = element_text(size = 14),
    legend.position = "bottom", 
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.5)
    #legend.justification = "left"
  )

forestplot_neg_f2

# Arrange with cowplot
# Arrange plots
combined_plot <- cowplot::plot_grid(
  forestplot_pos_f2, forestplot_neg_f2,
  ncol = 1,
  align = "v", # horizontal alignment (v=stacked)
  axis = "l", # align the left (y-axis) of the plots
  labels = c("C", "D"),  # Auto-label plots
  label_size = 16
)

combined_plot

# Save
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
ggsave(paste0(plots_dir,"/MR_WM_forestplot_posNegFacet_filteredGI_updatedColors_v2.png"),
       plot = combined_plot,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 10,
       height = 10,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)
