# Author: Kristen Sutton
# Date: September 30, 2024
# Goal: With all filtering methods complete, look at summary data together

# Libraries 
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(VennDiagram)

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

# Path to save plots
vd_path <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/venndiagram"

# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

# Helper function to save Venn diagram
save_venn <- function(x, name, ...){
  library(VennDiagram)
  venn.plot <- venn.diagram(x,
               filename = paste0(vd_path, "/", name, "_venn.png"),
               output = TRUE,
               # Circlues
               lwd = 2, 
               lty = 'blank',
               fill = c("#d73027", "#fee090", "#4575b4"),
               # Numbers
               cex = 2,
               fontfamily = "sans",
               # Set names
               cat.fontface = "bold",
               cat.fontfamily = "sans",
               cat.default.pos = "outer",
               cat.pos = c(-27, 27, 135)
               )
}

# MRobj
mrobj <- readRDS("interim_data/merged_GI/MRobject_steiger_PheWAStt_PheWASPVE_241120.rds")

# Make mrstg column
# Using `lapply()` to add a new stgfilt_keep column
mrobj <- lapply(mrobj, function(df) {
  df$stgfilt_keep <- ifelse(df$steiger_dir == "TRUE" & df$steiger_pval <= 0.05, TRUE, FALSE )
  return(df)
})

# Rename for backup
my_list = mrobj

# Columns to use in sets
columns_to_use <- c("PheWAS_PVE_keep", "PheWAS_ttest_keep", "stgfilt_keep")

# Function to create a Venn diagram for a single data frame
create_venn_diagram <- function(df, name) {
  # Create a list of sets based on TRUE values in the specified columns
  sets <- list(
    PheWAS_PVE = which(df[[columns_to_use[1]]] == TRUE), #PVE
    PheWAS_TTest = which(df[[columns_to_use[2]]] == TRUE), #PheWASClust
    Steiger = which(df[[columns_to_use[3]]] == TRUE) #Steiger
  )
  
  # Print a message to confirm creation
  print(paste("Venn diagram for", name, "created."))
  
  # Display venn diagram in R
  display_venn(sets)
  # Save plots
  save_venn(sets, names(my_list)[i])
}

for (i in seq_along(my_list)) {
  create_venn_diagram(my_list[[i]], names(my_list)[i])
}

################
################
################
# Run again without names on the circles (for publication - I want to put the name where I want in adobe)

# Path to save plots
#vd_path <- "/pl/active/colelab/users/kjames/enviroMR/results/julie_collab/venn_diagram"

# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

# Helper function to save Venn diagram
save_venn <- function(x, name, ...){
  library(VennDiagram)
  venn.plot <- venn.diagram(x,
                            filename = paste0(vd_path, "/", name, "_venn_nonames.png"), #added nonames extension
                            output = TRUE,
                            # Circlues
                            lwd = 2, 
                            lty = 'blank',
                            fill = c("#d73027", "#fee090", "#4575b4"),
                            # Numbers
                            cex = 2,
                            fontfamily = "sans",
                            # Set names
                            cat.cex = 0,
                            #cat.fontface = "bold",
                            #cat.fontfamily = "sans",
                            cat.default.pos = "outer",
                            cat.pos = c(-27, 27, 135)
  )
}




# MRobj
mrobj <- readRDS("interim_data/merged_GI/MRobject_steiger_PheWAStt_PheWASPVE_241120.rds")

# Make mrstg column
# Using `lapply()` to add a new stgfilt_keep column
mrobj <- lapply(mrobj, function(df) {
  df$stgfilt_keep <- ifelse(df$steiger_dir == "TRUE" & df$steiger_pval <= 0.05, TRUE, FALSE )
  return(df)
})

# Rename for backup
my_list = mrobj

# Columns to use in sets
columns_to_use <- c("PheWAS_PVE_keep", "PheWAS_ttest_keep", "stgfilt_keep")

# Function to create a Venn diagram for a single data frame
create_venn_diagram <- function(df, name) {
  # Create a list of sets based on TRUE values in the specified columns
  sets <- list(
    PheWAS_PVE = which(df[[columns_to_use[1]]] == TRUE), #PVE
    PheWAS_Ttest = which(df[[columns_to_use[2]]] == TRUE), #PheWASClust
    Steiger = which(df[[columns_to_use[3]]] == TRUE) #Steiger
  )
  
  # Print a message to confirm creation
  print(paste("Venn diagram for", name, "created."))
  
  # Display venn diagram in R
  display_venn(sets)
  # Save plots
  save_venn(sets, names(my_list)[i])
}

for (i in seq_along(my_list)) {
  create_venn_diagram(my_list[[i]], names(my_list)[i])
}

