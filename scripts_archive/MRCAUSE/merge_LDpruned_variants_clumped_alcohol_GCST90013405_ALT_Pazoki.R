# Combine the LD pruned individual chromosome data into 1 file

#Merge the files together
library(data.table)

# USER INPUT
exposure <- "alcohol"
outcome <- "GCST90013405_ALT_Pazoki" 

# GENERIC CODE BELOW
#Read files with full path
files <- list.files(paste0("/pl/active/colelab/users/kjames/enviroMR/LD_1KG_pruned/",exposure,"_",outcome), full.names = TRUE)

#Combine all data
data <- rbindlist(lapply(files, fread), fill = TRUE)

#Save
write.csv(data, paste0("/pl/active/colelab/users/kjames/enviroMR/MR_CAUSE/pruned_data/LD_1KG_pruned_",exposure,"_",outcome,".csv"))

