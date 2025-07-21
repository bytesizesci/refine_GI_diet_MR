#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=05:00:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --job-name=LDprune_oilyfish_height
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristen.sutton@cuanschutz.edu
#SBATCH --output=/pl/active/colelab/users/kjames/enviroMR/config2/LDprune.%j.out
#SBATCH --error=/pl/active/colelab/users/kjames/enviroMR/config2/LDprune.%j.err
#SBATCH --mem=30G 

# Be in "/pl/active/colelab/users/kjames/enviroMR" when you run this sbatch
cd /pl/active/colelab/users/kjames/enviroMR
 
# Step 0
# Exporting the apptainer tmp directories so that
# we do not fill /tmp
export ALPINE_SCRATCH=/gpfs/alpine1/scratch/$USER
export APPTAINER_TMPDIR=$ALPINE_SCRATCH/singularity/tmp
export APPTAINER_CACHEDIR=$ALPINE_SCRATCH/singularity/cache
mkdir -pv $APPTAINER_CACHEDIR $APPTAINER_TMPD
 
# Step 1 – to use with R interactively (install packages or test code line by line)
apptainer shell --fakeroot --overlay /projects/$USER/.rstudioserver/rstudio-server-4.2.2_overlay.img /curc/sw/containers/open_ondemand/rstudio-server-4.2.2.sif
 
# Uncomment the next line and delete "Apptainer>" if you wanted to use R interactively. Saved here as a reminder of other ways to use the apptainer. 
#Apptainer> R
# make sure to “exit” “exit” when you are finished
 
# Step 2 – launch an R script
# Replace test_script.R with your script!
apptainer exec --bind /projects/$USER,/scratch/alpine/$USER --overlay /projects/$USER/.rstudioserver/rstudio-server-4.2.2_overlay.img:ro /curc/sw/containers/open_ondemand/rstudio-server-4.2.2.sif Rscript scripts/MRCAUSE_pruneWorkFlow_clumped_oilyfish_GIANT_height_Yengo.R
