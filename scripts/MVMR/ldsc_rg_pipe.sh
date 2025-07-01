#!/bin/sh
#
#SBATCH --nodes=1                   #NUMBER OF REQUESTED NODES
#SBATCH --ntasks=1                  #NUMBER OF REQUESTED CORES
#SBATCH --partition=amilan	   #PARTITION
#SBATCH --account=amc-general	   #ACCOUNT FOR ANSCHUTZ
#SBATCH --time=4:00:00             #WALL TIME
#SBATCH --mem=2gb                  #SPECIFY MEMORY IN MB
#SBATCH --array=1-10              #SPLIT THIS JOB INTO ARRAY FOR EACH TRAIT, but do not include the last file (82 instead of 83)
#SBATCH --output=/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/err_out/calc_LDSC_%A.%a.out
#SBATCH --error=/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/err_out/calc_LDSC_%A.%a.err


#load ldsc module
module use --append /pl/active/colelab/common/software/lmod-files
module load ldsc

#set directories
in_dir=/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PVE_descending_input

out_dir=/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/results

ref_dir=/pl/active/colelab/common/eur_w_ld_chr

#assign array job number to variable for use later
jobnum="${SLURM_ARRAY_TASK_ID}"
jobnum=$(($jobnum - 1))

#Run this if you have munged files and need to get list of mungeed GWAS file names
#only needs to be run once
#ls $in_dir/munged_stats | grep '.sumstats.gz' > $in_dir/munged_file_list.txt

# Run this if you already have a file with file names in it
file_list=$in_dir/whiteXwholegrain_GCST90239658_LDL_Graham_descendingPVE_filenames.txt

# Check file_list exists 
if [ ! -f "$file_list" ]; then
  echo "File list not found." 
  exit 1
fi

# For debugging
echo "Job number: $jobnum"
echo "Input directory: $in_dir"
echo "Inside file_list:" cat $file_list

# Get file name and trait
trait_array=(`cat "$file_list"`) #binary - can't look at this on command line
trait_file="${trait_array[$jobnum]}" #binary - can't look at this on command line
trait_name="${trait_file%.ldsc.imputed_v3.both_sexes.tsv.bgz}" #remove stuff after %
remove_path="/pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/LDSC/"
trait_name="${trait_name/$remove_path/}"
echo "Trait name: $trait_name"

# Get list of files for input to LDSC
len=${#trait_array[@]} # get length of file list
ldsc_files=${trait_array[@]:$jobnum:$len} # Make sub-arrays to do pairwise comparisons
ldsc_file_list=`cat $ldsc_files | tr '\n' ',' | sed 's/,$//'`
#ldsc_file_list=`cat $ldsc_files | sed -z 's/\n/,/g;s/,$/\n/'`
#echo "LDSC File list: $ldsc_file_list"

# Status update
echo "Starting with LDSC... "

#calc ldsc correlations
ldsc.py \
--rg $ldsc_file_list \
--ref-ld-chr $ref_dir \
--w-ld-chr $ref_dir \
--out $trait_name \


## command to pull out final correlation table and add it to overall results file
echo "finished with LDSC... adding correlation table to overall output table file"

if  [ $jobnum == 0 ]
then
	grep 'p1' $out_dir/"$trait_name".log > $out_dir/all_corrs.txt
fi

grep -A90 'p1' $out_dir/"$trait_name".log | grep '/pl' >> $out_dir/all_corrs.txt