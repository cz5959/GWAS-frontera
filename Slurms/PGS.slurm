#!/bin/bash

#SBATCH -J PGS            		        # Job name
#SBATCH -o PGS.%j.out           		# Name of stdout output file (%j expands to jobId)
#SBATCH -p flex   	      			    # Queue name
#SBATCH -N 1                  			# Total number of nodes requested
#SBATCH -n 1                 			# Total number of mpi tasks requested
#SBATCH -t 1:00:00           			# Run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=carriezhu@utexas.edu	# Address email notifications
#SBATCH --mail-type=all				# Email at begin and end of job

# Change directories to your $SCRATCH directory where your executable is
cd $SCRATCH/Scripts/PGS

#load modules and directories
module load Rstats

#launch script executable

declare -a arr=("height" "bmi" "testosterone" "waist_to_hip_bmi" "RBC_count" "IGF1" "creatinine" "weight" "calcium" \
"protein_total" "urea" "SHBG" "whole_body_fat_mass" "FVC_best" "HbA1c" "urate" "arm_fatfree_mass_L" "arm_fatfree_mass_R" \
"eosinophil_perc" "lymphocyte_perc" "waist_circ" "hip_circ" "waist_to_hip" "diastolicBP_auto" "systolicBP_auto" "albumin" "pulse_rate")

for i in "${arr[@]}"
do
        mkdir -p $SCRATCH/GWAS_Results/${i}/PGS
        ./PGS_testset_1.R -p $i
done
