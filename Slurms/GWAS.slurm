#!/bin/bash

#SBATCH -J GWAS           # Job name
#SBATCH -o GWAS.o%j       # Name of stdout output file
#SBATCH -e GWAS.e%j       # Name of stderr error file
#SBATCH -p small            # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 07:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH --mail-user=carriezhu@utexas.edu

# Any other commands must follow all #SBATCH directives...
cd $SCRATCH/Scripts

# phenotypes
declare -a arr=("creatinine" "RBC_count")

# load modules
module load python3
module load Rstats

# Launch script
for i in "${arr[@]}"
do
        ./Association_script.sh $i
        ./combine_results.sh $i
        ./Manhattan_plot.py $i
        ./QQ_plot.R $i
done

script="Association_script.sh"
./$script -p albumin & ./$script -p arm_fatfree_mass_L & ./$script -p arm_fatfree_mass_R & ./$script -p bmi & ./$script -p calcium &
./$script -p creatinine & ./$script -p diastolicBP_auto & ./$script -p eosinophil_perc & ./$script -p FVC_best & ./$script -p HbA1c &
./$script -p height & ./$script -p hip_circ & ./$script -p IGF1 & ./$script -p lymphocyte_perc & ./$script -p protein_total &
./$script -p pulse_rate & ./$script -p RBC_count & ./$script -p SHBG & ./$script -p systolicBP_auto & ./$script -p testosterone &
./$script -p urate & ./$script -p urea & ./$script -p waist_circ & ./$script -p waist_to_hip & ./$script -p weight &
./$script -p whole_body_fat_mass & ./$script -p wth_bmi_adj &
wait 
echo done 


PATTERN="*_pm.txt"
if ls $PATTERN 1> /dev/null 2>&1; then
    echo "Files exist"
fi