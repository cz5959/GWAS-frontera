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
