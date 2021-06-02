#!/bin/bash

#SBATCH -J linear_regression                           # Job name
#SBATCH -o OLS_job.%j.out                        # Name of stdout output file (%j expands to jobId)
#SBATCH -p normal                               # Queue name
#SBATCH -N 1                                    # Total number of nodes requested
#SBATCH -n 1                                    # Total number of mpi tasks requested
#SBATCH -t 08:00:00                             # Run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=carriezhu@utexas.edu        # Address email notifications
#SBATCH --mail-type=all                         # Email at begin and end of job

# Change directories to your $SCRATCH directory where your executable is
cd /scratch/08005/cz5959/gwas_practice/QC_2

#load modules and directories
module load python3
pip3 install --user sklearn
pip3 install --user statsmodels
pip3 install --user scipy 

#launch script executable
./linear_regression.py
