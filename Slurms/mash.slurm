#!/bin/bash

#SBATCH -J mash                           # Job name
#SBATCH -o mash_job.%j.out                        # Name of stdout output file (%j expands to jobId)
#SBATCH -p flex                               # Queue name
#SBATCH -N 1                                    # Total number of nodes requested
#SBATCH -n 1                                    # Total number of mpi tasks requested
#SBATCH -t 01:00:00                             # Run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=carriezhu@utexas.edu        # Address email notifications
#SBATCH --mail-type=all                         # Email at begin and end of job

# Change directories to your $SCRATCH directory where your executable is
cd $SCRATCH/Scripts

# load modules and directories
module load Rstats
# export paths
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SCRATCH/gsl/lib
export PATH=$PATH:$SCRATCH/gsl/bin 		# to make gsl-config avail

#launch script executable
./mash.r