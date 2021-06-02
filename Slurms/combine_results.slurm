#!/bin/bash

#SBATCH -J GWAS_combine_results            		# Job name
#SBATCH -o combine_results_job.%j.out           		# Name of stdout output file (%j expands to jobId)
#SBATCH -p flex   	      			# Queue name
#SBATCH -N 1                  			# Total number of nodes requested
#SBATCH -n 1                 			# Total number of mpi tasks requested
#SBATCH -t 0:30:00           			# Run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=carriezhu@utexas.edu	# Address email notifications
#SBATCH --mail-type=all				# Email at begin and end of job


# Change directories to your $SCRATCH directory where your executable is
cd $SCRATCH/Scripts

#load modules and directories

#launch script executable
#./Association_script.sh
./combine_results.sh
