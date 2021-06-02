#!/bin/bash
#----------------------------------------------------
# Example SLURM job script to run applications on 
# TACC's Lonestar 5 system.
#
# Example of serial job submission to normal queue
# Last revised: 14 Sep 2020
#----------------------------------------------------

#SBATCH -J GWAS_scatter_plot            		# Job name
#SBATCH -o scatter_plot_job.%j.out           		# Name of stdout output file (%j expands to jobId)
#SBATCH -p flex   	      			# Queue name
#SBATCH -N 1                  			# Total number of nodes requested
#SBATCH -n 1                 			# Total number of mpi tasks requested
#SBATCH -t 0:10:00           			# Run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=carriezhu@utexas.edu	# Address email notifications
#SBATCH --mail-type=all				# Email at begin and end of job


# Change directories to your $SCRATCH directory where your executable is
cd $SCRATCH/Scripts

#load modules and directories
module load python3
#module load Rstats
#module load RstatsPackages

#launch script executable
#./Manhattan_plot.py
./scatter_p.py
#./QQ_plot.R
