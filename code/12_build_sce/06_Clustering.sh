#!/bin/bash
#
#SBATCH --job-name=06_clustering
#SBATCH --output=logs/06_clustering.log
#SBATCH --error=logs/06_clustering.log
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=2
#
# Mimimum memory required per allocated  CPU
#SBATCH --mem-per-cpu=50G
#
# Send mail to the email address when the job fails
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date

#load R
module load r_nac

#run the Rjob
Rscript 06_Clustering.R 

echo "********* Job Ends *********"
date
