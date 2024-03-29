#!/bin/bash
#
#SBATCH --job-name=07_DEG
#SBATCH --output=logs/07_PrelimDEG.log
#SBATCH --error=logs/07_PrelimDEG.log
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=2
#
# Mimimum memory required per allocated  CPU
#SBATCH --mem-per-cpu=5G
#
# Send mail to the email address when the job fails
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date

#load R
module load r_nac
#run the Rjob
Rscript 07_PrelimDEGTesting.R

echo "********* Job Ends *********"
date
