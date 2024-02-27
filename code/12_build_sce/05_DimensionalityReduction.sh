#!/bin/bash
#
#SBATCH --job-name=05_DimRed
#SBATCH --output=logs/05_DimensionalityReduction.log
#SBATCH --error=logs/05_DimensionalityReduction.log
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=2
#
# Mimimum memory required per allocated  CPU
#SBATCH --mem-per-cpu=80G
#
# Send mail to the email address when the job fails
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date

#load R
module load r_nac

#run the Rjob
Rscript 05_DimensionalityReduction.R

echo "********* Job Ends *********"
date
