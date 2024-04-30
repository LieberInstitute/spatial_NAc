#!/bin/bash
#
#SBATCH --job-name=10_sub_DimRed
#SBATCH --output=logs/10_sub_DimRed.log
#SBATCH --error=logs/10_sub_DimRed.log
#
#
# Mimimum memory required per allocated  CPU
#SBATCH --mem=20G
#
# Send mail to the email address when the job fails
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date

#load R
module load r_nac
#run the Rjob
Rscript 10_neuronsubcluster_DimRed.R

echo "********* Job Ends *********"
date
