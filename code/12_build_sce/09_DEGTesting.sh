#!/bin/bash
#
#SBATCH --job-name=09_CT_DEG
#SBATCH --output=logs/09_CT_DEG.log
#SBATCH --error=logs/09_CT_DEG.log
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=2
#
# Mimimum memory required per allocated  CPU
#SBATCH --mem-per-cpu=25G
#
# Send mail to the email address when the job fails
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date

#load R
module load r_nac
#run the Rjob
Rscript 09_DEGTesting.R 

echo "********* Job Ends *********"
date
