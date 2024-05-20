#!/bin/bash
#
#SBATCH --job-name=09.5_CT_gmr
#SBATCH --output=logs/09.5_CT_gmr.log
#SBATCH --error=logs/09.5_CT_gmr.log
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=2
#
# Mimimum memory required per allocated  CPU
#SBATCH --mem-per-cpu=30G
#
# Send mail to the email address when the job fails
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date

#load R
module load r_nac
#run the Rjob
Rscript 09.5_meanratio_CellType.R

echo "********* Job Ends *********"
date
