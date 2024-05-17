#!/bin/bash
#
#SBATCH --job-name=08_anno
#SBATCH --output=logs/08_anno.log
#SBATCH --error=logs/08_anno.log
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=2
#
# Mimimum memory required per allocated  CPU
#SBATCH --mem-per-cpu=35G
#
# Send mail to the email address when the job fails
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date

#load R
module load r_nac
#run the Rjob
Rscript 08_annotation.R

echo "********* Job Ends *********"
date
