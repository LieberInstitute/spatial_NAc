#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=02_emptyDrops
#SBATCH --output=../../processed-data/21_Reiner_snRNAseq/logs/emptydrops.%a.log
#SBATCH --error=../../processed-data/21_Reiner_snRNAseq/logs/emptydrops.%a.log
#SBATCH --mem=40G
#SBATCH --array=1-22



echo "********* Job Starts *********"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load r_nac
Rscript 02_calculate_droplet_scores.R

echo "********* Job Ends *********"
date
