#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=06-check_counts
#SBATCH -o ../../processed-data/VistoSeg/logs/06-check_counts.log
#SBATCH -e ../../processed-data/VistoSeg/logs/06-check_counts.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.3
Rscript 06-check_counts.R

echo "**** Job ends ****"
date
