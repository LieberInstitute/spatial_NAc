#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=60G
#SBATCH --job-name=01-run_precast
#SBATCH -o ../../processed-data/10_precast/logs/01-run_precast_%a.log
#SBATCH -e ../../processed-data/10_precast/logs/01-run_precast_%a.log
#SBATCH --array=3-28%3

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.3
Rscript 01-run_precast.R

echo "**** Job ends ****"
date
