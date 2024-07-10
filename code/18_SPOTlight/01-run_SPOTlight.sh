#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --time=72:0:0
#SBATCH --mem=200G
#SBATCH --job-name=01-run_SPOTlight
#SBATCH -o ../../processed-data/18_SPOTlight/logs/01-run_SPOTlight.log
#SBATCH -e ../../processed-data/18_SPOTlight/logs/01-run_SPOTlight.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load r_nac
Rscript 01-run_SPOTlight.R

echo "**** Job ends ****"
date
