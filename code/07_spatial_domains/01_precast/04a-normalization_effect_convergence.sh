#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --time=48:0:0
#SBATCH --mem=250G
#SBATCH --job-name=04a_convergence_normalization
#SBATCH -o ../../../processed-data/07_spatial_domains/01_precast/logs/04a_normalization_convergence.log
#SBATCH -e ../../../processed-data/07_spatial_domains/01_precast/logs/04a_normalization_convergence.log

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
ulimit -s unlimited
Rscript 04a-normalization_effect_convergence.R

echo "**** Job ends ****"
date
