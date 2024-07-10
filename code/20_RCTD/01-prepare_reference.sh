#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --time=48:0:0
#SBATCH --mem=60G
#SBATCH --job-name=01-select_k_CV
#SBATCH -o ../../processed-data/20_RCTD/logs/01-prepare_reference.log
#SBATCH -e ../../processed-data/20_RCTD/logs/01-prepare_reference.log

set -e
marker_genes=TRUE

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load r_nac
Rscript 01-prepare_reference.R -m marker_genes

echo "**** Job ends ****"
date
