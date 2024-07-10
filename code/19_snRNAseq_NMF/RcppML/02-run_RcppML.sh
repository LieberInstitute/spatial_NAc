#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --time=48:0:0
#SBATCH --mem=60G
#SBATCH --job-name=01-select_k_CV
#SBATCH -o ../../../processed-data/19_snRNAseq_NMF/RcppML/logs/01-select_k_CV.log
#SBATCH -e ../../../processed-data/19_snRNAseq_NMF/RcppML/logs/01-select_k_CV.log

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
Rscript 01-select_k_CV.R 

echo "**** Job ends ****"
date
