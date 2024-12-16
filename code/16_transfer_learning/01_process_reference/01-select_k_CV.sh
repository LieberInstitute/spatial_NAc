#!/bin/bash

#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=10-00:00:00
#SBATCH --mem=80G
#SBATCH --job-name=01-select_k_CV
#SBATCH -o ../../../processed-data/16_transfer_learning/01_process_reference/logs/01_cross_validation_rat_all_genes.log
#SBATCH -e ../../../processed-data/16_transfer_learning/01_process_reference/logs/01_cross_validation_rat_all_genes.log

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
gene_selection="all_genes"
data_set="rat_case_control"
Rscript 01-select_k_CV.R -g $gene_selection -d $data_set

echo "**** Job ends ****"
date
