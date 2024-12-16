#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=48:0:0
#SBATCH --mem=60G
#SBATCH --job-name=02-run_RcppML
#SBATCH -o ../../../processed-data/16_transfer_learning/01_process_reference/logs/02-run_RcppML_all_genes_rat.log
#SBATCH -e ../../../processed-data/16_transfer_learning/01_process_reference/logs/02-run_RcppML_all_genes_rat.log

set -e
gene_selection="all_genes"
dataName="rat_case_control"
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load r_nac
Rscript 02-NMF_reference.R -g $gene_selection -d $dataName

echo "**** Job ends ****"
date
