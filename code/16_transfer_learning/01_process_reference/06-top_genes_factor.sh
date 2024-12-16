#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=48:0:0
#SBATCH --mem=250G
#SBATCH --job-name=06-run_top_genes
#SBATCH -o ../../../processed-data/16_transfer_learning/01_process_reference/logs/06-top_genes_human_all_genes.log
#SBATCH -e ../../../processed-data/16_transfer_learning/01_process_reference/logs/06-top_genes_human_all_genes.log

set -e
data="human_NAc"
gene_selection="all_genes"
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load r_nac
Rscript 06-top_genes_factor.R -d $data -g $gene_selection

echo "**** Job ends ****"
date
