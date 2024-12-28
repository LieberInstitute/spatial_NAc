#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=12:0:0
#SBATCH --mem=60G
#SBATCH --job-name=05-run_GSEA
#SBATCH -o ../../../processed-data/16_transfer_learning/01_process_reference/logs/05-gene_set_enrichment_human_all_genes.log
#SBATCH -e ../../../processed-data/16_transfer_learning/01_process_reference/logs/05-gene_set_enrichment_human_all_genes.log

set -e
gene_selection="all_genes"
data="human_NAc"
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load r_nac
Rscript 05-gene_set_enrichment.R -d $data -g $gene_selection

echo "**** Job ends ****"
date
