#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=12:0:0
#SBATCH --mem=60G
#SBATCH --job-name=05-run_GSEA
#SBATCH -o ../../../processed-data/19_snRNAseq_NMF/RcppML/logs/05-run_GSEA.log
#SBATCH -e ../../../processed-data/19_snRNAseq_NMF/RcppML/logs/05-run_GSEA.log

set -e
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
Rscript 05-gene_set_enrichment.R -g $gene_selection

echo "**** Job ends ****"
date
