#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=48:0:0
#SBATCH --mem=250G
#SBATCH --mail-user=pravich2@jh.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=01-project_sn_nmf
#SBATCH -o ../../../processed-data/16_transfer_learning/02_target_projections/logs/01-project_sn_nmf_human_NAc_all_genes.log
#SBATCH -e ../../../processed-data/16_transfer_learning/02_target_projections/logs/01-project_sn_nmf_human_NAc_all_genes.log

set -e
gene_selection="all_genes"
dataName="human_NAc"
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load r_nac
Rscript 01-project_sn_nmf.R -g $gene_selection -d $dataName

echo "**** Job ends ****"
date