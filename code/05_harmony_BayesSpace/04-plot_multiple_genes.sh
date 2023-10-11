#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=15G
#SBATCH --job-name=04-plot_multiple_genes
#SBATCH -o ../../processed-data/05_harmony_BayesSpace/logs/04-plot_multiple_genes.log
#SBATCH -e ../../processed-data/05_harmony_BayesSpace/logs/04-plot_multiple_genes.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.3
Rscript 04-plot_multiple_genes.R

echo "**** Job ends ****"
date
