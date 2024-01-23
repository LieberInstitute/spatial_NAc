#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --job-name=05-run_nnSVG
#SBATCH -o ../../processed-data/05_harmony_BayesSpace/logs/05-run_nnSVG_precast_%a.log
#SBATCH -e ../../processed-data/05_harmony_BayesSpace/logs/05-run_nnSVG_precast_%a.log
#SBATCH --array=1-10%5

#   'TRUE' or 'FALSE': if TRUE, find SVGs within k=2 PRECAST clusters. If FALSE,
#   run nnSVG without covariates
USE_PRECAST=TRUE

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.3
Rscript 05-run_nnSVG.R -p $USE_PRECAST

echo "**** Job ends ****"
date
