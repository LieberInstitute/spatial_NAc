#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --job-name=07-run_nnSVG
#SBATCH -o ../../processed-data/05_harmony_BayesSpace/logs/07-run_nnSVG_precast_%a.log
#SBATCH -e ../../processed-data/05_harmony_BayesSpace/logs/07-run_nnSVG_precast_%a.log
#SBATCH --array=1-10%10

#   'TRUE' or 'FALSE': if TRUE, find SVGs within k=2 PRECAST clusters. If FALSE,
#   run nnSVG without covariates
USE_PRECAST=FALSE

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load r_nac
Rscript 07-run_nnSVG.R -p $USE_PRECAST

echo "**** Job ends ****"
date
