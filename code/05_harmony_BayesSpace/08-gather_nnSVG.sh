#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --job-name=08-gather_nnSVG
#SBATCH -o ../../processed-data/05_harmony_BayesSpace/logs/08-gather_nnSVG.log
#SBATCH -e ../../processed-data/05_harmony_BayesSpace/logs/08-gather_nnSVG.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load r_nac
Rscript 08-gather_nnSVG.R

echo "**** Job ends ****"
date
