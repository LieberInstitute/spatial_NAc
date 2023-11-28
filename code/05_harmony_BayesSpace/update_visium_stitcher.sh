#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=update_visium_stitcher
#SBATCH -o ../../processed-data/05_harmony_BayesSpace/logs/update_visium_stitcher.log
#SBATCH -e ../../processed-data/05_harmony_BayesSpace/logs/update_visium_stitcher.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.3
Rscript update_visium_stitcher.R

echo "**** Job ends ****"
date
