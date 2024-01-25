#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --job-name=06-gather_nnSVG
#SBATCH -o ../../processed-data/05_harmony_BayesSpace/logs/06-gather_nnSVG_precast.log
#SBATCH -e ../../processed-data/05_harmony_BayesSpace/logs/06-gather_nnSVG_precast.log

#   'TRUE' or 'FALSE': if TRUE, explore SVG results from nnSVG using k=2
#   PRECAST clusters as a covariate
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
Rscript 06-gather_nnSVG.R -p $USE_PRECAST

echo "**** Job ends ****"
date
