#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH --job-name=02-explore_results
#SBATCH -o ../../../processed-data/07_spatial_domains/01_precast/logs/02-explore_results_precast_%a.log
#SBATCH -e ../../../processed-data/07_spatial_domains/01_precast/logs/02-explore_results_precast_%a.log
#SBATCH --array=1-5%5

#   'TRUE' or 'FALSE': if TRUE, use nnSVGs which were obtained after controlling for precast k = 2 clusters, else use default
nnSVG_TYPE=TRUE
use_random_start=TRUE
set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load r_nac
Rscript 02-explore_results.R -n $nnSVG_TYPE -r $use_random_start

echo "**** Job ends ****"
date
