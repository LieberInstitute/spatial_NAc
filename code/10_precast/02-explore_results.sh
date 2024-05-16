#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --job-name=02-explore_results
#SBATCH -o ../../processed-data/10_precast/logs/02-explore_results_precast.log
#SBATCH -e ../../processed-data/10_precast/logs/02-explore_results_precast.log


#   'TRUE' or 'FALSE': if TRUE, use nnSVGs which were obtained after controlling for precast k = 2 clusters, else use default
nnSVG_TYPE=TRUE

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load r_nac
Rscript 02-explore_results.R -n $nnSVG_TYPE

echo "**** Job ends ****"
date
