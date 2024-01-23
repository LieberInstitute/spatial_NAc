#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --job-name=02-explore_results
#SBATCH -o ../../processed-data/10_precast/logs/02-explore_results.log
#SBATCH -e ../../processed-data/10_precast/logs/02-explore_results.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.3
Rscript 02-explore_results.R

echo "**** Job ends ****"
date
