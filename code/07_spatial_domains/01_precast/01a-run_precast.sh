#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --time=56:0:0
#SBATCH --mem=80G
#SBATCH --job-name=01-run_precast_BIC
#SBATCH -o ../../processed-data/07_spatial_domains/01_precast/logs/01-run_precast_BIC.log
#SBATCH -e ../../processed-data/07_spatial_domains/01_precast/logs/01-run_precast_BIC.log


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
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load r_nac
ulimit -s unlimited
Rscript 01a-run_precast.R -n $nnSVG_TYPE

echo "**** Job ends ****"
date
