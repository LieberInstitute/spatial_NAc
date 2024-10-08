#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --time=48:0:0
#SBATCH --mem=60G
#SBATCH --job-name=01-run_precast_random
#SBATCH -o ../../processed-data/10_precast/logs/01-run_precast_random3_%a.log
#SBATCH -e ../../processed-data/10_precast/logs/01-run_precast_random3_%a.log
#SBATCH --array=3-15%15

#   'TRUE' or 'FALSE': if TRUE, use nnSVGs which were obtained after controlling for precast k = 2 clusters, else use default
nnSVG_TYPE=TRUE
SPECIFY_K=TRUE
RANDOM_SEED=3

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
Rscript 01-run_precast.R -n $nnSVG_TYPE -k $SPECIFY_K -r $RANDOM_SEED

echo "**** Job ends ****"
date
