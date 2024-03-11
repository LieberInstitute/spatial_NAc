#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -t 1:00:00
#SBATCH --job-name=05-precast_figures
#SBATCH -o ../../processed-data/10_precast/logs/05-precast_figures.log
#SBATCH -e ../../processed-data/10_precast/logs/05-precast_figures.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.3.x
Rscript 05-precast_figures.R

echo "**** Job ends ****"
date
