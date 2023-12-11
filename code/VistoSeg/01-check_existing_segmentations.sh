#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH --job-name=01-check_existing_segmentations
#SBATCH -o ../../processed-data/VistoSeg/logs/01-check_existing_segmentations.log
#SBATCH -e ../../processed-data/VistoSeg/logs/01-check_existing_segmentations.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.3
Rscript 01-check_existing_segmentations.R

echo "**** Job ends ****"
date
