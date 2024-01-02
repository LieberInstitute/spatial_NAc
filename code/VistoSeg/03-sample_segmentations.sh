#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --job-name=03-sample_segmentations
#SBATCH -o ../../processed-data/VistoSeg/logs/03-sample_segmentations.log
#SBATCH -e ../../processed-data/VistoSeg/logs/03-sample_segmentations.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load samui/1.0.0-next.45
python 03-sample_segmentations.py

echo "**** Job ends ****"
date
