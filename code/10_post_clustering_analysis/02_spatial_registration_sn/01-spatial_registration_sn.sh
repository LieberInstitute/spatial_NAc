#!/bin/bash

#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=80G
#SBATCH -t 3-00:00
#SBATCH --job-name=01-spatial_registration_sn
#SBATCH -o ../../../processed-data/16_spatial_registration_sn/logs/01-spatial_registration_sn.log
#SBATCH -e ../../../processed-data/16_spatial_registration_sn/logs/01-spatial_registration_sn.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load r_nac
Rscript 01-spatial_registration_sn.R 

echo "**** Job ends ****"
date
