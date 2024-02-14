#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -t 1-00:00
#SBATCH --job-name=03-pseudobulk_k2
#SBATCH -o ../../processed-data/10_precast/logs/03-pseudobulk_k2.log
#SBATCH -e ../../processed-data/10_precast/logs/03-pseudobulk_k2.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.3
Rscript 03-pseudobulk.R -k 2

echo "**** Job ends ****"
date
