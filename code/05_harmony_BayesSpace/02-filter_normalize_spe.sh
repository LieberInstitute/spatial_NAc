#!/bin/bash

#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=80G
#SBATCH --job-name=filter_normalize_spe
#SBATCH -o ../../processed-data/05_harmony_BayesSpace/logs/02-filter_normalize_spe.log
#SBATCH -e ../../processed-data/05_harmony_BayesSpace/logs/02-filter_normalize_spe.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Using $SLURM_CPUS_ON_NODE cores"

module load conda_R/4.3
Rscript 02-filter_normalize_spe.R

echo "**** Job ends ****"
date
