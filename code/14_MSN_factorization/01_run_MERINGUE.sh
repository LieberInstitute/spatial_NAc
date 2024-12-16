#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=60G
#SBATCH --time=48:0:0
#SBATCH --job-name=01-run_MERINGUE
#SBATCH -o ../../processed-data/14_MSN_factorization/01_meringue/logs/01-run_meringue_%a.log
#SBATCH -e ../../processed-data/14_MSN_factorization/01_meringue/logs/01-run_meringue_%a.log
#SBATCH --array=10-10%10

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load r_nac
Rscript 01_run_MERINGUE.R

echo "**** Job ends ****"
date
