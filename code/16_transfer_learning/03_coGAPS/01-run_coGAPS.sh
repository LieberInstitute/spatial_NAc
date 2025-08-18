#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=48:0:0
#SBATCH --mem=250G
#SBATCH --mail-user=pravich2@jh.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=01-run_coGAPS
#SBATCH -o ../../../processed-data/16_transfer_learning/03_coGAPS/logs/01-run_coGAPS_rat_case_control_repeated.log
#SBATCH -e ../../../processed-data/16_transfer_learning/03_coGAPS/logs/01-run_coGAPS_rat_case_control_repeated.log

set -e
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load r_nac
Rscript 01-run_coGAPS.R

echo "**** Job ends ****"
date