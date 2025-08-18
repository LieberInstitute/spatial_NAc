#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=48:0:0
#SBATCH --mem=60G
#SBATCH --job-name=01-run_NMF
#SBATCH -o ../../../processed-data/16_transfer_learning/04_cellType_NMF/logs/01-run_NMF_rat_case_control_cocaine_acute_30.log
#SBATCH -e ../../../processed-data/16_transfer_learning/04_cellType_NMF/logs/01-run_NMF_rat_case_control_cocaine_acute_30.log

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
Rscript 01-run_NMF.R -d rat_case_control_cocaine_acute -n 30 -s FALSE

echo "**** Job ends ****"
date
