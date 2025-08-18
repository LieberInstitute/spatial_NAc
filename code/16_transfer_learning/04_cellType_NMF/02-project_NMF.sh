#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=8:0:0
#SBATCH --mem=80G
#SBATCH --job-name=02-project_NMF
#SBATCH -o ../../../processed-data/16_transfer_learning/04_cellType_NMF/logs/02-project_NMF_rat_case_control_morphine_repeated_50.log
#SBATCH -e ../../../processed-data/16_transfer_learning/04_cellType_NMF/logs/02-project_NMF_rat_case_control_morphine_repeated_50.log

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
Rscript 02-project_NMF.R -d rat_case_control_morphine_repeated -n 50 -s FALSE

echo "**** Job ends ****"
date
