#!/bin/bash

#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=7-00:00:00
#SBATCH --mem=100G
#SBATCH --job-name=01-select_k_CV
#SBATCH -o ../../../processed-data/16_transfer_learning/01_process_reference/logs/01_cross_validation_rat_case_control_morphine_acute.log
#SBATCH -e ../../../processed-data/16_transfer_learning/01_process_reference/logs/01_cross_validation_rat_case_control_morphine_acute.log

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
data_set="rat_case_control_morphine_acute"
Rscript 01-select_k_CV.R -d $data_set

echo "**** Job ends ****"
date
