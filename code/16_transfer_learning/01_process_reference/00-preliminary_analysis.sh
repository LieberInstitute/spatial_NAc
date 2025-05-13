#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=48:0:0
#SBATCH --mem=60G
#SBATCH --job-name=00-preliminary_analysis
#SBATCH -o ../../../processed-data/16_transfer_learning/01_process_reference/logs/00-prelim_analysis_all_genes_rat_morphine_repeated.log
#SBATCH -e ../../../processed-data/16_transfer_learning/01_process_reference/logs/00-prelim_analysis_all_genes_rat_morphine_repeated.log

set -e
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
data_set="rat_case_control_morphine_repeated"
module load r_nac
Rscript 00-preliminary_analysis.R -d $data_set

echo "**** Job ends ****"
date

