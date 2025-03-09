#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=48:0:0
#SBATCH --mem=60G
#SBATCH --job-name=02_sn_NHP_comparison
#SBATCH -o ../../processed-data/17_cross_species_comparison/logs/02-sn_NHP_comparison.log
#SBATCH -e ../../processed-data/17_cross_species_comparison/logs/02-sn_NHP_comparison.log

set -e
#   'TRUE' or 'FALSE': if TRUE, only use marker genes for the known cell types
subset_neurons=FALSE

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load r_nac
Rscript 02-sn_NHP_comparison.R -n $subset_neurons

echo "**** Job ends ****"
date

