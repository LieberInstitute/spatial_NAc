#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=48:0:0
#SBATCH --mem=60G
#SBATCH --job-name=03_spatial_rat_comparison
#SBATCH -o ../../processed-data/17_cross_species_comparison/logs/03-spatial_rat_comparison.log
#SBATCH -e ../../processed-data/17_cross_species_comparison/logs/03-spatial_rat_comparison.log

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
Rscript 03-spatial_rat_comparison.R

echo "**** Job ends ****"
date

