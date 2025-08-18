#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=8:0:0
#SBATCH --mem=250G
#SBATCH --job-name=01_CCC_w_LIANA
#SBATCH -o /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/22_gene_risk_LR_analysis/01_liana/logs/01-run_LIANA.log
#SBATCH -e /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/22_gene_risk_LR_analysis/01_liana/logs/01-run_LIANA.log

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
Rscript 01_CCC_with_LIANA.R

echo "**** Job ends ****"
date

