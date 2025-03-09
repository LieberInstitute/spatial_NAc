#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=60G
#SBATCH --time=48:0:0
#SBATCH --job-name=03-runRCTD
#SBATCH -o ../../../processed-data/08_spot_deconvo/01_RCTD/logs/03-runRCTD_markers_%a.log
#SBATCH -e ../../../processed-data/08_spot_deconvo/01_RCTD/logs/03-runRCTD_markers_%a.log
#SBATCH --array=1-10%10

#   'TRUE' or 'FALSE': if TRUE, use only the marker genes for cell types to run RCTD
marker_genes=TRUE

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load r_nac
Rscript 03-runRCTD.R -m $marker_genes

echo "**** Job ends ****"
date
