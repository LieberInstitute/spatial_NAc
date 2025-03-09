#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=60G
#SBATCH --job-name=02-prepare_myRCTD
#SBATCH -o ../../../processed-data/08_spot_deconvo/01_RCTD/logs/02-prepare_myRCTD_all_genes_%a.log
#SBATCH -e ../../../processed-data/08_spot_deconvo/01_RCTD/logs/02-prepare_myRCTD_all_genes_%a.log
#SBATCH --array=1-10%10

#   'TRUE' or 'FALSE': if TRUE, only use marker genes for the known cell types
marker_genes=FALSE

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load r_nac
Rscript 02-prepare_myRCTD.R -m $marker_genes

echo "**** Job ends ****"
date
