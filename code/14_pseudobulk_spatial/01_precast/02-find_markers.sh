#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH -t 2-00:00
#SBATCH --job-name=02-find_markers
#SBATCH -o ../../../processed-data/14_pseudobulk_spatial/01_precast/logs/02-find_markers_donor_k%a.log
#SBATCH -e ../../../processed-data/14_pseudobulk_spatial/01_precast/logs/02-find_markers_donor_k%a.log
#SBATCH --array=3-18%15

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load r_nac
agg_level="sample_id"
Rscript 02-find_markers.R -a $agg_level

echo "**** Job ends ****"
date
