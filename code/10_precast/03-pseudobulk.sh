#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -t 1-00:00
#SBATCH --job-name=03-pseudobulk
#SBATCH -o ../../processed-data/10_precast/logs/03-pseudobulk_donor_k%a.log
#SBATCH -e ../../processed-data/10_precast/logs/03-pseudobulk_donor_k%a.log
#SBATCH --array=3-28%15

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

agg_level="sample_id"
module load r_nac
Rscript 03-pseudobulk.R -a $agg_level

echo "**** Job ends ****"
date
