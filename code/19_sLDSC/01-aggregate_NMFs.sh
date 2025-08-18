#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=10:0:0
#SBATCH --mem=80G
#SBATCH --mail-user=pravich2@jh.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=01-aggregate_NMF
#SBATCH -o ../../processed-data/19_sLDSC/human_NAc_NMF/logs/01-aggregate_NMFs.log
#SBATCH -e ../../processed-data/19_sLDSC/human_NAc_NMF/logs/01-aggregate_NMFs.log

set -e
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load r_nac
Rscript 01-aggregate_NMFs.R

echo "**** Job ends ****"
date