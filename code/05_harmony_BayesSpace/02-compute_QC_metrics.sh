#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=250G
#SBATCH --time=24:0:0
#SBATCH --job-name=compute_QC_metrics
#SBATCH --mail-user=pravich2@jh.edu
#SBATCH --mail-type=ALL
#SBATCH -o ../../processed-data/05_harmony_BayesSpace/logs/02-compute_QC_metrics.log
#SBATCH -e ../../processed-data/05_harmony_BayesSpace/logs/02-compute_QC_metrics.log
set -e
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Using $SLURM_CPUS_ON_NODE cores"
module load r_nac
Rscript 02-compute_QC_metrics.R
echo "**** Job ends ****"
date
