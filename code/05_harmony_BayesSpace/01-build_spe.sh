#!/bin/bash

#$ -cwd
#$ -N "build_spe"
#$ -o ../../processed-data/05_harmony_BayesSpace/logs/01-build_spe.log
#$ -e ../../processed-data/05_harmony_BayesSpace/logs/01-build_spe.log
#$ -pe local 4
#$ -l mf=10G,h_vmem=10G,h_fsize=50G

#SBATCH -q shared
#SBATCH -c 4
#SBATCH --mem=40G
#SBATCH --job-name=build_spe
#SBATCH -o ../../processed-data/05_harmony_BayesSpace/logs/01-build_spe.log
#SBATCH -e ../../processed-data/05_harmony_BayesSpace/logs/01-build_spe.log

USE_SLURM=1

if [[ $USE_SLURM -eq 1 ]]; then
    job_id=$SLURM_JOB_ID
    job_name=$SLURM_JOB_NAME
else
    job_id=$JOB_ID
    job_name=$JOB_NAME
fi

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.3
Rscript 01-build_spe.R

echo "**** Job ends ****"
date
