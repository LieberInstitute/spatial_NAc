#!/bin/bash

#SBATCH -q shared
#SBATCH --mem=70G
#SBATCH --job-name=r_to_python
#SBATCH -o ../../processed-data/08_spot_deconvo/logs/01-r_to_python.log
#SBATCH -e ../../processed-data/08_spot_deconvo/logs/01-r_to_python.log

if [[ ! -z $SLURMD_NODENAME ]]; then
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

module load r_nac
Rscript 01-r_to_python.R

echo "**** Job ends ****"
date
