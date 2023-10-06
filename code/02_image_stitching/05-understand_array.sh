#!/bin/bash -l

#$ -cwd
#$ -N "understand_array"
#$ -o ../../processed-data/02_image_stitching/05-understand_array.log
#$ -e ../../processed-data/02_image_stitching/05-understand_array.log
#$ -l mf=5G,h_vmem=5G

#SBATCH -p shared
#SBATCH --mem=5G
#SBATCH --job-name=understand_array
#SBATCH -o ../../processed-data/02_image_stitching/05-understand_array.log
#SBATCH -e ../../processed-data/02_image_stitching/05-understand_array.log

if [[ ! -z $SLURMD_NODENAME ]]; then
    job_id=$SLURM_JOB_ID
    job_name=$SLURM_JOB_NAME
    node_name=$SLURMD_NODENAME
else
    job_id=$JOB_ID
    job_name=$JOB_NAME
    node_name=$HOSTNAME
fi

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Node name: ${node_name}"

module load conda_R/4.3
Rscript 05-understand_array.R

echo "**** Job ends ****"
date
