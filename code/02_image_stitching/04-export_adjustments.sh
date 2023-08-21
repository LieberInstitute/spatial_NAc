#!/bin/bash -l

#$ -cwd
#$ -N "export_adjustments"
#$ -o ../../processed-data/02_image_stitching/04-export_adjustments.log
#$ -e ../../processed-data/02_image_stitching/04-export_adjustments.log
#$ -l mf=5G,h_vmem=5G

#SBATCH -q shared
#SBATCH --mem=5G
#SBATCH --job-name=export_adjustments
#SBATCH -o ../../processed-data/02_image_stitching/04-export_adjustments.log
#SBATCH -e ../../processed-data/02_image_stitching/04-export_adjustments.log

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

module load loopy/1.0.0-next.24
python 04-export_adjustments.py

echo "**** Job ends ****"
date
