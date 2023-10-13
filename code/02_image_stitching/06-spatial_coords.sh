#!/bin/bash -l

#$ -cwd
#$ -N "spatial_coords"
#$ -o /dev/null
#$ -e /dev/null
#$ -l mf=5G,h_vmem=5G,h_fsize=50G

#SBATCH -p shared
#SBATCH --mem=5G
#SBATCH --job-name=spatial_coords
#SBATCH -o /dev/null
#SBATCH -e /dev/null

donor=Br6432

if [[ ! -z $SLURMD_NODENAME ]]; then
    job_id=$SLURM_JOB_ID
    job_name=$SLURM_JOB_NAME
    node_name=$SLURMD_NODENAME
else
    job_id=$JOB_ID
    job_name=$JOB_NAME
    node_name=$HOSTNAME
fi

log_path="../../processed-data/02_image_stitching/06-spatial_coords_${donor}.log"

{
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Node name: ${node_name}"

module load conda_R/4.3
Rscript 06-spatial_coords.R -d $donor

echo "**** Job ends ****"
date
} > $log_path 2>&1
