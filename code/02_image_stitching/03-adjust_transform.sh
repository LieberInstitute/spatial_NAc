#!/bin/bash
#$ -cwd
#$ -N "adjust_transform"
#$ -o /dev/null
#$ -e /dev/null
#$ -l mf=5G,h_vmem=5G

#SBATCH -p shared
#SBATCH --mem=5G
#SBATCH --job-name=adjust_transform
#SBATCH -o /dev/null
#SBATCH -e /dev/null

donor=Br3942

if [[ ! -z $SLURMD_NODENAME ]]; then
    job_id=$SLURM_JOB_ID
    job_name=$SLURM_JOB_NAME
    node_name=$SLURMD_NODENAME
    module_name=samui
else
    job_id=$JOB_ID
    job_name=$JOB_NAME
    node_name=$HOSTNAME
    module_name=loopy
fi

log_path="../../processed-data/02_image_stitching/03-adjust_transform_${donor}.log"

{
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Node name: ${node_name}"

module load ${module_name}/1.0.0-next.24
python 03-adjust_transform.py $donor

echo "**** Job ends ****"
date
} > $log_path 2>&1
