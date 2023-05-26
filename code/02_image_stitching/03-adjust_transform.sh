#!/bin/bash
#$ -cwd
#$ -N "adjust_transform"
#$ -o /dev/null
#$ -e /dev/null
#$ -l mf=5G,h_vmem=5G

slide="V11U23-406"

log_path="../../processed-data/02_image_stitching/03-adjust_transform_${slide}.log"

{
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load loopy/1.0.0-next.24
python 03-adjust_transform.py $slide

echo "**** Job ends ****"
date
} > $log_path 2>&1
