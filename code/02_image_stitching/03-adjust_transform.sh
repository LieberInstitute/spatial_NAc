#!/bin/bash
#$ -cwd
#$ -N "adjust_transform"
#$ -o /dev/null
#$ -e /dev/null
#$ -l mf=5G,h_vmem=5G

slide="V12D07-078"
arrays="A1_B1_C1_D1"

log_path="../../processed-data/02_image_stitching/03-adjust_transform_${slide}_${arrays}.log"

{
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load loopy/1.0.0-next.24
python 03-adjust_transform.py $slide $arrays

echo "**** Job ends ****"
date
} > $log_path 2>&1
