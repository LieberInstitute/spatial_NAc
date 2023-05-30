#!/bin/bash
#$ -cwd
#$ -N "export_adjustments"
#$ -o ../../processed-data/02_image_stitching/04-export_adjustments.log
#$ -e ../../processed-data/02_image_stitching/04-export_adjustments.log
#$ -l mf=5G,h_vmem=5G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load loopy/1.0.0-next.24
python 04-export_adjustments.py

echo "**** Job ends ****"
date
