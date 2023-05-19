#!/bin/bash
#$ -cwd
#$ -N "samui_test"
#$ -o ../../processed-data/02_image_stitching/01-samui_test_V11U08-082_adjusted.log
#$ -e ../../processed-data/02_image_stitching/01-samui_test_V11U08-082_adjusted.log
#$ -l mf=40G,h_vmem=40G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load loopy/1.0.0-next.24
python 01-samui_test.py

echo "**** Job ends ****"
date
