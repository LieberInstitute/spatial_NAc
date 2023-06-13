#!/bin/bash
#$ -cwd
#$ -N "understand_array"
#$ -o ../../processed-data/02_image_stitching/05-understand_array.log
#$ -e ../../processed-data/02_image_stitching/05-understand_array.log
#$ -l mf=5G,h_vmem=5G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.3
Rscript 05-understand_array.R

echo "**** Job ends ****"
date
