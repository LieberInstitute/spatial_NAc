#!/bin/bash
#$ -cwd
#$ -N "spatial_coords"
#$ -o /dev/null
#$ -e /dev/null
#$ -l mf=5G,h_vmem=5G,h_fsize=50G

slide="V12D07-074"
arrays="A1_B1_C1_D1"

log_path="../../processed-data/02_image_stitching/06-spatial_coords_${slide}_${arrays}.log"

{
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.3
Rscript 06-spatial_coords.R -s $slide -a $arrays

echo "**** Job ends ****"
date
} > $log_path 2>&1
