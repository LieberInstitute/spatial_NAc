#!/bin/bash
#$ -cwd
#$ -N "spatial_coords"
#$ -o /dev/null
#$ -e /dev/null
#$ -l mf=5G,h_vmem=5G,h_fsize=50G

donor=

log_path="../../processed-data/02_image_stitching/06-spatial_coords_${donor}.log"

{
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.3
Rscript 06-spatial_coords.R -d $donor

echo "**** Job ends ****"
date
} > $log_path 2>&1
