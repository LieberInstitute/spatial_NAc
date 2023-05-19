#!/bin/bash
#$ -cwd
#$ -N "samui_test"
#$ -o /dev/null
#$ -e /dev/null
#$ -l mf=40G,h_vmem=40G,h_fsize=50G

slide="V11U23-404"
mode="initial"

log_path="../../processed-data/02_image_stitching/01-samui_test_${slide}_${mode}.log"

{
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load loopy/1.0.0-next.24
python 01-samui_test.py $slide $mode

echo "**** Job ends ****"
date
} > $log_path 2>&1
