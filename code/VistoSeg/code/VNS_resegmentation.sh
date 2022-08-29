#!/bin/bash
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G,h_stack=256M,h_fsize=100G
#$ -o /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/code/VistoSeg/code/logs/$TASK_ID_VNS_reseg.txt
#$ -e /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/code/VistoSeg/code/logs/$TASK_ID_VNS_reseg.txt
#$ -m e
#$ -M heenadivecha@gmail.com
#$ -t 1-16
#$ -tc 1

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: $(cat /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/code/VistoSeg/code/VNS_list.txt | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")"
echo "****"

module load matlab/R2019a


toolbox='/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/code/VistoSeg/code'
fname=$(cat /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/code/VistoSeg/code/VNS_list.txt | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")

matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), VNS_resegmentation('$fname',5)"

echo "**** Job ends ****"
date



