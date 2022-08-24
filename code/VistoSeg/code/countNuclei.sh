#!/bin/bash
#$ -cwd
#$ -l mem_free=8G,h_vmem=8G,h_fsize=100G
#$ -pe local 8
#$ -N spatial_NAc
#$ -o /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/code/VistoSeg/code/logs/countNuclei.$TASK_ID.txt
#$ -e /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/code/VistoSeg/code/logs/countNuclei.$TASK_ID.txt
#$ -m e
#$ -t 2-16
#$ -tc 10

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: $(cat /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/code/VistoSeg/code/countNuclei_list.txt | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")"
echo "****"


module load matlab/R2019a


toolbox='/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/code/VistoSeg/code'

## Read parameters
mask=$(awk 'BEGIN {FS="\t"} {print $1}' countNuclei_list.txt | awk "NR==${SGE_TASK_ID}")
jsonname=$(awk 'BEGIN {FS="\t"} {print $2}' countNuclei_list.txt | awk "NR==${SGE_TASK_ID}")
posname=$(awk 'BEGIN {FS="\t"} {print $3}' countNuclei_list.txt | awk "NR==${SGE_TASK_ID}")


matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei('$mask','$jsonname', '$posname')"

echo "**** Job ends ****"
date

