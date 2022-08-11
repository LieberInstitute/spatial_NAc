#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -pe local 10
#$ -o logs/06-02-2022_V11U08-082_Br8325_Round1_A1_VNS.txt
#$ -e logs/06-02-2022_V11U08-082_Br8325_Round1_A1_VNS.txt 
#$ -m e
#$ -M heenadivecha@gmail.com
#$ -t 1
#$ -tc 2

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

matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), VNS('$fname',5)"

echo "**** Job ends ****"
date



