#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -pe local 8
#$ -o logs/07-27-2022_V11U23-405_Br6522_Round3_splitSlide.txt
#$ -e logs/07-27-2022_V11U23-405_Br6522_Round3_splitSlide.txt 
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
echo "Sample id: 07-27-2022_V11U23-405_Br6522_Round3"
echo "****"

module load matlab/R2019a

toolbox='/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/code/VistoSeg/code'
fname='/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/raw-data/images/CS2/round3/07-27-2022_V11U23-405_Br6522_Round3.tif'


matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide('$fname',0,0,0,0)"

echo "**** Job ends ****"
date
