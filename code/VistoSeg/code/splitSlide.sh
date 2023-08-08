#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -pe local 8
#$ -o logs/08-10-2022_V11D01-385_Br3942_Round4_splitSlide.txt
#$ -e logs/08-10-2022_V11D01-385_Br3942_Round4_splitSlide.txt 
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
echo "Sample id: 08-10-2022_V11D01-385_Br3942_Round4"
echo "****"

module load matlab/R2019a

toolbox='/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/code/VistoSeg/code'
fname='/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/raw-data/images/CS2/round4/08-10-2022_V11D01-385_Br3942_Round4.tif'


matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide('$fname',0,0,0,0)"

echo "**** Job ends ****"
date
