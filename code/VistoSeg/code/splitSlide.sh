#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -pe local 7
#$ -o logs/06-02-2022_V11U08-083_Br8492_Round1_splitSlide.txt
#$ -e logs/06-02-2022_V11U08-083_Br8492_Round1_splitSlide.txt 
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
echo "Sample id: 06-02-2022_V11U08-083_Br8492_Round1"
echo "****"

module load matlab/R2019a

toolbox='/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/code/VistoSeg/code'
fname='//dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/raw-data/images/CS2/round1/06-02-2022_V11U08-083_Br8492_Round1.tif'


matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide('$fname',0,0,0,0)"

echo "**** Job ends ****"
date
