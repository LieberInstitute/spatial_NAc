#!/bin/bash
#$ -cwd
#$ -N "build_spe"
#$ -o ../../processed-data/05_harmony_BayesSpace/logs/01-build_spe.log
#$ -e ../../processed-data/05_harmony_BayesSpace/logs/01-build_spe.log
#$ -l mf=20G,h_vmem=20G,h_fsize=50G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.3
Rscript 01-build_spe.R

echo "**** Job ends ****"
date
