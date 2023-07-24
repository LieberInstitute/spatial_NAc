#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -pe local 4
#$ -N preprocess_and_harmony
#$ -o ../../processed-data/05_harmony_BayesSpace/logs/02-preprocess_and_harmony.log
#$ -e ../../processed-data/05_harmony_BayesSpace/logs/02-preprocess_and_harmony.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load conda_R/4.3
module list

Rscript 02-preprocess_and_harmony.R

echo "**** Job ends ****"
date
