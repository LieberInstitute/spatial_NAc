#!/bin/bash
#$ -cwd
#$ -l mem_free=32G,h_vmem=32G,h_fsize=100G
#$ -N BayesSpace_k_search
#$ -o ../../processed-data/05_harmony_BayesSpace/logs/03-BayesSpace_k_search.log
#$ -e ../../processed-data/05_harmony_BayesSpace/logs/03-BayesSpace_k_search.log
#$ -t 2-28
#$ -tc 20

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

Rscript 03-BayesSpace_k_search.R

echo "**** Job ends ****"
date
