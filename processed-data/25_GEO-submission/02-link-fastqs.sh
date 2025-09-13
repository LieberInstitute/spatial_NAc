#!/bin/bash
#SBATCH --array=1-38
#SBATCH -o logs/submissionLinks-250908_%a.o.txt

# /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/01_spaceranger_reorg
# /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/25_GEO-submission/

NAMESPACE=$(awk '{print $1}' 00-list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
SAMPLE=$(awk '{print $2}' 00-list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
SRLOC=$(awk '{print $3}' 00-list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")


cd ./fastqs/

for file in ../../../raw-data/FASTQ_reorg/${SAMPLE}/*.fastq.gz; do ln -rsf $file ./${NAMESPACE}-$(basename $file) ; done

