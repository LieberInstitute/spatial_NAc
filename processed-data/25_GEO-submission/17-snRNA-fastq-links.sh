#!/bin/bash
#SBATCH --array=1-20
#SBATCH -o logs/submissionLinks-250908_%a.o.txt

# /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/25_GEO-submission/
# /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/25_GEO-submission/snRNA-data/fastq/

NAMESPACE=$(awk '{print $1}' 10-snRNA-list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
SAMPLE=$(awk '{print $2}' 10-snRNA-list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
SRLOC=$(awk '{print $3}' 10-snRNA-list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")


cd ./snRNA-data/fastq/

ln -rsf ../../../../raw-data/FASTQ_reorg/${SAMPLE}/*.fastq.gz ./

# for file in ../../../../raw-data/FASTQ_reorg/${SAMPLE}/*.fastq.gz; do ln -rsf $file ./${NAMESPACE}-$(basename $file) ; done

