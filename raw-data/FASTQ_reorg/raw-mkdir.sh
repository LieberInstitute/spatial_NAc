#!/bin/bash
#SBATCH --array=1-3

SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" raw-mkdir.txt)

mkdir -p ./logs/
mkdir -p /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/raw-data/FASTQ_reorg/${SAMPLE}/

mv slurm-*.out ./logs