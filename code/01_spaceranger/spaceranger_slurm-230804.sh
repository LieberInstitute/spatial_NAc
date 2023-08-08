#!/bin/bash
#SBATCH --mem=80G
#SBATCH -n 8
#SBATCH --job-name=spaceranger
#SBATCH -o logs/spaceranger_slurm_230804o.txt
#SBATCH -e logs/spaceranger_slurm_230804e.txt
#SBATCH --array=1-16

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load SpaceRanger
module load spaceranger/2.1.0

## List current modules for reproducibility
module list

## Locate file
SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" samples_2023-08-04.txt)
echo "Processing sample ${SAMPLE}"
date

## Get slide and area
SLIDE=$(echo ${SAMPLE} | cut -d "_" -f 1)
CAPTUREAREA=$(echo ${SAMPLE} | cut -d "_" -f 2)
SAM=$(paste <(echo ${SLIDE}) <(echo "-") <(echo ${CAPTUREAREA}) -d '')
echo "Slide: ${SLIDE}, capture area: ${CAPTUREAREA}"

## Find FASTQ file path
FASTQPATH=$(ls -d /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/raw-data/FASTQ_reorg/${SAMPLE}/)

## Hank from 10x Genomics recommended setting this environment
export NUMBA_NUM_THREADS=1

## Run SpaceRanger
spaceranger count \
    --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=${FASTQPATH}\
    --image=/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/Images/VistoSeg/Capture_areas/newrounds/${SAMPLE}.tif \
    --slide=${SLIDE} \
    --area=${CAPTUREAREA} \
    --loupe-alignment=/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/Images/loupe-alignment/newrounds/${SAM}.json \
    --jobmode=local \
    --localcores=8 \
    --localmem=64

## Move output
echo "Moving results to new location"
date
mkdir -p /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/01_spaceranger_reorg/${SAMPLE}
mv ${SAMPLE} /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/01_spaceranger_reorg/

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
