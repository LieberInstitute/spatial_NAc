#!/bin/bash
#SBATCH --mem=64G
#SBATCH -n 8
#SBATCH --job-name=NAc-cellranger
#SBATCH -o logs/NAc_cellranger.txt
#SBATCH --array=1-6

# Files to process
# 15c_NAc_SVB - 20c_NAc_SVB 

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load CellRanger
module load cellranger/7.2.0

## List current modules for reproducibility
module list

## Locate file
SAMPLE=$(awk 'BEGIN {FS="\t"} {print $1}' 03_cellranger.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
SAMPLEID=$(awk 'BEGIN {FS="\t"} {print $2}' 03_cellranger.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
## SAMPLE=$(awk "NR==${SGE_TASK_ID}" 03_cellranger.txt)
echo "Processing sample ${SAMPLE}"
date

## Run CellRanger
cellranger count --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/raw-data/FASTQ_reorg/${SAMPLE} \
    --sample=${SAMPLEID} \
    --jobmode=local \
    --localcores=8 \
    --localmem=64

## Move output
echo "Moving data to new location"
date
mkdir -p /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/11_cellranger/
mv ${SAMPLE} /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/11_cellranger/

echo "**** Job ends ****"
date