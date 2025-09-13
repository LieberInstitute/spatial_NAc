#!/bin/bash
#SBATCH --array=1-20
#SBATCH -o logs/submissionLinks-250908_%a.o.txt

# /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/raw-data/FASTQ_reorg
# /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/25_GEO-submission/snRNA-data

SAMPLE=$(awk '{print $1}' 10-snRNA-list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
SRLOC=$(awk '{print $2}' 10-snRNA-list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
NAMESPACE=$(awk '{print $3}' 10-snRNA-list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")


mkdir -p /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/25_GEO-submission/snRNA-data/processed

cd ./snRNA-data/processed/

ln -rsf ../../../11_cellranger/${SAMPLE}/outs/raw_feature_bc_matrix/barcodes.tsv.gz ./${NAMESPACE}-barcodes.tsv.gz
ln -rsf ../../../11_cellranger/${SAMPLE}/outs/raw_feature_bc_matrix/features.tsv.gz ./${NAMESPACE}-features.tsv.gz
ln -rsf ../../../11_cellranger/${SAMPLE}/outs/raw_feature_bc_matrix/matrix.mtx.gz ./${NAMESPACE}-matrix.mtx.gz
