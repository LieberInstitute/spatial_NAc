#!/bin/bash
#SBATCH --array=1-38
#SBATCH -o logs/submissionLinks-250908_%a.o.txt

# /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/01_spaceranger_reorg
# /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/25_GEO-submission/

NAMESPACE=$(awk '{print $1}' 00-list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
SAMPLE=$(awk '{print $2}' 00-list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
SRLOC=$(awk '{print $3}' 00-list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")




mkdir -p /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/25_GEO-submission/spatial-data

cd ./spatial-data/

ln -rsf ../../${SRLOC}/${SAMPLE}/outs/spatial/aligned_fiducials.jpg ./${NAMESPACE}-aligned_fiducials.jpg
ln -rsf ../../${SRLOC}/${SAMPLE}/outs/spatial/detected_tissue_image.jpg ./${NAMESPACE}-detected_tissue_image.jpg
ln -rsf ../../${SRLOC}/${SAMPLE}/outs/spatial/scalefactors_json.json ./${NAMESPACE}-scalefactors_json.json
ln -rsf ../../${SRLOC}/${SAMPLE}/outs/spatial/tissue_hires_image.png ./${NAMESPACE}-tissue_hires_image.png
ln -rsf ../../${SRLOC}/${SAMPLE}/outs/spatial/tissue_lowres_image.png ./${NAMESPACE}-tissue_lowres_image.png
ln -rsf ../../${SRLOC}/${SAMPLE}/outs/spatial/tissue_positions_list.csv ./${NAMESPACE}-tissue_positions_list.csv
ln -rsf ../../${SRLOC}/${SAMPLE}/outs/raw_feature_bc_matrix/barcodes.tsv.gz ./${NAMESPACE}_barcodes.tsv.gz
ln -rsf ../../${SRLOC}/${SAMPLE}/outs/raw_feature_bc_matrix/features.tsv.gz ./${NAMESPACE}_features.tsv.gz
ln -rsf ../../${SRLOC}/${SAMPLE}/outs/raw_feature_bc_matrix/matrix.mtx.gz ./${NAMESPACE}_matrix.mtx.gz