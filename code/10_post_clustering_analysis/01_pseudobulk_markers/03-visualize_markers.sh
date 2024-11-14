#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH -t 1-00:00
#SBATCH --job-name=03-visualize_markers
#SBATCH -o ../../../processed-data/10_post_clustering_analysis/01_pseudobulk_markers/01_precast/logs/03-visualize_markers_rs5_k%a.log
#SBATCH -e ../../../processed-data/10_post_clustering_analysis/01_pseudobulk_markers/01_precast/logs/03-visualize_markers_rs5_k%a.log
#SBATCH --array=3-15%15

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load r_nac
capture_area_path="01_precast/pseudobulk_capture_area/random_start_5"
donor_path="01_precast/pseudobulk_donor/random_start_5"
cluster_col="precast_k${SLURM_ARRAY_TASK_ID}"
spot_path="01_precast/nnSVG_precast/random_start_5/PRECAST_k${SLURM_ARRAY_TASK_ID}_marker_genes.csv"
Rscript 03-visualize_markers.R -c $capture_area_path -d $donor_path -p $cluster_col -s $spot_path

echo "**** Job ends ****"
date
