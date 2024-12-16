#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -t 1-00:00
#SBATCH --job-name=03-visualize_markers
#SBATCH -o ../../../processed-data/10_post_clustering_analysis/01_pseudobulk_markers/02_BayesSpace/logs/03-visualize_markers_k%a.log
#SBATCH -e ../../../processed-data/10_post_clustering_analysis/01_pseudobulk_markers/02_BayesSpace/logs/03-visualize_markers_k%a.log
#SBATCH --array=3-28%15

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load r_nac
capture_area_path="02_BayesSpace/pseudobulk_capture_area"
donor_path="NA"
cluster_col="BayesSpace_harmony_k${SLURM_ARRAY_TASK_ID}"
spot_path="NA"
Rscript 03-visualize_markers.R -c $capture_area_path -d $donor_path -p $cluster_col -s $spot_path

echo "**** Job ends ****"
date
