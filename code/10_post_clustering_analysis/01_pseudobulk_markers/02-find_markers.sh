#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH -t 2-00:00
#SBATCH --job-name=02-find_markers
#SBATCH -o ../../../processed-data/10_post_clustering_analysis/01_pseudobulk_markers/02_BayesSpace/logs/02-find_markers_BayesSpace_k%a.log
#SBATCH -e ../../../processed-data/10_post_clustering_analysis/01_pseudobulk_markers/02_BayesSpace/logs/02-find_markers_BayesSpace_k%a.log
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
cluster_col="BayesSpace_harmony_k${SLURM_ARRAY_TASK_ID}"
pseudo_path="02_BayesSpace/pseudobulk_capture_area"
agg_level="sample_id_original"
Rscript 02-find_markers.R -c $cluster_col -p $pseudo_path -a $agg_level

echo "**** Job ends ****"
date
