#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH -t 2-00:00
#SBATCH --job-name=01-pseudobulk
#SBATCH -o ../../../processed-data/10_post_clustering_analysis/01_pseudobulk_markers/01_precast/logs/01-pseudobulk_BayesSpace_k%a.log
#SBATCH -e ../../../processed-data/10_post_clustering_analysis/01_pseudobulk_markers/01_precast/logs/01-pseudobulk_BayesSpace_k%a.log
#SBATCH --array=2-28%15

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load r_nac
cluster_dir="05_harmony_BayesSpace/05-BayesSpace_k_search/BayesSpace_harmony_k${SLURM_ARRAY_TASK_ID}"
cluster_file="clusters.csv"
agg_level="sample_id_original"
out_path="02_BayesSpace/pseudobulk_capture_area"
Rscript 01-compute_pseudobulk.R -d $cluster_dir -f $cluster_file -a $agg_level -o $out_path

echo "**** Job ends ****"
date
