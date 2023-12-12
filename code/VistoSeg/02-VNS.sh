#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH --job-name=02-VNS
#SBATCH -o ../../processed-data/VistoSeg/logs/02-VNS_%a.log
#SBATCH -e ../../processed-data/VistoSeg/logs/02-VNS_%a.log
#SBATCH --array=2-38%5

num_clusters=5

#   Paths as R code
toolbox_dir="here::here('code', 'VistoSeg', 'code')"
inputs_csv="here::here('processed-data', 'VistoSeg', 'VistoSeg_inputs.csv')"

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load matlab/R2023a
module load conda_R/4.3
module list

#   Evaluate variable values in R and for this array task
toolbox=$(Rscript -e "cat($toolbox_dir)")
image_path=$(Rscript -e "
sample_info = read.csv($inputs_csv)
cat(sample_info[$SLURM_ARRAY_TASK_ID, 'raw_image_path'])
")
sample_id=$(Rscript -e "
sample_info = read.csv($inputs_csv)
cat(sample_info[$SLURM_ARRAY_TASK_ID, 'sample_id'])
")

echo "Running VNS with sample ID ${sample_id}..."

#   Run VNS
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), VNS('$image_path',$num_clusters)"

echo "**** Job ends ****"
date
