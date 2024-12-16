#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=80G
#SBATCH --time=24:0:0
#SBATCH --job-name=BayesSpace_k_search
#SBATCH -o ../../processed-data/05_harmony_BayesSpace/logs/05-BayesSpace_k_search_%a.log
#SBATCH -e ../../processed-data/05_harmony_BayesSpace/logs/05-BayesSpace_k_search_%a.log
#SBATCH --array=4-28%15

if [[ ! -z $SLURMD_NODENAME ]]; then
    job_id=$SLURM_JOB_ID
    job_name=$SLURM_JOB_NAME
    node_name=$SLURMD_NODENAME
else
    job_id=$JOB_ID
    job_name=$JOB_NAME
    node_name=$HOSTNAME
fi

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Node name: ${node_name}"

## List current modules for reproducibility
module load r_nac
module list

Rscript 05-BayesSpace_k_search.R

echo "**** Job ends ****"
date
