#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=8:0:0
#SBATCH --mem=80G
#SBATCH --mail-user=pravich2@jh.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=03-visualize_sn_nmf
#SBATCH -o ../../../processed-data/16_transfer_learning/02_target_projections/logs/03-visualize_sn_nmf_rat_case_control_human_NAc.log
#SBATCH -e ../../../processed-data/16_transfer_learning/02_target_projections/logs/03-visualize_sn_nmf_rat_case_control_human_NAc.log

set -e
dataName="human_NAc"
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load r_nac
Rscript 03-visualize_sn_nmf.R -d $dataName

echo "**** Job ends ****"
date