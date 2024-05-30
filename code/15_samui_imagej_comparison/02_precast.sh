#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=30G
#SBATCH --job-name=02_precast
#SBATCH -c 1
#SBATCH -t 2-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-8%20

## Define loops and appropriately subset each variable for the array task ID
all_final_step=(imagej samui)
final_step=${all_final_step[$(( $SLURM_ARRAY_TASK_ID / 4 % 2 ))]}

all_k=(2 10 19 28)
k=${all_k[$(( $SLURM_ARRAY_TASK_ID / 1 % 4 ))]}

## Explicitly pipe script output to a log
log_path=../../processed-data/15_samui_imagej_comparison/logs/02_precast_${final_step}_${k}_${SLURM_ARRAY_TASK_ID}.log

{
set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.4

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 02_precast.R --final_step ${final_step} --k ${k}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.2
## available from http://research.libd.org/slurmjobs/

