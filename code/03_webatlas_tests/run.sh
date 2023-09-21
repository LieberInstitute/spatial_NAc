#   First, had to go into 'nextflow.config' and set singularity to true instead
#   of docker. First run of the pipeline needs a huge amount of memory (e.g.
#   ~60G) to build the singularity image from the docker version

#!/bin/bash

#$ -cwd
#$ -N "run"
#$ -o ../../processed-data/03_webatlas_tests/run.log
#$ -e ../../processed-data/03_webatlas_tests/run.log
#$ -pe local 4
#$ -l mf=10G,h_vmem=10G,h_fsize=50G

#SBATCH -q shared
#SBATCH -c 4
#SBATCH --mem=40G
#SBATCH --job-name=run
#SBATCH -o ../../processed-data/03_webatlas_tests/run.log
#SBATCH -e ../../processed-data/03_webatlas_tests/run.log

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

module load nextflow/22.10.7
module load singularity/3.2.1

WORK_DIR=$(git rev-parse --show-toplevel)/processed-data/03_webatlas_tests/work

nextflow run webatlas-pipeline/main.nf \
    -params-file ./config.yaml \
    -entry Full_pipeline \
    -w $WORK_DIR 

echo "**** Job ends ****"
