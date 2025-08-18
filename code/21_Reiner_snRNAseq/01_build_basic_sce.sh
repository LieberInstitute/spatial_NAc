#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=200GB
#SBATCH --job-name=01_buildbasicSCE
#SBATCH --output=../../processed-data/21_Reiner_snRNAseq/logs/01_build_basic_sce.log
#SBATCH --error=../../processed-data/21_Reiner_snRNAseq/logs/01_build_basic_sce.log


echo "********* Job Starts *********"
date


module load r_nac
Rscript 01_build_basic_sce.R

echo "********* Job Ends *********"
date
