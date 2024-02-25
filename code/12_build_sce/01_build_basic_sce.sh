#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=40GB
#SBATCH --job-name=01_buildbasicSCE
#SBATCH --output=logs/01_build_basic_sce.log
#SBATCH --error=logs/01_build_basic_sce.log


echo "********* Job Starts *********"
date


module load r_nac
Rscript 01_build_basic_sce.R


echo "********* Job Ends *********"
date
