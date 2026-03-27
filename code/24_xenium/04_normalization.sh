#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=50GB
#SBATCH --job-name=04_norm
#SBATCH --output=logs/04_norm.log
#SBATCH --error=logs/04_norm.log


echo "********* Job Starts *********"
date

module load r_nac
Rscript 04_normalization.R

echo "********* Job Ends *********"
date
