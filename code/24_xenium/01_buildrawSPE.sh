#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=50GB
#SBATCH --job-name=01_rawSPE
#SBATCH --output=logs/01_rawSPE.log
#SBATCH --error=logs/01_rawSPE.log


echo "********* Job Starts *********"
date


module load r_nac
Rscript 01_buildrawSPE.R


echo "********* Job Ends *********"
date
