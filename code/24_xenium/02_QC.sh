#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=50GB
#SBATCH --job-name=02_QC
#SBATCH --output=logs/02_QC.log
#SBATCH --error=logs/02_QC.log


echo "********* Job Starts *********"
date


module load r_nac
Rscript 02_QC.R


echo "********* Job Ends *********"
date
