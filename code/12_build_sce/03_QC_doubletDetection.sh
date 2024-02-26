#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=03_QC
#SBATCH --output=logs/03_QC_doubletdetection.log 
#SBATCH --error=logs/03_QC_doubletdetection.log 
#SBATCH --mem=50G
#SBATCH --mail-type=END
#SBATCH --mail-user=Robert.Phillips@libd.org

echo "********* Job Starts *********"
date

module load r_nac
Rscript 03_QC_doubletDetection.R 

echo "********* Job Ends *********"
date




