#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=50GB
#SBATCH --job-name=07_xenium_seurat
#SBATCH --output=logs/07_xenium_seurat.log
#SBATCH --error=logs/07_xenium_seurat.log


echo "********* Job Starts *********"
date

module load r_nac
Rscript 07_convert_xenium_seurat.R

echo "********* Job Ends *********"
date
