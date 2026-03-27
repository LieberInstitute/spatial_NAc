#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=50GB
#SBATCH --job-name=06_seurat
#SBATCH --output=logs/06_visium_seurat.log
#SBATCH --error=logs/06_visium_seurat.log


echo "********* Job Starts *********"
date

module load r_nac
Rscript 06_convert_visium_seurat.R

echo "********* Job Ends *********"
date
