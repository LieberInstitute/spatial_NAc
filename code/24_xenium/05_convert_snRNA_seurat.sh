#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=50GB
#SBATCH --job-name=05_snRNA_seurat
#SBATCH --output=logs/05_snRNA_seurat.log
#SBATCH --error=logs/05_snRNA_seurat.log


echo "********* Job Starts *********"
date

module load r_nac
Rscript 05_convert_snRNA_seurat.R

echo "********* Job Ends *********"
date
