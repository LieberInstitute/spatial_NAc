#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=250GB
#SBATCH --job-name=Fig1_Heatmap
#SBATCH --output=logs/Fig1_Heatmap.log
#SBATCH --error=logs/Fig1_Heatmap.log


echo "********* Job Starts *********"
date


module load r_nac
Rscript Fig1_Heatmap.R


echo "********* Job Ends *********"
date
