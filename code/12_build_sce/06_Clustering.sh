#!/bin/bash
#
#SBATCH --job-name=06_clustering
#SBATCH --output=logs/06_clustering.log
#SBATCH --error=logs/06_clustering.log
#SBATCH -p shared
#SBATCH --mem=25G

echo "********* Job Starts *********"
date

#load R
module load r_nac
Rscript 06_Clustering.R 

echo "********* Job Ends *********"
date
