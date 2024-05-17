#!/bin/bash
#
#SBATCH --job-name=purity
#SBATCH --output=logs/purity.log
#SBATCH --error=logs/purity.log
#SBATCH -p shared
#SBATCH --mem=70G

echo "********* Job Starts *********"
date

#load R
module load r_nac
Rscript neighborhood_purity.R

echo "********* Job Ends *********"
date
