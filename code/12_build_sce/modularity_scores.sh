#!/bin/bash
#
#SBATCH --job-name=mod_scores
#SBATCH --output=logs/mod_scores.log
#SBATCH --error=logs/mod_scores.log
#SBATCH -p shared
#SBATCH --mem=55G

echo "********* Job Starts *********"
date

#load R
module load r_nac
Rscript modularity_scores.R

echo "********* Job Ends *********"
date
