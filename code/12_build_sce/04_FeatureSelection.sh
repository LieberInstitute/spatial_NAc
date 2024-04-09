#!/bin/bash

#SBATCH --job-name=04_featureslxn
#SBATCH --output=logs/04_featureslxn.log 
#SBATCH --error=logs/04_featureslxn.log
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=90G
#SBATCH --mail-type=END
#SBATCH --mail-user=Robert.Phillips@libd.org

echo "********* Job Starts *********"
date

module load r_nac
Rscript 04_FeatureSelection.R

echo "********* Job Ends *********"
date




