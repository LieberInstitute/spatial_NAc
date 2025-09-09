#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=250GB
#SBATCH --job-name=ppp1r1b_featureplot
#SBATCH --output=logs/ppp1r1b_featureplot.log
#SBATCH --error=logs/ppp1r1b_featureplot.log


echo "********* Job Starts *********"
date


module load r_nac
Rscript PPP1R1B_FeaturePlot.R


echo "********* Job Ends *********"
date
