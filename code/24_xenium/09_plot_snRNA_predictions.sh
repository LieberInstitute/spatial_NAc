#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=50GB
#SBATCH --job-name=09_plot
#SBATCH --output=logs/09_plot.log
#SBATCH --error=logs/09_plot.log


echo "********* Job Starts *********"
date

module load r_nac
Rscript 09_plot_snRNA_predictions.R

echo "********* Job Ends *********"
date
