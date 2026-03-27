#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=50GB
#SBATCH --job-name=03_discard
#SBATCH --output=logs/03_discard.log
#SBATCH --error=logs/03_discard.log


echo "********* Job Starts *********"
date


module load r_nac
Rscript 03_discard.R


echo "********* Job Ends *********"
date
