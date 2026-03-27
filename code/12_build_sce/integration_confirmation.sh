#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=50GB
#SBATCH --job-name=integration_confirmation
#SBATCH --output=logs/integration_confirmation.log
#SBATCH --error=logs/integration_confirmation.log


echo "********* Job Starts *********"
date


module load r_nac
Rscript integration_confirmation.R


echo "********* Job Ends *********"
date
