#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=50GB
#SBATCH --job-name=08_transfer
#SBATCH --output=logs/08_transfer.log
#SBATCH --error=logs/08_transfer.log


echo "********* Job Starts *********"
date

module load r_nac
Rscript 08_snRNA_xenium_transfer.R

echo "********* Job Ends *********"
date
