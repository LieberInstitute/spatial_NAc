#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=25GB
#SBATCH --job-name=Table_S2
#SBATCH --output=logs/Table_S2.log
#SBATCH --error=logs/Table_S2.log


echo "********* Job Starts *********"
date


module load r_nac
Rscript make_supp_table_2.R


echo "********* Job Ends *********"
date
