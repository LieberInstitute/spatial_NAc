#!/bin/bash
#
#SBATCH --mail-type=END
#SBATCH --mail-user=Robert.Phillips@libd.org
#SBATCH --job-name=01_NMF
#SBATCH --output=logs/01_NMF.log
#SBATCH --error=logs/01_NMF.log
#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=15G
#SBATCH --time=10-00:00:00

echo "********* Job Starts *********"
date

#load R
module load r_nac
Rscript 01_NMF_Factor.R

echo "********* Job Ends *********"
date
