#!/bin/bash
#
#SBATCH --job-name=DimRed
#SBATCH --output=DimRed.out
#SBATCH --error=DimRed.err
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=2
#
# Mimimum memory required per allocated  CPU
#SBATCH --mem-per-cpu=40G
#
# Send mail to the email address when the job fails
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=robert.phillips@libd.org

cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc

echo "********* Job Starts *********"
date

#load R
module load r_nac

#list modules for reproducibility purposes
module list

#run the Rjob
R CMD BATCH /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/code/12_build_sce/05_DimensionalityReduction.R

echo "********* Job Ends *********"
date
