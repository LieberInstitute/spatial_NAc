#!/bin/bash
#
#SBATCH --job-name=11_neuron_clustering
#SBATCH --output=logs/11_neuron_clustering.log
#SBATCH --error=logs/11_neuron_clustering.log
#SBATCH -p shared
#SBATCH --mem=25G

echo "********* Job Starts *********"
date

#load R
module load r_nac
Rscript 11_neuronsubcluster_clustering.R

echo "********* Job Ends *********"
date
