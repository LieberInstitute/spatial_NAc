#!/bin/bash
#SBATCH
#SBATCH --job-name=run_vs_stitcher
#SBATCH --time=0:15:0
#SBATCH --partition=defq
#SBATCH --mem=20G

module load anaconda 
conda activate vs_stitcher

python vs_stitcher.py > vs_stitcher.log


