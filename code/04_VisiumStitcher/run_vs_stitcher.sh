#!/bin/bash

#$ -cwd
#$ -N "run_visium_stitcher"
#$ -o ../../processed-data/04_VisiumStitcher/logs/01-vs_stitcher.log
#$ -e ../../processed-data/04_VisiumStitcher/logs/01-vs_stitcher.log
#$ -l mf=40G,h_vmem=40G,h_fsize=50G

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH --job-name=run_vs_stitcher
#SBATCH -o ../../processed-data/04_VisiumStitcher/logs/01-vs_stitcher.log
#SBATCH -e ../../processed-data/04_VisiumStitcher/logs/01-vs_stitcher.log


source ~/.bashrc
module load anaconda 
conda activate vs_stitcher

python vs_stitcher.py


