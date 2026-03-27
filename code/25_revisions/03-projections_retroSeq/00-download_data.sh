#!/bin/bash

#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=10:0:0
#SBATCH --mem=80G
#SBATCH --job-name=00-mch_download
#SBATCH -o ../../../processed-data/25_revisions/03-projections_retroSeq/00_download_data.log
#SBATCH -e ../../../processed-data/25_revisions/03-projections_retroSeq/00_download_data.log

cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/25_revisions/03-projections_retroSeq
wget -O  CEMBA.epiretro.mcds.tar.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE230nnn/GSE230782/suppl/GSE230782%5FCEMBA%2Eepiretro%2Emcds%2Etar%2Egz"
gunzip CEMBA.epiretro.mcds.tar.gz
tar -xvf CEMBA.epiretro.mcds.tar

