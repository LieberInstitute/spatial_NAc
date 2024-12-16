#!/bin/bash

#SBATCH -p shared
#SBATCH -c 4
#SBATCH --time=8:0:0
#SBATCH --mem=30G
#SBATCH --job-name=move_files

cp -r /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/10_precast/nnSVG_default /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/07_spatial_domains/01_precast/
cp -r /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/10_precast/nnSVG_precast /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/07_spatial_domains/01_precast/
 