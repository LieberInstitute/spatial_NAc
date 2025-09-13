#!/bin/bash
#SBATCH -o logs/submissionLinks-snRNA-250909_%a.o.txt


cd ./snRNA-data/processed/

md5sum * >> ../../16-processed-checksum-snRNA

