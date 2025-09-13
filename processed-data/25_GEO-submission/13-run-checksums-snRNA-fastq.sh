#!/bin/bash
#SBATCH -o logs/submissionLinks-250909_%a.o.txt


cd ./snRNA-data/fastq/

md5sum * >> ../../15-raw-checksum-snRNA.txt

