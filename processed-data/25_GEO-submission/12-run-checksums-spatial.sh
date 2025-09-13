#!/bin/bash
#SBATCH -o logs/submissionLinks-250908_%a.o.txt


cd ./spatial-data/

md5sum * >> ../07-processed-checksum.txt

