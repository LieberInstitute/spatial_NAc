#!/bin/bash
#SBATCH -o logs/submissionLinks-250908_%a.o.txt


cd ./fastqs/

md5sum * >> ../08-raw-checksum.txt

