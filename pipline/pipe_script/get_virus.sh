#! /bin/bash

/home/meta/anaconda3/envs/metagenomic/bin/bowtie2 -x $1 -1 $2 -2 $3 -S $4 -p 50 2>$5
# -f 4 means retain the reads that are not mapping to the reference genome 
/home/meta/anaconda3/envs/samtools1.9/bin/samtools fastq -@ 50 -F 4 -1 $6 -2 $7 -s $8 $4
