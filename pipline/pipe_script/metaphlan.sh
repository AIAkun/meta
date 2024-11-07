#! /bin/bash

zcat $1 $2|/home/meta/anaconda3/envs/metagenomic/bin/metaphlan --input_type fastq --bowtie2out $3 --output_file $4 --nproc 50 
