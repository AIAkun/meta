#! /bin/bash

/home/meta/anaconda3/envs/metagenomic/bin/fastp -i $1 -o $2 -I $3 -O $4 -w 50 -h $5 -j $6

# 加fastqc
