#! /bin/bash

/home/meta/anaconda3/envs/metagenomic/bin/bwa mem -t 50 $1 $2 $3 | /home/meta/anaconda3/envs/samtools1.9/bin/samtools view -bS - | /home/meta/anaconda3/envs/samtools1.9/bin/samtools sort - > $4
# filter out the reads that are not mapping to the gene set (Flag=4), 
# filter out the reads for the secondary mapping (Flag=256), 
# filter out the chimeric reads (Flag=2048), 
# filter out the reads with an alignment number less than 2
/home/meta/anaconda3/envs/samtools1.9/bin/samtools view -F 4 -F 256 -F 2048 $4|awk '{if($3!="*") print $3}'|sort| uniq -c|awk 'BEGIN {FS=" ";OFS=","} {print $2,$1}' | awk 'BEGIN {FS=",";OFS=","} {if ($2 > 1) print $1"\t"$2; else print $1"\t0"}'|sed '1i gene\t'$5'' > $6
