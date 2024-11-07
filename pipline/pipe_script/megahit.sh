#! /bin/bash

/home/meta/anaconda3/envs/metagenomic/bin/megahit -1 $1 -2 $2 -o $3 --out-prefix $4 -t 50
/home/meta/anaconda3/envs/metagenomic/bin/seqkit seq -m 500 $5 --remove-gaps > $6
sed -i 's/>/>'$4'_/g' $6


# 加quast
/home/meta/anaconda3/envs/graphlan/bin/quast -o $3'/quast_eval' $6


