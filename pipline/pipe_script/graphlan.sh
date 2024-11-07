#! /bin/bash

#/home/meta/anaconda3/envs/graphlan/bin/export2graphlan.py -i $1 --tree $2 --annotation $3 --most_abundant 100 --abundance_threshold 1 --least_biomarkers 10 --annotations 5,6 --external_annotations 7 --min_clade_size 1
/home/meta/anaconda3/envs/graphlan/bin/export2graphlan.py --skip_rows 1 -i $1 --tree $2 --annotation $3 --most_abundant 50 --abundance_threshold 0.01 --least_biomarkers 10 --external_annotations 7 --min_clade_size 1 --annotation_legend_font_size 10
# 绘图
/home/meta/anaconda3/envs/graphlan/bin/graphlan_annotate.py --annot $3 $2 $4
/home/meta/anaconda3/envs/graphlan/bin/graphlan.py --dpi 300 $4 $5 
