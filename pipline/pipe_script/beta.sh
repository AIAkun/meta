#! /bin/bash


mpa_path=/home/meta/anaconda3/envs/metagenomic/lib/python3.7/site-packages/metaphlan/utils
output=$1/beta_diversity
mkdir -p ${output}
 
for metric in bray-curtis jaccard weighted-unifrac unweighted-unifrac clr aitchison; do
    ${mpa_path}/calculate_diversity.R \
        -f $2 \
        -o ${output} \
        -d beta \
        -p beta \
        -t /home/meta/anaconda3/envs/metagenomic/lib/python3.7/site-packages/metaphlan/utils//mpa_vJun23_CHOCOPhlAnSGB_202403.nwk \
        -m ${metric}
done

# import pandas as pd
# import numpy as np
# from sklearn.manifold import TSNE
# import seaborn as sns; sns.set(style='ticks')
# import matplotlib
# import matplotlib.pyplot as plt
 
 
# # load beta diversity distance-matric (Bray-Curtis)
# bc_distmat = pd.read_csv('metaphlan_output/beta_diversity/beta_bray-curtis.tsv', sep='\t', header=0, index_col=0)
 
# # calculate dimensionality reduction using t-SNE (t-Stochastic Neighbour Embedding)
# coords_tsne = TSNE(n_components=2, metric="precomputed", random_state=42).fit_transform(bc_distmat)
 
# bc_points = pd.DataFrame(bc_distmat.columns, columns=['sampleID'])
# bc_points['tSNE1'] = [i[0] for i in coords_tsne]
# bc_points['tSNE2'] = [i[1] for i in coords_tsne]
# bc_points.to_csv('metaphlan_output/beta_diversity/beta_bray-curtis_tSNE.tsv', sep='\t', index=False)
# # df_meta = pd.read_csv('metaphlan_output/metadata_samples.tsv', header=0, sep='\t')
 
# bc_points = bc_points.merge(df_meta[['sampleID', 'bs_ls' ]], on='sampleID', how='left')
 
# plt.figure(figsize=(6, 5), dpi=150)
 
# ax = sns.scatterplot(data=bc_points, x=bc_points.columns[1], y=bc_points.columns[2], hue='bs_ls', alpha=0.8)
# ax.set_title('t-SNE Bray-Curtis')
 
# plt.savefig('metaphlan_output/beta_diversity/beta_bray-curtis_tSNE.png')
