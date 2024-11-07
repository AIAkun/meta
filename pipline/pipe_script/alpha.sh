#! /bin/bash


mpa_path=/home/meta/anaconda3/envs/metagenomic/lib/python3.7/site-packages/metaphlan/utils
output=$1/alpha_diversity
mkdir ${output}
 
for metric in richness shannon simpson gini; do
    ${mpa_path}/calculate_diversity.R \
        -f $2 \
        -o ${output} \
        -d alpha \
        -p alpha \
        -m ${metric}
done

# ls metaphlan_output/alpha_diversity/
# # df_meta = pd.read_csv('/home/hutlab_public/Tutorials/metaphlan4/input/metadata_samples.tsv', header=0, sep='\t')
# df_alpha = df_alpha.read_csv('metaphlan_output/alpha_diversity/alpha_shannon.tsv', sep='\t', header=0, index_col=0)
# df_alpha = df_alpha.reset_index().rename(columns={'index':'sampleID'}).merge(df_meta[['sampleID', 'bs_ls' ]], on='sampleID', how='left')
 
# plt.figure(figsize=(6, 5), dpi=150)
 
# ax = sns.boxplot(data=df_alpha, y="diversity_shannon", x='bs_ls', hue='bs_ls')
# ax.set_title('Shannon-Index')
 
# plt.savefig('metaphlan_output/alpha_diversity/alpha_shannon.png', dpi=150)
