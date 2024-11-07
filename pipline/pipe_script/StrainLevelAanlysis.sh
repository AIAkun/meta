#! /bin/bash

%%bash
 
mkdir -p strainphlan_output/consensus_marker
 
for s in $(ls -d ${input}/*/); do
    sample=$(basename ${s})
    sample2markers.py -i metaphlan_output/${sample}/${sample}.sam.bz2 \
                      -f bz2 \
                      -o strainphlan_output/consensus_marker \
                      -n 4 \
                      --clades t__SGB4933
done


%%bash
 
mkdir strainphlan_output/db_markers
 
extract_markers.py -c t__SGB4933 \
                   -o strainphlan_output/db_markers


%%bash
 
mkdir strainphlan_output/t__SGB4933
 
strainphlan -s strainphlan_output/consensus_marker/*.pkl \
            -m strainphlan_output/db_markers/t__SGB4933.fna \
            -o strainphlan_output/t__SGB4933 \
            -n 4 \
            -c t__SGB4933 \
            --marker_in_n_samples 4 \
            --sample_with_n_markers 30 \
            --sample_with_n_markers_after_filt 30

%%bash
 
for var_annot in bs_ls sam; do
    echo "Plotting StrainPhlAn4 t__SGB4933 tree with ${var_annot}"
    add_metadata_tree.py -t strainphlan_output/t__SGB4933/RAxML_bestTree.t__SGB4933.StrainPhlAn4.tre \
                         -f /home/hutlab_public/Tutorials/metaphlan4/input/metadata_samples.tsv \
                         -m ${var_annot} \
                         --string_to_remove _filtered
 
    plot_tree_graphlan.py -t strainphlan_output/t__SGB4933/RAxML_bestTree.t__SGB4933.StrainPhlAn4.tre.metadata \
                          -m ${var_annot}
    
    mv strainphlan_output/t__SGB4933/RAxML_bestTree.t__SGB4933.StrainPhlAn4.tre.metadata.png \
       strainphlan_output/t__SGB4933/RAxML_bestTree.t__SGB4933.StrainPhlAn4.tre.metadata_${var_annot}.png
done                               