nucl_nonerude=$1
db=$2
virus=$3
blast_best='/home/meta/metagenomic/pipeline/pipe_script/blast_best.py'
merge='/home/meta/metagenomic/pipeline/pipe_script/merge2.py'
merge_virus='/home/meta/metagenomic/pipeline/pipe_script/merge_virus.py'
blastn -db $db -query $nucl_nonerude -perc_identity 80 -qcov_hsp_perc 80 -max_target_seqs 5 -num_threads 16 -evalue 1e-5 -outfmt 6 -out $virus/blast_virus.txt
cat $virus/blast_virus.txt | cut -f 1 | sort -u > $virus/virus.map.id
seqkit grep -f $virus/map.id $nucl_nonerude > $virus/virus.fasta
seqkit grep --invert-match -f $virus/map.id $nucl_nonerude > $virus/virus_host.fasta

input=$virus/blast_virus.txt
sed -i "1i query\tsubject\tidentity_percent\talignment length\tmismatches\tgap opens\tq_start\tq_end\ts_start\ts_end\tevalue\tbit score" $input
python $blast_best $input  $virus/virus.best_blast.xls
annotation='/home/meta/metagenomic/database/virus_blastn/virus-blastn_anno.txt'
python $merge --nt $annotation --i $virus/best_blast.xls --o $virus/blast_anno_virus.xls

# python /home/meta/metagenomic/softs/CrisprOpenDB/CL_Interface.py -i $virus/map.fasta -n 50 > $virus/vh.txt

# cat $virus/vh.txt | grep -v 'No hits found. Sorry'|grep '(' |sed 's/(//g'|sed 's/)//g'| sed 's/,/\t/g'|sed 's/'"'"/' ''/g' |sed 's/ //g'> $virus/vh_result.txt

# sed -i "1i query\thost_bacteria\tlevel" $virus/vh_result.txt

# python $merge_virus --nt $virus/vh_result.txt --i $virus/blast_anno_virus.xls --o $virus/final_anno.xls