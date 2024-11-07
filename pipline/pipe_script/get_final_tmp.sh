output=$1
level=$2
##转换为biom格式
kraken-biom  --max D --min $level $output/bracken/kreport/*.bracken.S.kreport  -o $output/$level.biom
biom convert  -i $output/$level.biom  --to-tsv --header-key taxonomy -o $output/$level.count.tsv.all.tmp
 
##删除第一行
sed -i '1d' $output/$level.count.tsv.all.tmp

cat $output/$level.count.tsv.all.tmp| grep -v 'k__Archaea'|grep -v 'k__Eukaryota'|grep -v 'k__Bacteria' > $output/$level.count.tsv.Viruses.tmp
cat $output/$level.count.tsv.all.tmp| grep -v 'k__Viruses'| grep -v 'k__Eukaryota'|grep -v 'k__Bacteria' > $output/$level.count.tsv.Archaea.tmp
cat $output/$level.count.tsv.all.tmp| grep -v 'k__Viruses'| grep -v 'k__Archaea'|grep -v 'k__Bacteria' > $output/$level.count.tsv.Eukaryota.tmp
cat $output/$level.count.tsv.all.tmp| grep -v 'k__Viruses'| grep -v 'k__Archaea'|grep -v 'k__Eukaryota' > $output/$level.count.tsv.Bacteria.tmp
