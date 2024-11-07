kraken2 --db /home/meta/metagenomic/database/k2_virus/  --threads 56  \
--report $1 \
--output $2  \
--paired $3 $4

cat 03-kraken2/all/V350092088_L03_100_CONTROL.output | awk '$1 == "C"' | cut -f 2 