grep "@" dat/chromosome_chunks.fasta.fq | awk -F '.' '{ print $1 }' > dat/chr_tax.txt
grep "@" dat/plasmid_chunks.fasta.fq | awk -F '.' '{ print $1 }' > dat/pla_tax.txt

cat dat/pla_tax.txt dat/chr_tax.txt > dat/train.tax.txt