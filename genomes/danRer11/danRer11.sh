wget http://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.fa.gz
gunzip danRer11.fa.gz 
faidx -x danRer11.fa
rm chr*_*.fa
cat chr*.fa > danRer11.fa
faidx danRer11.fa -i chromsizes > danRer11.chrom.sizes
