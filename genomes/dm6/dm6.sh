wget http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz 
gunzip dm6.fa.gz 
faidx -x dm6.fa
rm chr*_*.fa
cat chr*.fa > dm6.fa
faidx dm6.fa -i chromsizes > dm6.chrom.sizes
