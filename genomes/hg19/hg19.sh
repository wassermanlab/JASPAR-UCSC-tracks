wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar xvfz chromFa.tar.gz 
rm chr*_*.fa chromFa.tar.gz
cat chr*.fa > hg19.fa
faidx hg19.fa -i chromsizes > hg19.chrom.sizes
