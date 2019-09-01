wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar xvfz chromFa.tar.gz 
rm *_*.fa
cat *.fa > hg19.fa
rm chromFa.tar.gz
faidx hg19.fa -i chromsizes > hg19.chrom.sizes