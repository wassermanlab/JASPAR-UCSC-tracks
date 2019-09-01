wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz
tar xvfz chromFa.tar.gz 
rm *_*.fa
cat *.fa > mm9.fa
rm chromFa.tar.gz
faidx mm9.fa -i chromsizes > mm9.chrom.sizes