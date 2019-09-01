wget http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz
tar xvfz chromFa.tar.gz 
cat chr*.fa > sacCer3.fa
rm chromFa.tar.gz
faidx sacCer3.fa -i chromsizes > sacCer3.chrom.sizes
