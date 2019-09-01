wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
tar xvfz chromFa.tar.gz 
rm *_*.fa
cat *.fa > mm10.fa
rm chromFa.tar.gz
faidx mm10.fa -i chromsizes > mm10.chrom.sizes