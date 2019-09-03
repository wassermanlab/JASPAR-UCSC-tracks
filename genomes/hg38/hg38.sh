wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
tar xvfz hg38.chromFa.tar.gz 
find chroms -type f | grep -v "_" | perl -e 'while(<>){chomp;system("mv $_ .");}'
rm -r chroms/ hg38.chromFa.tar.gz
cat *.fa > hg38.fa
faidx hg38.fa -i chromsizes > hg38.chrom.sizes
