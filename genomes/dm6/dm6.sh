wget http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz 
gunzip dm6.fa.gz 
cat dm6.fa | perl -e 'while(<>){chomp;if(substr($_,0,1) eq ">"){$chr=substr($_,1);}system("echo \"$_\" >> $chr.fa");}'
rm *_*.fa
faidx dm6.fa -i chromsizes > dm6.chrom.sizes