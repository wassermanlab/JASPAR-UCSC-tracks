wget http://hgdownload.soe.ucsc.edu/goldenPath/danRer10/bigZips/danRer10.fa.gz
gunzip danRer10.fa.gz 
cat danRer10.fa | perl -e 'while(<>){chomp;if(substr($_,0,1) eq ">"){$chr=substr($_,1);}system("echo \"$_\" >> $chr.fa");}'
rm *_*.fa
faidx danRer10.fa -i chromsizes > danRer10.chrom.sizes