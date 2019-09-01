wget http://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.fa.gz
gunzip danRer11.fa.gz 
cat danRer11.fa | perl -e 'while(<>){chomp;if(substr($_,0,1) eq ">"){$chr=substr($_,1);}system("echo \"$_\" >> $chr.fa");}'
rm *_*.fa
faidx danRer11.fa -i chromsizes > danRer11.chrom.sizes