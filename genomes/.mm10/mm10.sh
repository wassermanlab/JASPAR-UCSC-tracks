#!/usr/bin/env bash

wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
tar xvfz chromFa.tar.gz 
rm chr*_*.fa chromFa.tar.gz
cat chr*.fa > mm10.fa
faidx mm10.fa -i chromsizes > mm10.chrom.sizes
