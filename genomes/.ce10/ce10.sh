#!/usr/bin/env bash

wget http://hgdownload.cse.ucsc.edu/goldenPath/ce10/bigZips/chromFa.tar.gz
tar xvfz chromFa.tar.gz
rm chromFa.tar.gz
cat *.fa > ce10.fa
faidx ce10.fa -i chromsizes > ce10.chrom.sizes
