#!/bin/bash

CHROM_SIZES_FILE=$1
SCANS_DIR=$2

cut -f 1 $CHROM_SIZES_FILE | LC_COLLATE=C sort -k1,1 -k2,2n |
{
while read i; do zless $SCANS_DIR/*.bed.gz | grep "^$i[[:space:]]" | cut -f 1-4,6,7 | sort -k1,1 -k2,2n | awk -v max=1000 '{if($5>max){print $1"\t"$2"\t"$3"\t"$4"\t"max"\t"$6}else{print $0}}'; done
}
