#!/bin/bash

SCANS_DIR=$1
OUTPUT_BED=$2

UNSORTED_BED=.unsorted.bed

zless $SCANS_DIR/*.bed.gz | cut -f 1-4,6,7 > $UNSORTED_BED
LC_ALL=C sort --parallel=16 --buffer-size=1G -k1,1 -k2,2n $UNSORTED_BED > $OUTPUT_BED
