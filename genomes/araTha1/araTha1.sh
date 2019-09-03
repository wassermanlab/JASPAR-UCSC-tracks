wget https://genome-test.gi.ucsc.edu/~hiram/hubs/Plants/araTha1/araTha1.2bit
twoBitToFa araTha1.2bit araTha1.fa
faidx araTha1.fa -i chromsizes > araTha1.chrom.sizes
