# JASPAR UCSC tracks
For the 2018 release of JASPAR, we have performed TFBS predictions on the human genome using the CORE vertebrates TF binding profiles which are publicly available as UCSC Genome Browser track data hubs:
* [hg19](http://www.google.com)
* [hg38](http://www.google.com)

# Usage
We generated custom UCSC Genome Browser track data hubs containing genome-wide TFBS predictions for TF binding profiles in the JASPAR CORE vertebrates collection. Specifically, for each profile, the human genome assemblies [hg19](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/) and [hg38](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/) were scanned in parallel using the [TFBS Perl module](http://tfbs.genereg.net) and [FIMO](http://meme-suite.org/doc/fimo.html), as distributed within the [MEME suite](http://meme-suite.org/meme-software/4.11.2/meme_4.11.2_2.tar.gz) (version 4.11.2).

## Step 1. Reformat JASPAR profiles
For scanning the human genome with the BioPerl TFBS module, we converted profiles to [PWM](https://en.wikipedia.org/wiki/Position_weight_matrix)s using the `jaspar2pfm.py` script.

`./jaspar2pfm.py -b ./files/JASPAR2018_CORE_vertebrates.txt -o $PROFILES_DIR`

For the FIMO scan, profiles were reformatted to [MEME motifs](http://meme-suite.org/doc/meme-format.html) using the `jaspar2meme.py` script.

`./jaspar2meme.py -b ./files/JASPAR2018_CORE_vertebrates.txt -m $MEME_DIR -o $PROFILES_DIR`

## Step 2. Scanning of the human genome
For each TF binding profile, the human genome was scanned and matches with a relative score ≥ 0.8 and with a <i>p</i>-value < 0.05 were kept using the `jaspar_search.py` script (<i>i.e.</i> TFBS predictions that were not consistent between the TFBS Perl module and FIMO were filtered out.)

`./jaspar_search.py -f $GENOME_FASTA -j $JASPAR_MATRIX_ID -m $MEME_DIR -o $SCANS_DIR -p $PROFILES_DIR`

## Step 3. Create a sorted BED file

 The remaining TFBS predictions were converted to genome tracks and colored according to their FIMO p-value (scaled between 0-1000, where 0 corresponds to a p-value of 1 and 1000 to a p-value ≤ 10-10) to allow for comparison of prediction confidence between different profiles. The tracks are collected as a data hub that can be visualized in the UCSC Genome Browser or downloaded for custom analysis (URL). Code and data used to create the UCSC tracks are available at https://github.com/oriolfornes/jaspar-ucsc-tracks. Individual matches for each TF binding profile on each chromosome of the human genome (hg19 and hg38 genome assemblies) are available at $URL.

`./fetch_binding_sites.py -i $SCANS_DIR -p $PROFILES_DIR | sort -k1,1 -k2,2n > $BED_FILE`

## Step 4. Create a UCSC Genome Browser bigBed track file
`bedToBigBed -type=bed6 -tab $BED_FILE $CHROM_SIZES $BIGBED_FILE`
