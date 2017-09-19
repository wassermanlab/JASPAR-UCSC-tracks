# JASPAR UCSC tracks
For the 2018 release of [JASPAR](http://jaspar2018.genereg.net/), we have scanned the human genome for TF binding sites (TFBS) using the [CORE vertebrates TF binding profiles](http://jaspar2018.genereg.net/collection/core/). The predictions are publicly available as UCSC Genome Browser track data hubs:
* [hg19 assembly](http://www.google.com)
* [hg38 assembly](http://www.google.com)

# Usage
We generated custom UCSC Genome Browser track data hubs containing genome-wide TFBS predictions for TF binding profiles in the JASPAR CORE vertebrates collection. Specifically, for each profile, the human genome assemblies [hg19](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/) and [hg38](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/) were scanned in parallel using the [TFBS Perl module](http://tfbs.genereg.net) and [FIMO](http://meme-suite.org/doc/fimo.html), as distributed within the [MEME suite](http://meme-suite.org/meme-software/4.11.2/meme_4.11.2_2.tar.gz) (version 4.11.2). The data are purely computational, and as such not all predicted binding sites will be functional. Please refer to the JASPAR2018 paper for more details.

## Step 1. Reformat JASPAR profiles
For scanning the human genome with the BioPerl TFBS module, we converted profiles to [PWMs](https://en.wikipedia.org/wiki/Position_weight_matrix) using the `jaspar2pfm.py` script.

`./jaspar2pfm.py -b ./files/JASPAR2018_CORE_vertebrates.txt -o $PROFILES_DIR`

For the FIMO scan, profiles were reformatted to [MEME motifs](http://meme-suite.org/doc/meme-format.html) using the `jaspar2meme.py` script.

`./jaspar2meme.py -b ./files/JASPAR2018_CORE_vertebrates.txt -m $MEME_DIR -o $PROFILES_DIR`

## Step 2. Scanning of the human genome
For each TF binding profile, the human DNA sequence (in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format) was scanned using the `jaspar_search.py` script, and matches with a relative score ≥ 0.8 and with a *p*-value < 0.05 were kept (*i.e.* TFBS predictions that were not consistent between the TFBS Perl module and FIMO were filtered out.)

`./jaspar_search.py -f $GENOME_FASTA -j $JASPAR_MATRIX_ID -m $MEME_DIR -o $SCANS_DIR -p $PROFILES_DIR`

## Step 3. Create a sorted BED file
TFBS predictions were converted to [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). As scores (column 5), we used FIMO *p*-values (scaled between 0-1000, where 0 corresponds to *p*-value = 1, 600 to *p*-value = 10<sup>-6</sup>, and 1000 to *p*-value ≤ 10<sup>-10</sup>) to allow for comparison of prediction confidence between different profiles.

`./fetch_binding_sites.py -i $SCANS_DIR -p $PROFILES_DIR | sort -k1,1 -k2,2n > $BED_FILE`

## Step 4. Create a UCSC Genome Browser bigBed track file
Finally, BED files were converted to [bigBed format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1.5) for visualization in the UCSC Genome Browser using [UCSC binaries for standalone command-line use](http://hgdownload.cse.ucsc.edu/admin/exe/). The genome browser translates BED score values into shades of [gray](https://genome.ucsc.edu/FAQ/FAQformat.html#format1).

`bedToBigBed -type=bed6 -tab -extraIndex=name $BED_FILE $CHROM_SIZES $BIGBED_FILE`
