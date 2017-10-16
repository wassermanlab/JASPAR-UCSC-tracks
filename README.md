# JASPAR UCSC tracks
For the 2018 release of [JASPAR](http://jaspar2018.genereg.net/), we have performed TFBS predictions on the human genome ([hg19](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/) and [hg38](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/) assemblies) using the [CORE vertebrates profiles](http://jaspar2018.genereg.net/collection/core/), which are publicly available as a UCSC Genome Browser track data hub:
* [UCSC tracks for the hg19 assembly](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hubUrl=http://expdata.cmmt.ubc.ca/JASPAR/UCSC_tracks/hub.txt)
* [UCSC tracks for the hg38 assembly](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubUrl=http://expdata.cmmt.ubc.ca/JASPAR/UCSC_tracks/hub.txt)

For other popular organisms, including *Arabidopsis thaliana* ([araTha1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001735.3/)), *Drosophila melanogaster* ([dm6](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001215.4/)), mouse ([mm10](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.20/)), worm ([ce10](https://www.ncbi.nlm.nih.gov/assembly/GCF_000002985.5/)), yeast ([sacCer3](https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2/)), and zebrafish ([danRer10](https://www.ncbi.nlm.nih.gov/assembly/GCF_000002035.5/)), genome tracks are currently in preparation. We will make these tracks available as soon as they are computed. Please do not hesitate to contact us for requesting genome tracks for your popular organisms of interest. 

## Dependencies
The scripts for creating the JASPAR UCSC tracks require the following dependencies:
* [MEME](http://meme-suite.org/doc/overview.html) suite (≥4.12.0)
* Perl (>5.14) with the [BioPerl](http://bioperl.org) and [TFBS](http://tfbs.genereg.net) packages
* Python (>2.7) with the *Biopython* (≥1.65) and *numpy* (≥1.8.2) libraries
* [UCSC binaries](http://hgdownload.cse.ucsc.edu/admin/exe/) for standalone command-line use

## Usage
We generated a custom UCSC Genome Browser track data hub containing genome-wide binding site predictions for TF profiles in the JASPAR CORE vertebrates collection. For each profile, the human genome assemblies [hg19](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/) and [hg38](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/) were scanned in parallel using the [TFBS](http://tfbs.genereg.net) Perl module and [FIMO](http://meme-suite.org/doc/fimo.html), as distributed within the [MEME](http://meme-suite.org/meme-software/4.11.2/meme_4.11.2_2.tar.gz) suite. Please refer to the JASPAR 2018 manuscript for more details.

### Reformat JASPAR profiles
For scanning the human genome with the BioPerl TFBS module, we converted profiles to [PWMs](https://en.wikipedia.org/wiki/Position_weight_matrix) using the `jaspar2pfm.py` script.

`./jaspar2pfm.py -b ./files/JASPAR2018_CORE_vertebrates.txt -o $PROFILES_DIR`

For the FIMO scan, profiles were reformatted to [MEME motifs](http://meme-suite.org/doc/meme-format.html) using the `jaspar2meme.py` script.

`./jaspar2meme.py -b ./files/JASPAR2018_CORE_vertebrates.txt -m $MEME_DIR -o $PROFILES_DIR`

Individual PFMs in JASPAR/MEME format can also be retrieved from the JASPAR downloads page [here](http://jaspar.genereg.net/downloads/).

### Scanning of the human genome
For each TF binding profile, the human DNA sequence (in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format) was scanned using the `jaspar_search.py` script, and matches with a relative score ≥0.8 and with a *p*-value <0.05 were kept (*i.e.* TFBS predictions that were not consistent between the TFBS Perl module and FIMO were filtered out).

`./jaspar_search.py -f $GENOME_FASTA -j $JASPAR_MATRIX_ID -m $MEME_DIR -o $SCANS_DIR -p $PROFILES_DIR`

### Create a sorted BED file
TFBS predictions were converted to [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). As scores (column 5), we used FIMO *p*-values (scaled between 0-1000, where 0 corresponds to *p*-value = 1 and 1000 to *p*-value ≤ 10<sup>-10</sup>) to allow for comparison of prediction confidence between different profiles.

`./fetch_binding_sites.py -i $SCANS_DIR -p $PROFILES_DIR | sort -k1,1 -k2,2n > $BED_FILE`

### Create a UCSC Genome Browser bigBed track file
Finally, BED files were converted to [bigBed format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1.5) for visualization in the UCSC Genome Browser using [bedToBigBed](http://hgdownload.cse.ucsc.edu/admin/exe/), as distributed within the UCSC binaries for standalone command-line use.

`bedToBigBed -type=bed6 -tab -extraIndex=name $BED_FILE $CHROM_SIZES $BIGBED_FILE`
