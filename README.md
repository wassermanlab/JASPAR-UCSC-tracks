# JASPAR UCSC tracks
**TO BE UPDATED**

## News
01/07/2018 The JASPAR UCSC tracks module now uses [`PWMScan`](http://ccg.vital-it.ch/pwmscan), instead of [`MEME`](http://meme-suite.org/doc/overview.html) and the `Perl` [`TFBS`](http://tfbs.genereg.net) package, for speeding-up genome-wide TFBS predictions.

## Content
The repository is organized as follows:
* The `profiles` folder contains the output from `make_files.py`: *i.e.* the JASPAR profiles from the different CORE collections
* The scripts `functions.py`, `jaspar_search.py` and `make_files.py`

The original scripts used for the publication of [`JASPAR 2018`](https://doi.org/10.1093/nar/gkx1126) have been placed in the `version-1.0` folder.

## Dependencies
The scripts for creating the JASPAR UCSC tracks require the following dependencies:
* [`PWMScan`](http://ccg.vital-it.ch/pwmscan)
* [`Python 3.7`](https://www.python.org/download/releases/3.7/) with the [`Biopython`](http://biopython.org), [`NumPy`](http://www.numpy.org), [`pyfaidx`](https://peerj.com/preprints/970/) and [`tqdm`](https://tqdm.github.io) libraries
* [`UCSC binaries`](http://hgdownload.cse.ucsc.edu/admin/exe/) for standalone command-line use

## TFBS prediction availability
Pre-calculated genome-wide TFBS predictions for human and various model organisms are available through [`http://jaspar.genereg.net/genome-tracks/#ucsc_tracks`](http://jaspar.genereg.net/genome-tracks/#ucsc_tracks) and [`http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/`](http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/).

## Usage
To illustrate the generation of the JASPAR genome tracks, we provide an example for [yeast](https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2/):
* Scan the yeast genome using all profiles from fungi
`python scan_sequence.py --fasta-file ./genomes/sacCer3/sacCer3.fa --profiles-dir ./profiles/ --dummy-dir ./tmp/ --output-dir ./scans/sacCer3/ --threads 32 --latest --taxon fungi`
* Create a BED file
For each TF binding profile, the human DNA sequence (in FASTA format) was scanned using the jaspar_search.py script, and matches with a relative score ≥0.8 and with a p-value <0.05 were kept (i.e. TFBS predictions that were not consistent between the TFBS Perl module and FIMO were filtered out).

./jaspar_search.py -f $GENOME_FASTA -j $JASPAR_MATRIX_ID -m $MEME_DIR -o $SCANS_DIR -p $PROFILES_DIR

Create a sorted BED file
TFBS predictions were converted to BED format. As scores (column 5), we used FIMO p-values (scaled between 0-1000, where 0 corresponds to p-value = 1 and 1000 to p-value ≤ 10-10) to allow for comparison of prediction confidence between different profiles.

./fetch_binding_sites.py -i $SCANS_DIR -p $PROFILES_DIR | sort -k1,1 -k2,2n > $BED_FILE

Create a UCSC Genome Browser bigBed track file
Finally, BED files were converted to bigBed format for visualization in the UCSC Genome Browser using bedToBigBed, as distributed within the UCSC binaries for standalone command-line use.

bedToBigBed -type=bed6 -tab -extraIndex=name $BED_FILE $CHROM_SIZES $BIGBED_FILE
*
*
`bash scans2bed.sh genomes/sacCer3/sacCer3.chrom.sizes scans/sacCer3/ > tracks/sacCer3/sacCer3.bed`

