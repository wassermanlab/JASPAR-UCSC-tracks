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
For yeast, this step should not take longer than a minute. For human (and other genomes of similar size), this step be completed within a few hours (total time depends on the number of threads specified).
`python scan_sequence.py --fasta-file ./genomes/sacCer3/sacCer3.fa --profiles-dir ./profiles/ --dummy-dir ./tmp/ --output-dir ./scans/sacCer3/ --threads 32 --latest --taxon fungi`
* Create a sorted BED file
TFBS predictions from `scan_sequence.py` are converted to BED format using the script `scans2bed.sh`. As scores (column 5), we use p-values from PWMScan (scaled between 0-1000, where 0 corresponds to p-value = 1 and 1000 to p-value ≤ 10-10) to allow for comparison of prediction confidence between different profiles. Again, for yeast this step should finish in a few minutes, while for larger genomes it can take hours.
`bash scans2bed.sh genomes/sacCer3/sacCer3.chrom.sizes scans/sacCer3/ > tracks/sacCer3/sacCer3.bed`
* Create a bigBed track file
Finally, BED files were converted to bigBed format for visualization in the UCSC Genome Browser using bedToBigBed, as distributed within the UCSC binaries for standalone command-line use. As in the previous steps, the time consumed by this step will increase with genome size.
`bedToBigBed -type=bed6 -tab -extraIndex=name tracks/sacCer3/sacCer3.bed genomes/sacCer3/sacCer3.chrom.sizes tracks/sacCer3/sacCer3.bb`