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
* [`GNU parallel`](https://www.gnu.org/software/parallel/)
* [`Python 3.7`](https://www.python.org/download/releases/3.7/) with the [`Biopython`](http://biopython.org), [`NumPy`](http://www.numpy.org), [`pyfaidx`](https://peerj.com/preprints/970/) and [`tqdm`](https://tqdm.github.io) libraries
* [`PWMScan`](http://ccg.vital-it.ch/pwmscan)
* [`UCSC binaries`](http://hgdownload.cse.ucsc.edu/admin/exe/) for standalone command-line use

## TFBS prediction availability
Predictions for human and 6 other model organisms are available online:
* [`http://jaspar.genereg.net/genome-tracks/#ucsc_tracks`](http://jaspar.genereg.net/genome-tracks/#ucsc_tracks)
* [`http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/`](http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/).

## Usage
To illustrate the generation of the JASPAR genome tracks, we provide an example for the [yeast genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2/):
1) Download the genome and chromosome sizes
We provide a [script](https://github.com/wassermanlab/JASPAR-UCSC-tracks/blob/master/genomes/sacCer3/sacCer3.sh) that automatically downloads and uncompresseses the genome of yeast, and calculates the size of each chromosome. 
2) Scan the genome using [all fungi profiles from the JASPAR CORE](http://jaspar.genereg.net/search?q=&collection=CORE&tax_group=fungi)
For yeast, this step should not take longer than a minute. For human (and other genomes of similar size), this step should be completed within a few hours (the total amount of time will depend on the number of `--threads` specified).
```
./scan_sequence.py --fasta-file ./genomes/sacCer3/sacCer3.fa --profiles-dir ./profiles/ --output-dir ./scans/sacCer3/ --threads 16 --latest --taxon fungi
```
3) Create a [genome browser bigBed track](https://genome.ucsc.edu/goldenPath/help/bigBed.html)
TFBS predictions from the Python script [`scan_sequence.py`](https://github.com/wassermanlab/JASPAR-UCSC-tracks/blob/master/scan_sequence.py) are merged into a bigBed track file using the bash script [`scans2bigBed`](https://github.com/wassermanlab/JASPAR-UCSC-tracks/blob/master/scans2bigBed). As scores (column 5), we use <i>p</i>-values from PWMScan (scaled between 0-1000, where 0 corresponds to <i>p</i>-value = 1 and 1000 to <i>p</i>-value ≤ 10-10) to allow for comparison of prediction confidence across TFBSs. Again, for yeast this step should finish within a few minutes, while for larger genomes it could take a few hours.
```
./scans2bigBed -c ./genomes/sacCer3/sacCer3.chrom.sizes -i ./scans/sacCer3/ -o ./tracks/sacCer3.bed -t 16
```