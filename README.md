# JASPAR UCSC tracks
This repository contains the data and code used to generate the JASPAR UCSC Genome Browser track data hub. For more information visit the [JASPAR website](http://jaspar.genereg.net/genome-tracks/#ucsc_tracks).

## News
01/07/2018 To speed-up TFBS predictions, we switched from [`MEME`](http://meme-suite.org/doc/overview.html) and the [`Perl TFBS`](http://tfbs.genereg.net) package to [`PWMScan`](http://ccg.vital-it.ch/pwmscan).

## Content
The repository is organized as follows:
* The folder `genomes` contains scripts to download and process different genome assemblies
* The folder `profiles` contains the output from the script [`get_profiles.py`](https://github.com/wassermanlab/JASPAR-UCSC-tracks/blob/master/profiles/get_profiles.py), which downloads JASPAR CORE profiles for different taxons
* The script [`scan_sequence.py`](https://github.com/wassermanlab/JASPAR-UCSC-tracks/blob/master/scan_sequence.py) takes as input the `profiles` folder and a nucleotide sequence, in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format) (*e.g.* a genome), and provides TFBS predictions
* The script [`scans2bigBed`](https://github.com/wassermanlab/JASPAR-UCSC-tracks/blob/master/scans2bigBed) creates a [bigBed track file](https://genome.ucsc.edu/goldenPath/help/bigBed.html) from TFBS predictions
* The file [`environment.yml`](https://github.com/wassermanlab/JASPAR-UCSC-tracks/blob/master/environment.yml) contains the conda environment used to generate the genomic tracks for JASPAR 2020 (it installs most of the dependencies described below)

The original scripts used for the publication of [`JASPAR 2018`](https://doi.org/10.1093/nar/gkx1126) have been placed in the folder [`version-1.0`](https://github.com/wassermanlab/JASPAR-UCSC-tracks/tree/master/version-1.0).

## Dependencies
The scripts for creating the genomic tracks require the following dependencies:
* [`GNU parallel`](https://www.gnu.org/software/parallel/)
* [`Python 3.7`](https://www.python.org/download/releases/3.7/) with the [`Biopython`](http://biopython.org) (<1.74), [`NumPy`](http://www.numpy.org), [`pyfaidx`](https://peerj.com/preprints/970/) and [`tqdm`](https://tqdm.github.io) libraries
* [`PWMScan`](http://ccg.vital-it.ch/pwmscan)
* [`UCSC binaries`](http://hgdownload.cse.ucsc.edu/admin/exe/) for standalone command-line use

Note that for running `scan_sequence.py`, only the `Python` dependencies and `PWMScan` are required.

## Installation
Except for `PWMScan`, which has to be [downloaded](https://sourceforge.net/projects/pwmscan/), installed, and appended to your `PATH` manually, the remaining dependencies can be installed through the [`conda`](https://docs.conda.io/en/latest/) package manager:
```
conda env create -f ./environment.yml
```

## Availability
Genomic tracks and TFBS predictions for human and 6 other model organisms are available online:
* [http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2020/](http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2020/)

## Usage
To illustrate the generation of genomic tracks, we provide an example for the [baker's yeast genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2/):
* Download the genome sequence and chromosome sizes (automated in this [script](https://github.com/wassermanlab/JASPAR-UCSC-tracks/blob/master/genomes/sacCer3/sacCer3.sh))
* Scan the genome sequence using **all** [fungi profiles from the JASPAR CORE](http://jaspar.genereg.net/search?q=&collection=CORE&tax_group=fungi)
```
./scan_sequence.py --fasta-file ./genomes/sacCer3/sacCer3.fa --profiles-dir ./profiles/ --output-dir ./tracks/sacCer3/ --threads 4 --latest --taxon fungi
```
For this example, this step should not take longer than a minute. For human (and for other similar genomes), this step should be completed within a few hours (the final amount of time will depend on the number of `--threads` specified).
* Create the genomic track
```
./scans2bigBed -c ./genomes/sacCer3/sacCer3.chrom.sizes -i ./tracks/sacCer3/ -o ./tracks/sacCer3.bb -t 4
```
TFBS predictions from the previous step are merged into a [bigBed track file](https://genome.ucsc.edu/goldenPath/help/bigBed.html). As scores (column 5), we use <i>p</i>-values from `PWMScan` (scaled between 0-1000, where 0 corresponds to <i>p</i>-value = 1 and 1000 to <i>p</i>-value ≤ 10-10). This allows for comparison of prediction confidence across TFBSs. Again, for this example, this step should be completed within a few minutes, while for larger genomes it can take a few hours.

**Important note:** both disk space and memory requirements for larger genomes (*i.e.* danRer11, hg19, hg38 and mm10) are substantial. In these cases, we highly recommend allocating at least 1Tb of disk space and 512Gb of ram.
