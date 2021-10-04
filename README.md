# JASPAR UCSC tracks
This repository contains the data and code used to generate the JASPAR UCSC Genome Browser track data hub.</br>For more information visit the [JASPAR website](http://jaspar.genereg.net/genome-tracks/#ucsc_tracks).

## News
01/07/2018 To speed-up TFBS predictions, we switched from [MEME](http://meme-suite.org/doc/overview.html) and the [Perl TFBS](http://tfbs.genereg.net) package to [PWMScan](http://ccg.vital-it.ch/pwmscan).

## Content
* The `genomes` folder contains scripts to download and process different genome assemblies
* The `profiles` folder contains the output from the script [`get-profiles.py`](https://github.com/wassermanlab/JASPAR-UCSC-tracks/blob/master/profiles/get_profiles.py), which downloads the JASPAR CORE profiles for different taxons
* The file [`environment.yml`](https://github.com/wassermanlab/JASPAR-UCSC-tracks/blob/master/environment.yml), within the `conda` folder, contains the conda environment used to generate the genomic tracks for JASPAR 2022 (see installation)
* The script [`install-pwmscan.sh`](https://github.com/wassermanlab/JASPAR-UCSC-tracks/blob/master/install-pwmscan.sh) downloads and installs PWMscan and places its binaries in the in the `bin` folder.
* The script [`scan-sequence.py`](https://github.com/wassermanlab/JASPAR-UCSC-tracks/blob/master/scan_sequence.py) takes as its input the `profiles` folder and a nucleotide sequence in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format)</br>(*e.g.* a genome), and outputs TFBS predictions
* The script [`scans2bigBed`](https://github.com/wassermanlab/JASPAR-UCSC-tracks/blob/master/scans2bigBed) creates a [bigBed track file](https://genome.ucsc.edu/goldenPath/help/bigBed.html) from TFBS predictions

The original scripts used for the publication of [JASPAR 2018](https://doi.org/10.1093/nar/gkx1126) have been placed in the folder [`version-1.0`](https://github.com/wassermanlab/JASPAR-UCSC-tracks/tree/master/version-1.0).

## Dependencies
* [Python 3.7](https://www.python.org/download/releases/3.7/) with the following libraries: [Biopython](http://biopython.org) (<1.74), [NumPy](http://www.numpy.org), [pyfaidx](https://peerj.com/preprints/970/) and [tqdm](https://tqdm.github.io)
* [PWMScan](http://ccg.vital-it.ch/pwmscan)
* [UCSC binaries](http://hgdownload.cse.ucsc.edu/admin/exe/) for standalone command-line use

Note that for running `scan_sequence.py`, only the Python dependencies and PWMScan are required.

## Installation
To install PWMScan, execute the script [`install-pwmscan.sh`](https://github.com/wassermanlab/JASPAR-UCSC-tracks/blob/master/install-pwmscan.sh).

The remaining dependencies can be installed through the [conda](https://docs.conda.io/en/latest/) package manager:
```
conda env create -f ./conda/environment.yml
```

## Availability
Genomic tracks and TFBS predictions for human and **seven** other model organisms, covering **11** genome assemblies, are available online:
* [http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2022/](http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2022/)

## Usage
To illustrate how the genomic tracks are generated, we provide an example for the [baker's yeast genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2/):
* Download the genome sequence and chromosome sizes (automated in this [script](https://github.com/wassermanlab/JASPAR-UCSC-tracks/blob/master/genomes/sacCer3/sacCer3.sh))
* Scan the genome sequence using [**all** fungi profiles from the JASPAR CORE](http://jaspar.genereg.net/search?q=&collection=CORE&tax_group=fungi)
```
./scan-sequence.py --fasta-file ./genomes/sacCer3/sacCer3.fa --profiles-dir ./profiles/ \
    --output-dir ./tracks/sacCer3/ --threads 4 --latest --taxon fungi
```
For this example, the scanning step should take no longer than a minute. For human and other similar genomes, this step is usually finished within a few hours (the final amount of time will depend on the number of `--threads` specified).
* Create the genomic track
```
./scans2bigBed -c ./genomes/sacCer3/sacCer3.fa.sizes -i ./tracks/sacCer3/ -o ./tracks/sacCer3.bb -t 4
```
TFBS predictions from the previous step are merged into a [bigBed track file](https://genome.ucsc.edu/goldenPath/help/bigBed.html). In column five, we use as scores the <i>p</i>-values from PWMScan (scaled between 0-1000, where 0 corresponds to <i>p</i>-value = 1 and 1000 to <i>p</i>-value ≤ 10-10). This allows for comparison of prediction confidence across TFBSs. Again, for this example, this step should be completed within a few minutes, while for larger genomes it can take a few hours.

**Important note:** disk space requirements for large genomes (*i.e.* danRer11, hg19, hg38, mm10, and mm39) are substantial. In these cases, we highly recommend allocating at least 1Tb of disk space.