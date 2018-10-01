# JASPAR UCSC tracks
**TO BE UPDATED**

## News
01/07/2018 The JASPAR UCSC tracks module now uses [`PWMScan`](http://ccg.vital-it.ch/pwmscan), instead of [`MEME`](http://meme-suite.org/doc/overview.html) and the `Perl` [`TFBS`](http://tfbs.genereg.net) package, for speeding-up genome-wide TFBS predictions.

## Content
The repository is organized as follows:
* The `files` folder contains the output from `make_files.py`: *i.e.* the JASPAR profiles from the different CORE collections
* The scripts `functions.py`, `jaspar_search.py` and `make_files.py`

The original scripts used for the publication of [`JASPAR 2018`](https://doi.org/10.1093/nar/gkx1126) have been placed in the `version-1.0` folder.

## Dependencies
The scripts for creating the JASPAR UCSC tracks require the following dependencies:
* [`PWMScan`](http://ccg.vital-it.ch/pwmscan)
* [`Python 2.7`](https://www.python.org/download/releases/2.7/) with the [`Biopython`](http://biopython.org), [`NumPy`](http://www.numpy.org) and [`tqdm`](https://tqdm.github.io) libraries
* [`UCSC binaries`](http://hgdownload.cse.ucsc.edu/admin/exe/) for standalone command-line use

## TFBS prediction availability
Pre-calculated genome-wide TFBS predictions for human and various model organisms are available through [`http://jaspar.genereg.net/genome-tracks/#ucsc_tracks`](http://jaspar.genereg.net/genome-tracks/#ucsc_tracks) and [`http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/`](http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/).

## Usage
**TO BE UPDATED**
