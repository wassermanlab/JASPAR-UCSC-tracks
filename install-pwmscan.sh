#!/bin/bash

## JASPAR-UCSC-tracks
## Author: Oriol Fornes
## Contact: oriol@cmmt.ubc.ca

##
## Download and install PWMScan
##

# From https://stackoverflow.com/questions/59895/
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
mkdir -p $DIR/pwmscan
cd $DIR/pwmscan
curl -L -O https://sourceforge.net/projects/pwmscan/files/pwmscan/rel-1.1.9/pwmscan.1.1.9.tar.gz
tar xvfz pwmscan.1.1.9.tar.gz
cd pwmscan
mkdir -p bin
make clean && make cleanbin
make && make install
ln -s ./bin/matrix_prob ../
ln -s ./bin/matrix_scan ../
