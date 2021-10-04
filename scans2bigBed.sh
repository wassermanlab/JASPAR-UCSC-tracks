#!/usr/bin/env bash

./scans2bigBed -c ./genomes/hg38/hg38.fa.sizes -i ./tracks/hg38/ -d ./ -o ./tracks/hg38.bb -t 4
./scans2bigBed -c ./genomes/mm10/mm10.fa.sizes -i ./tracks/mm10/ -d ./ -o ./tracks/mm10.bb -t 4
./scans2bigBed -c ./genomes/mm39/mm39.fa.sizes -i ./tracks/mm39/ -d ./ -o ./tracks/mm39.bb -t 4