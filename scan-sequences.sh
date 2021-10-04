#!/usr/bin/env bash

./scan-sequence.py --threads 4 --latest \
    -o ./tracks/sacCer3/ --taxon fungi ./genomes/sacCer3/sacCer3.fa ./profiles/
./scan-sequence.py --threads 4 --latest \
    -o ./tracks/dm6/ --taxon insects ./genomes/dm6/dm6.fa ./profiles/
./scan-sequence.py --threads 4 --latest \
    -o ./tracks/ce10/ --taxon nematodes ./genomes/ce10/ce10.fa ./profiles/
./scan-sequence.py --threads 4 --latest \
    -o ./tracks/ce11/ --taxon nematodes ./genomes/ce11/ce11.fa ./profiles/
./scan-sequence.py --threads 4 --latest \
    -o ./tracks/araTha1/ --taxon plants ./genomes/araTha1/araTha1.fa ./profiles/
./scan-sequence.py --threads 4 --latest \
    -o ./tracks/ci3/ --taxon urochordates ./genomes/ci3/ci3.fa ./profiles/
./scan-sequence.py --threads 4 --latest \
    -o ./tracks/danRer11/ --taxon vertebrates ./genomes/danRer11/danRer11.fa \
    ./profiles/
./scan-sequence.py --threads 4 --latest \
    -o ./tracks/hg19/ --taxon vertebrates ./genomes/hg19/hg19.fa ./profiles/
./scan-sequence.py --threads 4 --latest \
    -o ./tracks/hg38/ --taxon vertebrates ./genomes/hg38/hg38.fa ./profiles/
./scan-sequence.py --threads 4 --latest \
    -o ./tracks/mm10/ --taxon vertebrates ./genomes/mm10/mm10.fa ./profiles/
./scan-sequence.py --threads 4 --latest \
    -o ./tracks/mm39/ --taxon vertebrates ./genomes/mm39/mm39.fa ./profiles/