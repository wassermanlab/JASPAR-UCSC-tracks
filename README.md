# UCSC Genome Browser tracks for JASPAR 2018

## Step 1. Reformat JASPAR profiles
### to PFMs:
`./jaspar2pfm.py -b ./files/JASPAR2018_CORE_vertebrates.txt -o $PROFILES_DIR`
### to MEME motifs:
`./jaspar2meme.py -b ./files/JASPAR2018_CORE_vertebrates.txt -m $MEME_DIR -o $PROFILES_DIR`

## Step 2. Scanning of the human genome
`./jaspar_search.py -f $GENOME_FASTA -j $JASPAR_MATRIX_ID -m $MEME_DIR -o $SCANS_DIR -p $PROFILES_DIR`

## Step 3. Create a BED file
`./fetch_binding_sites.py -i $SCANS_DIR -p $PROFILES_DIR > $BED_FILE`

## Step 4. Create a UCSC Genome Browser bigBed track file
`./fetch_binding_sites.py -i $SCANS_DIR -p $PROFILES_DIR > $BED_FILE`
