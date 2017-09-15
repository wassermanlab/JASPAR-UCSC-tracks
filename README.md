# 1. Reformat 
* Downloaded repository under `./scripts/`
* Configure `Parameters` and `Paths` in `./scripts/config.ini`

# 2. Download files
## PDB:
* Go to the [Advanced Search](http://www.rcsb.org/pdb/search/advSearch.do?search=new) tool
* Click on `Choose a Query Type` and select `All/Experimental Type/Molecule Type`
* Leave `Experimental Method` as `Ignore` and select `Protein/NA complex` for `Molecule Type`
* Click on `Submit Query`
* Click on `Download Files`
* Check options `PDB` and `Biological Assemblies` (uncheck the rest)
* Select `uncompressed` for `Compresion Type`
* Click on `Launch Download Application` and download `download_rcsb.jnlp`
* Execute `download_rcsb.jnlp` application and download PDB files under `./rcsb/`
## PBM:
* Go to "http://cisbp.ccbr.utoronto.ca/entireDownload.php".
* Download "E-scores", "TF Information" and "PWMs" and unzip downloaded files under "./cisbp/".
* Go to "http://cisbp.ccbr.utoronto.ca/bulk.php" and click on "Download MySQL tables".
* Unzip downloaded file and tables under "./cisbp/".
## UniProt:
* Go to "http://www.uniprot.org/uniprot/".
* Search for "proteome:(reference:yes AND taxonomy:"Eukaryota [2759]")".
* Click on "Download" and then "Go".
* As of March 2017, this downloads 878 eukaryotic reference proteomes (134,224 sequences from SwissProt and 12,693,240 sequences from TrEMBL).
* Gunzip downloaded file as "uniprot.fasta" under "./uniprot/".
* Format "uniprot.fasta" to be searched by blast: i.e. "./src/ncbi-blast-2.2.31+/bin/makeblastdb -in ./uniprot/uniprot.fasta -dbtype prot"
* Go to "http://www.uniprot.org/docs/speclist.txt".
* Right click and "Save Page As..." > "speclist.txt" under "./uniprot/".

# 3. Identify TF/DNA complex structures
`python ./scripts/tfinder.py -f ./cisbp/cisbp_1.02.tf_families.sql -o ./tfinder/ -p ./cisbp/cisbp_1.02.proteins.sql -r ./rcsb/ -t ./cisbp/cisbp_1.02.tfs.sql -u ./uniprot/uniprot.fasta -v`

# 4. Parse PDB data
`python ./scripts/parsers/pdb.py -o ./pdb/ -r ./rcsb/ -t ./tfinder/tfs.txt -v`

# 5. Parse PBM data
`python scripts/parsers/pbm.py -e ./cisbp/Escores.txt -f ./cisbp/cisbp_1.02.tf_families.sql -o ./pbm/ -m ./cisbp/cisbp_1.02.motifs.sql --pdb=./pdb/ --pwm=./cisbp/pwms/ -s ./cisbp/cisbp_1.02.motif_sources.sql -t ./cisbp/cisbp_1.02.tfs.sql -u uniprot/uniprot.fasta`
