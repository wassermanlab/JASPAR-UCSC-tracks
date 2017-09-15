import os, sys, re
from Bio.motifs.jaspar.db import JASPAR5
import gzip
import shutil

def get_jaspar_profiles(matrix_id=None):
    """
    @input:
    matrix_id {str} e.g. MA0002.2
    file_name {filename} e.g. pfm_vertebrates.txt
    @return:
    {motif} or {list} of {motifs} if matrix_id is None

    """

    JASPAR_DB_HOST = 'vm5.cmmt.ubc.ca'
    JASPAR_DB_NAME = 'JASPAR_2016'
    JASPAR_DB_USER = 'jaspar_r'
    JASPAR_DB_PASS = ''
    
    jaspar_db = JASPAR5(host=JASPAR_DB_HOST, name=JASPAR_DB_NAME, user=JASPAR_DB_USER, password=JASPAR_DB_PASS)

    if matrix_id is None:
        return jaspar_db.fetch_motifs(collection='CORE')
    
    return jaspar_db.fetch_motif_by_id(matrix_id)

def parse_file(file_name, gz=False):
    """
    This function parses any file and yields lines one by one.
    
    @input:
    file_name {string}
    @return:
    line {string}

    """
 
    if os.path.exists(file_name):
        # Initialize #
        f = None
        # Open file handle #
        if gz:
            try: f = gzip.open(file_name, "rt")
            except: raise ValueError("Could not open file %s" % file_name)
        else:
            try: f = open(file_name, "rt")
            except: raise ValueError("Could not open file %s" % file_name)
        # For each line... #
        for line in f:
            yield line.strip("\n")
        f.close()
    else:
        raise ValueError("File %s does not exist!" % file_name)

def parse_fasta_file(file_name, gz=False, clean=True):
    """
    This function parses any FASTA file and yields sequences one by one
    in the form header, sequence.

    @input:
    file_name {string}
    @return:
    line {list} header, sequence

    """

    # Initialize #
    header = ""
    sequence = ""
    # For each line... #
    for line in parse_file(file_name, gz):
        if len(line) == 0: continue
        if line.startswith("#"): continue
        if line.startswith(">"):
            if sequence != "":
                if clean:
                    sequence = re.sub("\W|\d", "X", sequence)
                yield header, sequence
            m = re.search("^>(.+)", line)
            header = m.group(1)
            sequence = ""
        else:
            sequence += line.upper()
    if clean:
        sequence = re.sub("\W|\d", "X", sequence)

    yield header, sequence        

def parse_tsv_file(file_name, gz=False):
    """
    This function parses any TSV file and yields lines as a list.

    @input:
    file_name {string}
    @return: {list}
    line {list}

    """

    # For each line... #
    for line in parse_file(file_name, gz):
        line = line.split("\t")
        yield line

def write(file_name=None, line=None):
    """
    This function adds a {line} to file (or to stdout if no file
    was provided).

    @input:
    file_name {string}
    line {string}

    """
    if file_name is None: sys.stdout.write("%s\n" % line)
    else:
        with open(file_name, "a") as out_file:
            out_file.write("%s\n" % line)

def compress(input_file, output_file):
    """
    This function compresses a file using gzip.

    @input:
    input_file {filename}
    output_file {filename}

    """

    with open(input_file, "rb") as in_file, gzip.open(output_file, "wb") as out_file:
        shutil.copyfileobj(in_file, out_file)
