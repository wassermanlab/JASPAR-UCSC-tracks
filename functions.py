import os, sys
import gzip
import shutil

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