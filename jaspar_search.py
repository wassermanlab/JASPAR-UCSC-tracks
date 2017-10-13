#!/usr/bin/env python2.7
import os
from Bio import motifs
from numpy import log10 as log
import optparse
import shutil
import subprocess

# Import my functions #
import functions

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("./%prog -f <fasta_file> -j <jaspar_matrix_id> -m <meme_dir> -p <profiles_dir> [--dummy=<dummy_dir> -o <output_dir> -r <results_dir> --pv-thresh=<p_value_thresh> --rs-thresh=<rel_score_thresh>]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="<dummy_dir>")
    parser.add_option("-f", action="store", type="string", dest="fasta_file", help="FASTA file (e.g. chr1.fa)", metavar="<fasta_file>")
    parser.add_option("-j", action="store", type="string", dest="matrix_id", help="JASPAR matrix ID (e.g. MA0002.2)", metavar="<jaspar_matrix_id>")
    parser.add_option("-m", action="store", type="string", dest="meme_dir", help="Full path to MEME bin directory (i.e. where all MEME executables are located; e.g. $MEME_PATH/bin)", metavar="<meme_dir>")
    parser.add_option("-o", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="<output_dir>")
    parser.add_option("-p", action="store", type="string", dest="profiles_dir", help="Profiles directory (from jaspar2meme.py and jaspar2pfm.py)", metavar="<profiles_dir>")
    parser.add_option("-r", default="./", action="store", type="string", dest="results_dir", help="Results directory (output directory from applying %prog on a previous JASPAR release; default = None)", metavar="<results_dir>")
    parser.add_option("--pv-thresh", default=0.05, action="store", type="float", dest="p_value_thresh", help="P-value threshold (default = 0.05)", metavar="<p_value_thresh>")
    parser.add_option("--rs-thresh", default=0.8, action="store", type="float", dest="rel_score_thresh", help="Relative score threshold (default = 0.8)", metavar="<rel_score_thresh>")

    (options, args) = parser.parse_args()

    if options.fasta_file is None or options.matrix_id is None or options.meme_dir is None or options.profiles_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def jaspar_search(matrix_file, fasta_file, thresh=0.75):
    """
    @input:
    matrix_file {str} e.g. MA0002.2.pfm
    fasta_file {filename} e.g. chr1.fa
    thresh {float} e.g. 0.75
    @yield:

    """   

    # Initialize #
    jaspar_search = os.path.join(os.path.abspath(os.path.dirname(__file__)), "jaspar_search.pl")
    if thresh <= 1: thresh = int(thresh * 100)
    # Exec process #
    process = subprocess.check_output([jaspar_search, "-f", fasta_file, "-m", matrix_file, "-th", str(thresh) + '%'], stderr=subprocess.STDOUT)
    # For each line... #
    for line in process.split("\n"):
        # If match... #
        if line.startswith("chr"):
            #printf $ofh "$chrom\t%d\t%d\t%s\t%4.0f\n"
            try:
                chromosome, start, end, strand, relative_score = line.split("\t")
                yield chromosome, int(start), int(end), strand, float(relative_score) 
            except: continue

def fimo_search(meme_dir, meme_file, fasta_file, thresh=0.05):
    """
    @input:
    meme_dir {directory} i.e. bin directory where all MEME executables are located
    meme_file {filename} e.g. MA0002.2.meme
    fasta_file {filename} e.g. chr1.fa
    thresh {float} e.g. 0.05
    @yield:

    """

    # Initialize #
    fimo = os.path.join(os.path.abspath(os.path.dirname(__file__)), "fimo")
    # Exec process #
    process = subprocess.check_output([os.path.join(meme_dir, "fimo"), "--thresh", str(thresh), "--text", meme_file, fasta_file], stderr=subprocess.STDOUT)
    # For each line... #
    for line in process.split("\n"):
        if line.startswith("#"): continue
        line = line.split("\t")
        if len(line) < 7: continue
        #pattern name, sequence name, start, stop, strand, score, p-value, q-value, matched sequence
        try: yield line[1], int(line[2]), int(line[3]), line[4], float(line[6])
        except: continue

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Create output directory #
    if not os.path.exists(os.path.abspath(options.output_dir)):
        os.makedirs(os.path.abspath(options.output_dir))

    # Initialize #
    dummy_fasta = os.path.join(os.path.abspath(options.dummy_dir), "%s.fa" % os.getpid())
    pfm_file = os.path.join(os.path.abspath(options.profiles_dir), "%s.pfm" % options.matrix_id)

    # Load profile #
    with open(pfm_file) as f:
        profile = motifs.read(f, "jaspar")

    # For each header, sequence... #
    for header, sequence in functions.parse_fasta_file(os.path.abspath(options.fasta_file)):
        # Initialize #
        dummy_file = os.path.join(os.path.abspath(options.dummy_dir), "%s.%s.tab" % (options.matrix_id, header))
        output_file = os.path.join(os.path.abspath(options.output_dir), "%s.%s.tab.gz" % (options.matrix_id, header))
        # Remove dummy file if exist #
        if os.path.exists(dummy_file): os.remove(dummy_file)
        # Create an empty dummy file #
        os.sytem("touch %s" % dummy_file)
        # Skip if output file already exists #
        if os.path.exists(output_file): continue
        # If output file exists in a previous JASPAR release... #
        if options.results_dir is not None:
            # Initialize #
            results_file = os.path.join(os.path.abspath(options.results_dir), "%s.%s.tab.gz" % (options.matrix_id, header))
            # If results file exists... #
            if os.path.exists(results_file):
                # Copy #
                shutil.copy(results_file, output_file)
                continue
        # Get chunks #
        m = 10000;
        n = len(profile) - 1
        chunks = [sequence[i:i+m] for i in range(0, len(sequence), m - n)]
        # If last chunk is too small, merge it to the previous #
        if len(chunks[-1]) <= n: last_chunk = chunks.pop(-1); chunks[-1] += last_chunk
        # For each chunk... #
        for i in range(len(chunks)):
            # Initialize #
            relative_scores = {}
            chunk_start = i * (m - n)
            # Remove dummy FASTA file if exist #
            if os.path.exists(dummy_fasta): os.remove(dummy_fasta)
            # Create dummy FASTA file #
            functions.write(dummy_fasta, ">%s\n%s" % (header, chunks[i]))
            # For each jaspar match... #
            for chromosome, start, end, strand, relative_score in jaspar_search(pfm_file, dummy_fasta, options.rel_score_thresh):
                # Add to relative scores #
                relative_scores.setdefault((start + chunk_start, strand), int(relative_score * 1000))
            # For each fimo match... #
            for chromosome, start, end, strand, p_value in fimo_search(os.path.abspath(options.meme_dir), os.path.join(os.path.abspath(options.profiles_dir), "%s.meme" % options.matrix_id), dummy_fasta, options.p_value_thresh):
                # If match in relative scores... #
                if (start + chunk_start, strand) in relative_scores:
                    functions.write(dummy_file, "%s\t%s\t%s\t%s" % (start + chunk_start, strand, relative_scores[(start + chunk_start, strand)], int(log(p_value) * 1000 / -10)))
            # Remove dummy FASTA file if exist #
            if os.path.exists(dummy_fasta): os.remove(dummy_fasta)
        # If dummy file exists... #
        if os.path.exists(dummy_file):
            # Compress #
            functions.compress(dummy_file, output_file)
            # Remove dummy file #
            os.remove(dummy_file)
            