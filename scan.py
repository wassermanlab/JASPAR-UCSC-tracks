#!/usr/bin/env python2.7
import os, re
from Bio import motifs
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import optparse

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

    parser = optparse.OptionParser("./%prog -f <fasta_file> -m <matrix_file> [-o <output_file> -t <threshold>]")

    parser.add_option("-f", action="store", type="string", dest="fasta_file", help="FASTA file (e.g. chr1.fa)", metavar="<fasta_file>")
    parser.add_option("-m", action="store", type="string", dest="matrix_file", help="JASPAR profile (e.g. MA0002.2.pfm)", metavar="<matrix_file>")
    parser.add_option("-o", default="./", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="<output_file>")
    parser.add_option("-t", default="80%", action="store", type="string", dest="threshold", help="Min. score for TFBS matches; specify as absolute, e.g. 14.1, or relative score, e.g. 80%; default = 80%", metavar="<threshold>")

    (options, args) = parser.parse_args()

    if options.fasta_file is None or options.matrix_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Load profile #
    with open(os.path.abspath(options.matrix_file)) as f:
        profile = motifs.read(f, "jaspar")
    # Get pseudocounts #
    profile.pseudocounts = motifs.jaspar.calculate_pseudocounts(profile)
    # Get score threshold #
    m = re.search("^(.+)\%$", options.threshold)
    if m: score_threshold = (profile.pssm.max - profile.pssm.min) * float(m.group(1))/100 + profile.pssm.min
    else: score_threshold = float(options.threshold)
    # For each header, sequence... #
    for header, sequence in functions.parse_fasta_file(os.path.abspath(options.fasta_file)):
        # Format sequence #
        sequence = Seq(sequence, IUPAC.unambiguous_dna)
        # For each match position, score... #
        for position, score in profile.pssm.search(sequence, threshold=score_threshold):
            if position < 0:
                position += len(profile)
                strand = "-"
            else:
                strand = "+"
            # Write #
            functions.write(options.output_file, "%s\t%s\t%s\t%s\t%.3f" % (header, position, position + len(profile), strand, (score - profile.pssm.min) / (profile.pssm.max - profile.pssm.min)))
            