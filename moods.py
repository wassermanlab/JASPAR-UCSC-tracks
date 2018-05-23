#!/usr/bin/env python2.7
import os, re
from Bio import motifs
import MOODS.scan
import MOODS.tools
import MOODS.parsers
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

    parser = optparse.OptionParser("./%prog -f <fasta_file> -m <matrix_file> [--dummy=<dummy_dir> -o <output_file> -p <pvalue_threshold> -s <score_threshold>]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="<dummy_dir>")
    parser.add_option("-f", action="store", type="string", dest="fasta_file", help="FASTA file (e.g. chr1.fa)", metavar="<fasta_file>")
    parser.add_option("-m", action="store", type="string", dest="matrix_file", help="JASPAR profile (e.g. MA0002.2.pfm)", metavar="<matrix_file>")
    parser.add_option("-o", default="./", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="<output_file>")
    parser.add_option("-p", default=0.05, action="store", type=float, dest="pvalue_threshold", help="P-value threshold for TFBS matches; default = 0.05", metavar="<pvalue_threshold>")
    parser.add_option("-s", default="80%", action="store", type="string", dest="score_threshold", help="Score threshold for TFBS matches; specify as absolute, e.g. 14.1, or relative score, e.g. 80%; default = 80%", metavar="<score_threshold>")

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
        profile.pseudocounts = motifs.jaspar.calculate_pseudocounts(profile)
    # Get score threshold #
    m = re.search("^(.+)\%$", options.score_threshold)
    if m: score_threshold = (profile.pssm.max - profile.pssm.min) * float(m.group(1))/100 + profile.pssm.min
    else: score_threshold = float(options.threshold)
    # Create dummy profile #
    dummy_file = os.path.join(os.path.abspath(options.dummy_dir), "%s.%s.pfm" % (os.path.basename(__file__), os.getpid()))
    functions.write(dummy_file, motifs.write([profile], "pfm"))
    # Initialize #
    bg = MOODS.tools.flat_bg(4)
    pseudocounts = 0.01
    profile = MOODS.parsers.pfm_to_log_odds(dummy_file, bg, pseudocounts)
    pvalue_threshold = MOODS.tools.threshold_from_p(profile, bg, options.pvalue_threshold)
    # For each header, sequence... #
    for header, sequence in functions.parse_fasta_file(os.path.abspath(options.fasta_file)):
        # Get results #
        print("here")
        exit(0)
        # results = MOODS.scan.scan_dna(sequence, profile, bg, pvalue_threshold, 7)
        # For each result... #
        for result in results:
            print(dir(result))
            exit(0)
            # Write #
            functions.write(options.output_file, "%s\t%s\t%s\t%s\t%.3f" % (header, position + 1, position + len(profile), strand, (score - profile.pssm.min) / (profile.pssm.max - profile.pssm.min)))
    # Remove dummy profile #
    os.remove(dummy_file)

