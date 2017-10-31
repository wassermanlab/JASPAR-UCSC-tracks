#!/usr/bin/env python2.7
import os, sys, re
from Bio import motifs
from numpy import log10 as log
import optparse
import shutil

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

    parser = optparse.OptionParser("./%prog -i <input_dir> -p <profiles_dir> [-c <chr> --dummy=<dummy_dir> -f <format> -m <matrix_id> -o <output_file> --pv-thresh=<p_value_thresh> --rs-thresh=<rel_score_thresh> -s <scores> -z]")

    parser.add_option("-c", action="store", type="string", dest="chr", help="Chromosome (reports TFBSs for the given chromosome; e.g 1-22, X, Y, M...; provide multiple chromosomes using commas e.g. \"chr2\"; default = None)", metavar="<chr>")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="<dummy_dir>")
    parser.add_option("-f", default="bed", action="store", type="string", dest="format", help="Format for output results (i.e. \"bed\", \"csv\" or \"tsv\"; default = bed)", metavar="<format>")
    parser.add_option("-i", action="store", type="string", dest="input_dir", help="Input directory with TFBSs matches (i.e. output directory from jaspar_search.py)", metavar="<input_dir>")
    parser.add_option("-m", action="store", type="string", dest="matrix_id", help="Matrix ID (reports TFBSs for the given matrix id; e.g. \"MA0002.2\"; provide multiple matrix ids using commas; default = None)", metavar="<matrix_id>")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="<output_file>")
    parser.add_option("-p", action="store", type="string", dest="profiles_dir", help="Profiles directory (from jaspar2pfm.py)", metavar="<profiles_dir>")
    parser.add_option("--pv-thresh", default=0.05, action="store", type="float", dest="p_value_thresh", help="P-value threshold (reports TFBSs over the given threshold; default = 0.05)", metavar="<p_value_thresh>")
    parser.add_option("--rs-thresh", default=0.8, action="store", type="float", dest="rel_score_thresh", help="Relative score threshold (reports TFBSs below the given threshold; default = 0.8)", metavar="<rel_score_thresh>")
    parser.add_option("-s", default="p_value", action="store", type="string", dest="scores", help="TFBS scores (reports TFBS scores as \"p_value\", \"rel_score\" or \"both\", this last option available only for \"csv\" and \"tsv\" formats; default = p_value)", metavar="<scores>")
    parser.add_option("-z", "--gzip", default=False, action="store_true", dest="compress", help="Compress output file using Gzip (default = False)")

    (options, args) = parser.parse_args()

    if options.input_dir is None or options.profiles_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

#    if options.chr is not None:
#        for chromosome in options.chr.split(","):
#            if not re.search("^chr[0-9XYM]{1,2}$", chromosome):
#                parser.error("invalid chromosome: %s\n\tvalid chromosomes include 1-22, X, Y and M: e.g. \"chr2\"" % chromosome)

    if options.scores != "p_value" and options.scores != "rel_score" and options.scores != "both":
        parser.error("invalid type of TBFS score:%s\n\tTFBS scores can only be reported as \"p_value\", \"rel_score\" or \"both\", this last option available only for \"csv\" and \"tsv\" formats" % options.scores)

    return options

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Initialize #
    profiles = {}
    delimiter = "\t"
    if options.format == "csv": delimiter = ","
    if options.output_file is None: dummy_file = None
    else: dummy_file = os.path.join(os.path.abspath(options.dummy_dir), "%s.txt" % os.getpid())
    rel_score_thresh = int(options.rel_score_thresh * 1000) # transform relative score threshold
    p_value_thresh = int(log(options.p_value_thresh) * 1000 / -10) # transform p-value threshold
    # Remove dummy file if exist #
    if dummy_file is not None:
        if os.path.exists(dummy_file): os.remove(dummy_file)
    # Write #
    if options.format != "bed":
        header = delimiter.join(["chr", "start (1-based)", "end"])
        if options.scores == "rel_score": header += delimiter + "rel_score * 1000"
        elif options.scores == "p_value": header += delimiter + "p_value"
        else: header += delimiter + "rel_score * 1000" + delimiter + "-1 * log10(p_value) * 100"
        # Write #
        functions.write(dummy_file, header + delimiter + "strand")
    # For each matrix id and for each chr file... #
    for file_name in os.listdir(os.path.abspath(options.input_dir)):
        # Initialize #
        m = re.search("^(MA\d+\.\d)\.(chr[0-9XYM]{1,2})\.tab\.gz$", file_name)
        if not m: continue
        matrix_id = m.group(1)
        chromosome = m.group(2)
        # Skip file if wrong matrix id #
        if options.matrix_id is not None:
            if matrix_id not in options.matrix_id.split(","): continue
        # Skip file if wrong chromosome #
        if options.chr is not None:
            if chromosome not in options.chr.split(","): continue
        # If no profile for matrix id... #
        if matrix_id not in profiles:
            # Initialize #
            if os.path.exists(os.path.join(os.path.abspath(options.profiles_dir), "%s.pfm" % matrix_id)):
                pfm_file = os.path.join(os.path.abspath(options.profiles_dir), "%s.pfm" % matrix_id)
            elif os.path.exists(os.path.join(os.path.abspath(options.profiles_dir), "%s.jaspar" % matrix_id)):
                pfm_file = os.path.join(os.path.abspath(options.profiles_dir), "%s.jaspar" % matrix_id)
            else: raise ValueError("PFM in JASPAR format does not exist!")
            # Load profile #
            with open(pfm_file) as f:
                # Add profile to profiles #
                profiles.setdefault(matrix_id, motifs.read(f, "jaspar"))
        # For each line...
        for line in functions.parse_tsv_file(os.path.join(os.path.abspath(options.input_dir), file_name), gz=True):
            # Initialize #
            position = int(line[0])
            start = position
            if options.format == "bed": start -= 1 # convert to BED format (0-based)
            end = position + profiles[matrix_id].length - 1
            strand = line[1]
            rel_score = int(line[2])
            p_value =  int(line[3])
            # Skip matches that do not pass any of the score thresholds #
            if rel_score < rel_score_thresh or p_value < p_value_thresh: continue
            # For BED files, cap scores at 1000 (UCSC Genome Browser does not allow scores >1000) #
            if options.scores == "rel_score":
                if options.format == "bed": score = min([int(line[2]), 1000])
                else: score = line[2]
            elif options.scores == "p_value":
                if options.format == "bed": score = min([int(line[3]), 1000])
                else: score = line[3]
            # If both relative scores and p-values are required... #
            else: score = delimiter.join([line[2], line[3]])
            # Write #
            functions.write(dummy_file, delimiter.join(map(str, [chromosome, start, end, profiles[matrix_id].name, score, strand])))
    # If dummy file exists... #
    if dummy_file is not None:
        if os.path.exists(dummy_file):
            # Initialize #
            output_file = os.path.abspath(options.output_file)
            # If compress... #
            if options.compress:
                # Initialize #
                if not output_file.endswith(".gz"): output_file += ".gz"
                # Compress #
                functions.compress(dummy_file, output_file)
            # ... Else... #
            else:
                # Copy #
                shutil.copy(dummy_file, output_file)
            # Remove dummy file #
            os.remove(dummy_file)
