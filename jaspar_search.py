#!/usr/bin/env python2.7

import os, re
from Bio import motifs
from functools import partial
from itertools import chain
from multiprocessing import Pool
from numpy import log10 as log
import optparse
import shutil
import subprocess
from tqdm import tqdm

# Import my functions #
import functions

# Authorship #
__author__ = "Oriol Fornes"
__credits__ = ["Aziz Khan", "Oriol Fornes", "Marius Gheorghe", "Robin van der Lee",
    "David J. Arenillas", "Wyeth W. Wasserman", "Anthony Mathelier"]
__organization__ = "The JASPAR Consortium"
__copyright__ = "Copyright 2018"
__license__ = "LGPL"
__version__ = "2.0.1"
__maintainer__ = "Oriol Fornes"
__email__ = "oriol@cmmt.ubc.ca"
__status__ = "Development"

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("./%prog -f <files_dir> -i <input_file> -s <pwmscan_dir> [-b <background> --dummy=<dummy_dir> -m <matrix_id> -o <output_dir> --pv-thresh=<p_value_thresh> --rs-thresh=<rel_score_thresh> -t <taxon>] [-l]")

    parser.add_option("-f", action="store", type="string", dest="files_dir", help="Files directory (output directory from make_files.py)", metavar="<files_dir>")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (i.e. one or more sequences in FASTA format)", metavar="<input_file>")
    parser.add_option("-s", action="store", type="string", dest="pwmscan_dir", help="Full path to PWMScan bin directory (i.e. where PWMScan executables are located; e.g. $PWMSCAN_PATH/bin)", metavar="<pwmscan_dir>")

    group = optparse.OptionGroup(parser, "Non-mandatory options")
    group.add_option("-b", "--background", default="0.25,0.25,0.25,0.25", action="store", dest="background", help="Background frequency for A, C, G, and T (e.g. for hg38 use \"0.29,0.21,0.21,0.29\"; default = \"0.25,0.25,0.25,0.25\")", metavar="<background>")
    group.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="<dummy_dir>")
    group.add_option("-m", action="store", type="string", dest="matrix_id", help="JASPAR matrix ID (more than one can be provided as comma-separated values; e.g. \"MA0002.2,MA0003.3\")", metavar="<matrix_id>")
    group.add_option("-o", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="<output_dir>")
    group.add_option("--pv-thresh", default=0.05, action="store", type="float", dest="p_value_thresh", help="P-value threshold (default = 0.05)", metavar="<p_value_thresh>")
    group.add_option("--rs-thresh", default=0.8, action="store", type="float", dest="rel_score_thresh", help="Relative score threshold (default = 0.8)", metavar="<rel_score_thresh>")
    group.add_option("-t", action="store", dest="taxon", help="Taxonomic group (i.e. \"fungi\", \"insects\", \"nematodes\", \"plants\", or \"vertebrates\"; default = None)", metavar="<taxon>")
    group.add_option("--threads", default=1, action="store", type="int", dest="threads", help="Total number of cores to be used for the computation (default = 1)", metavar="<threads>")
    parser.add_option_group(group)

    group = optparse.OptionGroup(parser, "Search modes")
    group.add_option("-l", "--latest", default=False, action="store_true", dest="latest", help="Latest mode (only use the latest version of each profile; default = False)")
    parser.add_option_group(group)

    (options, args) = parser.parse_args()

    if options.input_file is None or options.pwmscan_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    for i in options.background.split(","):
        try:
            float(i)
        except ValueError:
            parser.error("invalid value for background frequency: \"%s\"\n\tPlease specify A, C, G and T background frequencies as comma-separated values" % i)
    if len(options.background.split(",")) != 4: parser.error("invalid background frequencies: \"%s\"\n\tPlease specify A, C, G and T background frequencies as comma-separated values" % options.background)

    if options.matrix_id is not None:
        for i in options.matrix_id.split(","):
            if not re.search("^MA\d{4}.\d$", i):
                parser.error("invalid matrix id: \"%s\"" % i)

    if options.taxon is not None:
        if options.taxon not in ["fungi", "insects", "nematodes", "plants", "vertebrates"]:
            parser.error("invalid taxon: %s\n\tvalid taxons include \"fungi\", \"insects\", \"nematodes\", \"plants\", and \"vertebrates\"" % options.taxon)

    return options

def scan_profiles(profiles, input_file, pwmscan_dir, output_dir="./", background=[0.25, 0.25, 0.25, 0.25],
    p_value_thresh=0.05, rel_score_thresh=0.8, dummy_dir="/tmp/", threads=1):
    """
    """

    # Parallelize PWM scan #
    pool = Pool(threads)
    parallelize_scan_profile = partial(pwm_scan, input_file, pwmscan_dir, output_dir, background,
    p_value_thresh, rel_score_thresh, dummy_dir)
    for _ in tqdm(pool.imap(parallelize_scan_profile, iter(profiles)), desc="PWM scan", total=len(profiles)):
        pass
    pool.close()
    pool.join()

def scan_profile(input_file, pwmscan_dir, output_dir="./", background=[0.25, 0.25, 0.25, 0.25],
    p_value_thresh=0.05, rel_score_thresh=0.8, dummy_dir="/tmp/", profile=None):
    """
    """

    # Initialize #
    bed_file = os.path.join(dummy_dir, "%s.%s.%s.bed" % (os.path.basename(__file__), os.getpid(), profile.matrix_id))
    pwm_file = os.path.join(dummy_dir, "%s.%s.%s.pwm" % (os.path.basename(__file__), os.getpid(), profile.matrix_id))
    tsv_file = os.path.join(dummy_dir, "%s.%s.%s.tsv" % (os.path.basename(__file__), os.getpid(), profile.matrix_id))
    output_file = os.path.join(output_dir, "%s.bed.gz" % profile.matrix_id)

    # Add background #
    profile.background = {"A": background[0], "C": background[1], "G": background[2], "T": background[3]}
    # Add JASPAR pseudocounts #
    profile.pseudocounts = motifs.jaspar.calculate_pseudocounts(profile)
    # Convert to PWMScan format #
    for i in range(len(profile.pssm["A"])):
        functions.write(pwm_file, "\t".join([str(int(profile.pssm[j][i] * 100)) for j in "ACGT"]))
    # Calculate distribution of matrix scores #
    try:
        process = subprocess.check_output([os.path.join(pwmscan_dir, "matrix_prob"), pwm_file], stderr=subprocess.STDOUT)
        for line in process.split("\n"):
            m = re.search("(\S+)\s+(\S+)\s+(\S+)%", line)
            if m:
                score = m.group(1)
                p_value = float(m.group(2))
                perc = float(m.group(3))
                functions.write(tsv_file, "%s\t%s\t%s" % (score, int(perc * 10), int(log(p_value) * 1000 / -10)))
                if p_value < options.p_value_thresh and perc >= options.rel_score_thresh * 100:
                    cutoff = score
    except: raise ValueError("Could not calculate distribution of matrix scores!")
    # Scan DNA sequence for TFBS matches #
    try:
        bash_command = '''%s -m %s -c %s %s | awk -v score_tab="%s" -v name="%s" 'BEGIN { while((getline line < score_tab) > 0 ) {split(line,f," "); scores[f[1]]=f[2]; pvalues[f[1]]=f[3]} close(score_tab) } {print $1"\t"$2"\t"$3"\t"name"\t"scores[$5]"\t"pvalues[$5]"\t"$6}' | gzip > %s''' % (os.path.join(pwmscan_dir, "matrix_scan"), pwm_file, cutoff, os.path.abspath(options.input_file), tsv_file, profile.name, bed_file)
        process = subprocess.call(bash_command, shell=True, stderr=subprocess.STDOUT)
    except:
        raise ValueError("Could not scan DNA sequence file for TFBS matches!")
    # Write output #
    shutil.copy(bed_file, output_file)
    # Remove files #
    os.remove(bed_file)
    os.remove(pwm_file)
    os.remove(tsv_file)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Create output dir #
    if not os.path.exists(os.path.abspath(options.output_dir)):
        os.makedirs(os.path.abspath(options.output_dir))

    # Initialize #
    profiles = {}
    taxons = ["fungi", "insects", "nematodes", "plants", "vertebrates"]

    # For each taxon... #
    for taxon in taxons:
        # Skip if wrong taxon #
        if options.taxon is not None:
            if taxon != options.taxon: continue
        # For each profile... #
        for profile_file in sorted(os.listdir(os.path.join(os.path.abspath(options.files_dir), taxon)), reverse=True):
            # Skip if wrong matrix ID #
            if options.matrix_id is not None:
                if profile_file[:8] not in options.matrix_id: continue
            # Load profile #
            with open(os.path.join(os.path.abspath(options.files_dir), taxon, profile_file)) as f:
                profile = motifs.read(f, "jaspar")
            # Add profile #
            profiles.setdefault(profile_file[:6], [])
            if options.latest:
                if len(profiles[profile_file[:6]]) == 1: continue
            profiles[profile_file[:6]].append(profile)

    # Create a list of profiles to scan #
    profiles = list(chain.from_iterable(profiles.values()))

    # Scan profiles #
    scan_profiles(profiles, os.path.abspath(options.input_file), os.path.abspath(options.pwmscan_dir),
        os.path.abspath(options.output_dir), map(float, options.background.split(",")), options.p_value_thresh,
        options.rel_score_thresh, os.path.abspath(options.dummy_dir), options.threads)
