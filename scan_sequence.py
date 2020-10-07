#!/usr/bin/env python

import argparse
from Bio import motifs
from functools import partial
from itertools import chain
from multiprocessing import Pool
from numpy import log10 as log
import os
import re
import shutil
import subprocess
from tqdm import tqdm

# Authorship
__author__ = "Oriol Fornes"
__organization__ = "The JASPAR Consortium"
__version__ = "2.0.1"
__maintainer__ = "Oriol Fornes"
__email__ = "oriol@cmmt.ubc.ca"
__status__ = "Production"

# Globals
pid = os.getpid()
taxons = ["fungi", "insects", "nematodes", "plants", "vertebrates"]

#-------------#
# Functions   #
#-------------#

usage_msg = """
usage: %s --fasta-file FILE --files-dir DIR
""" % os.path.basename(__file__)

help_msg = """%s
  --fasta-file FILE   one or more sequences in FASTA format
  --profiles-dir DIR  output directory from get_profiles.py

optional arguments:
  -h, --help          show this help message and exit
  --dummy-dir DIR     dummy directory (default = /tmp/)
  --output-dir DIR    output directory (default = ./)
  --threads INT       threads to use (default = 1)

search arguments:
  -A FLOAT            background freq for A (default = 0.25)
  -C FLOAT            background freq for C (default = 0.25)
  -G FLOAT            background freq for G (default = 0.25)
  -T FLOAT            background freq for T (default = 0.25)
  -l, --latest        use the latest version of each profile
  --profile [STR ...] profile ID(s) to use (default = all)
  --pthresh FLOAT     p-value threshold (default = 0.05)
  --rthresh FLOAT     relative score threshold (default = 0.8)
  --taxon [STR ...]   taxon(s) to use (default = all)
""" % usage_msg

def parse_args():
    """
    This function parses arguments provided via the command line and returns an {argparse} object.
    """

    # Initialize
    parser = argparse.ArgumentParser(add_help=False)

    # Mandatory args
    parser.add_argument("--fasta-file")
    parser.add_argument("--profiles-dir")

    # Optional args
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-h", "--help", action="store_true")
    optional_group.add_argument("--dummy-dir", default="/tmp/")
    optional_group.add_argument("--output-dir", default="./")
    optional_group.add_argument("--threads", default=1)
    optional_group.add_argument("-A", default=0.25)
    optional_group.add_argument("-C", default=0.25)
    optional_group.add_argument("-G", default=0.25)
    optional_group.add_argument("-T", default=0.25)
    optional_group.add_argument("-l", "--latest", action="store_true")
    optional_group.add_argument("--profile", nargs="*", default=[])
    optional_group.add_argument("--pthresh", default=0.05)
    optional_group.add_argument("--rthresh", default=0.8)
    optional_group.add_argument("--taxon", nargs="*", default=taxons)

    args = parser.parse_args()

    check_args(args)

    return(args)

def check_args(args):
    """
    This function checks an {argparse} object.
    """

    # Print help
    if args.help:
        print(help_msg)
        exit(0)

    # Check mandatory arguments
    if not args.fasta_file or not args.profiles_dir:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "arguments \"--fasta-file\" \"--profiles-dir\" are required\n"]
        print(": ".join(error))
        exit(0)

    # Check "--threads" argument
    try:
        args.threads = int(args.threads)
    except:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--threads\"", "invalid value", "\"%s\"\n" % args.threads]
        print(": ".join(error))
        exit(0)

    # Check "-A", "-C", "-G", "-T"
    try:
        args.A = float(args.A)
        args.C = float(args.C)
        args.G = float(args.G)
        args.T = float(args.T)
    except:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--threads\"", "invalid value", "\"%s\"\n" % args.threads]
        print(": ".join(error))
        exit(0)

    # Check "--profile" argument
    for profile in args.profile:
        if not re.search("^MA\d{4}.\d$", profile):
            error = [os.path.basename(__file__), "error", "argument \"--profile\"", "invalid value", "\"%s\"" % profile]
            print(": ".join(error))
            exit(0)

    # Check "--pthresh" argument
    try:
        args.pthresh = float(args.pthresh)
    except:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--pthresh\"", "invalid value", "\"%s\"\n" % args.pthresh]
        print(": ".join(error))
        exit(0)

    # Check "--rthresh" argument
    try:
        args.rthresh = float(args.rthresh)
    except:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--rthresh\"", "invalid value", "\"%s\"\n" % args.rthresh]
        print(": ".join(error))
        exit(0)

    # Check "--taxon" argument
    for taxon in args.taxon:
        if taxon not in taxons:
            error = [os.path.basename(__file__), "error", "argument \"--taxon\"", "invalid value", "\"%s\"" % taxon]
            print(": ".join(error))
            exit(0)

    return(args) 

def main():

    # Parse arguments
    args = parse_args()

    # Scan sequence
    scan_sequence(args.fasta_file, args.profiles_dir, args.dummy_dir,
        args.output_dir, args.threads, args.A, args.C, args.G, args.T,
        args.latest, args.profile, args.pthresh, args.rthresh, args.taxon)

def scan_sequence(fasta_file, profiles_dir, dummy_dir="/tmp/", output_dir="./",
    threads=1, A=0.25, C=0.25, G=0.25, T=0.25, latest=False, profile=[], 
    pthresh=0.05, rthresh=0.8, taxon=taxons):

    # Create output directory
    if not os.path.exists(os.path.abspath(output_dir)):
        os.makedirs(os.path.abspath(output_dir))

    # Get profiles to scan
    profiles = _get_profiles(profiles_dir, latest, profile, taxon)

    # Scan profiles
    _scan_profiles(profiles, fasta_file, dummy_dir, output_dir, threads, A, C,
        G, T, pthresh, rthresh)

def _get_profiles(profiles_dir, latest=False, profile=[], taxon=taxons):

    # Initialize
    profiles = []
    profiles_dict = {}

    # For each taxon...
    for t in taxon:

        # Initialize
        taxon_dir = os.path.join(os.path.abspath(profiles_dir), t)

        # For each profile...
        for profile_file in sorted(os.listdir(taxon_dir), reverse=True):

            # Skip wrong profiles
            if len(profile) > 0:
                if profile_file[:8] not in profile:
                    continue

            # Load profile
            with open(os.path.join(taxon_dir, profile_file)) as f:
                p = motifs.read(f, "jaspar")

            # Initialize key
            key = profile_file[:6]
            profiles_dict.setdefault(key, [])

            # Skip profile if only using the latest version of each profile
            if latest:
                if len(profiles_dict[key]) == 1:
                    continue

            # Add profile
            profiles_dict[key].append(p)

    # Create list of profiles
    for value_list in profiles_dict.values():
        for p in value_list:
            profiles.append(p)

    return(profiles)

def _scan_profiles(profiles, fasta_file, dummy_dir="/tmp/", output_dir="./",
    threads=1, A=0.25, C=0.25, G=0.25, T=0.25, pthresh=0.05, rthresh=0.8):

    # Parallelize scan
    pool = Pool(threads)
    parallelized = partial(_scan_profile, fasta_file=fasta_file,
        dummy_dir=dummy_dir, output_dir=output_dir, threads=threads, A=A, C=C,
        G=G, T=T, pthresh=pthresh, rthresh=rthresh)
    for _ in tqdm(pool.imap(parallelized, profiles), total=len(profiles)):
        pass
    pool.close()
    pool.join()

def _scan_profile(profile, fasta_file, dummy_dir="/tmp/", output_dir="./",
    threads=1, A=0.25, C=0.25, G=0.25, T=0.25, pthresh=0.05, rthresh=0.8):

    # Initialize
    bin_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bin")
    cutoff = None
    dummy_file = os.path.join(dummy_dir, "%s.%s.%s" % (os.path.basename(__file__), pid, profile.matrix_id))
    pwm_file = "%s.pwm" % dummy_file
    tsv_file = "%s.tsv" % dummy_file
    gzipped_file = "%s.gz" % tsv_file
    output_file = os.path.join(output_dir, "%s.tsv.gz" % profile.matrix_id)

    # Skip if profile already scanned (i.e. for speed)
    output_file = os.path.join(output_dir, "%s.tsv.gz" % profile.matrix_id)
    if not os.path.exists(output_file):

        # Add background
        profile.background = {"A": A, "C": C, "G": G, "T": T}

        # Add JASPAR pseudocounts
        profile.pseudocounts = motifs.jaspar.calculate_pseudocounts(profile)

        # Write profile in PWMScan format
        with open(pwm_file, "w") as f:

            for i in range(len(profile.pssm["A"])):

                f.write("%s\n" % "\t".join([str(int(profile.pssm[j][i]*100)) for j in "ACGT"]))

        # Calculate distribution of PWM scores
        cmd = "%s %s" % (os.path.join(bin_dir, "matrix_prob"), pwm_file)
        process = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        with open(tsv_file, "w") as f:

            for line in process.stdout.decode("utf-8").split("\n"):

                matches = re.findall("(\S+)", line)

                if len(matches) == 3:

                    score = matches[0]
                    p_value = float(matches[1])
                    perc = float(matches[2][:-1])

                    f.write("%s\t%s\t%s\n" % (score, int(perc * 10), int(log(p_value) * 1000 / -10)))

                    # Get PWM score cutoff
                    if p_value < pthresh and perc >= rthresh * 100 or cutoff is None:
                        cutoff = score

        # Scan FASTA file
        cmd_1 = "%s -m %s -c %s %s" % (os.path.join(bin_dir, "matrix_scan"), pwm_file, cutoff, fasta_file)
        cmd_2 = "gzip > %s" % gzipped_file
        cmd = '''%s | awk -v score_tab="%s" -v name="%s" 'BEGIN { while((getline line < score_tab) > 0 ) {split(line,f," "); scores[f[1]]=f[2]; pvalues[f[1]]=f[3]} close(score_tab) } {print $1"\t"$2"\t"$3"\t"name"\t"scores[$5]"\t"pvalues[$5]"\t"$6}' | %s''' % (cmd_1, tsv_file, profile.name, cmd_2)
        subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT)

        # Write output
        shutil.copy(gzipped_file, output_file)

        # Remove dummy files
        os.remove(pwm_file)
        os.remove(tsv_file)
        os.remove(gzipped_file)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()
