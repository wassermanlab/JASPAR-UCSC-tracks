#!/usr/bin/env python

import argparse
from Bio import motifs
import json
import os
import subprocess as sp

# Globals
taxons = [
    "fungi",
    "insects",
    "nematodes",
    "plants",
    "urochordates",
    "vertebrates"
]

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns
    an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-o",
        default="./",
        help="output directory (default = ./)",
        metavar="DIR"
    )
    parser.add_argument(
        "-v",
        default=2024,
        help="JASPAR version (default = 2024)",
        metavar="INT"
    )

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Get profiles
    get_profiles(args.o, args.v)

    # Convert profiles to PWMs
    jaspar_to_pwm(args.o)

    # Get profile names
    get_names(args.o)

def get_profiles(output_dir="./", version=2024):
    """
    For each taxon, this function downloads all profiles from the JASPAR CORE
    in JASPAR format.
    """

    # Initialize
    cwd = os.getcwd()

    # Create output directory
    if not os.path.exists(os.path.abspath(output_dir)):
        os.makedirs(os.path.abspath(output_dir))

    # For each taxon...
    for taxon in taxons:

        # Skip if taxon directory already exists
        taxon_dir = os.path.join(os.path.abspath(output_dir), taxon)
        if not os.path.exists(taxon_dir):

            # Create taxon directory
            os.makedirs(taxon_dir)

            # Move to taxon directory
            os.chdir(taxon_dir)

            # Get JASPAR profiles
            file_name = "JASPAR%s_CORE_%s_redundant_pfms_jaspar.zip" % \
                (version, taxon)
            url = "https://testjaspar.uio.no/download/data/%s/CORE/%s" % \
                (version, file_name)
            os.system("curl --silent -O %s" % url)

            # Unzip
            os.system("unzip -qq %s" % file_name)

            # Remove zip files
            os.remove("%s" % file_name)

            # Return to original directory
            os.chdir(cwd)

def jaspar_to_pwm(output_dir="./"):
    """
    For each taxon, this function reformats all profiles from JASPAR to
    PWMScan format.
    """

    # Initialize
    # perl_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),
    #     "jasparconvert.pl")

    # For each taxon...
    for taxon in taxons:

        # Initialize
        taxon_dir = os.path.join(os.path.abspath(output_dir), taxon)

        # For each profile...
        for f in os.listdir(taxon_dir):

            # Skip non-JASPAR profiles
            if not f.endswith(".jaspar"):
                continue

            # JASPAR to PWMScan
            with open(os.path.join(taxon_dir, f)) as handle:
                m = motifs.read(handle, "jaspar")
            m.pseudocounts = motifs.jaspar.calculate_pseudocounts(m)
            pwm = list(map(list, zip(*[m.pssm[nt] for nt in "ACGT"])))
            pwm_file = os.path.join(taxon_dir, f"{f[:8]}.pwm")
            if not os.path.exists(pwm_file):
                with open(pwm_file, "w") as handle:
                    for i in pwm:
                        s = " ".join(["{:7d}".format(round(j*100)) for j in i])
                        handle.write("%s\n" % s)

def get_names(output_dir="./"):
    """
    This function extracts the name of each JASPAR profile and saves them
    in a JSON file.
    """

    # Initialize
    names = {}

    # Skip if already done
    json_file = os.path.join(output_dir, "names.json")
    if not os.path.exists(json_file):

        # For each taxon...
        for taxon in taxons:

            # Initialize
            taxon_dir = os.path.join(os.path.abspath(output_dir), taxon)

            # For each profile...
            for f in os.listdir(taxon_dir):

                # Skip non-JASPAR profiles
                if not f.endswith(".jaspar"):
                    continue

                # Get profile name
                with open(os.path.join(taxon_dir, f)) as handle:
                    m = motifs.read(handle, "jaspar")
                if m.name.startswith(m.matrix_id):
                    name = m.name[len(m.matrix_id)+1:]
                else:
                    name = m.name
                names.setdefault(m.matrix_id, name)

        # Write JSON
        with open(json_file, "w") as handle:
            json.dump(names, handle, sort_keys=True, indent=4)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()
