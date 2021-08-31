#!/usr/bin/env python

import argparse
import os
import subprocess as sp

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

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Get profiles
    get_profiles(args.o)

    # To PWM
    jaspar_to_pwm(args.o)

def get_profiles(output_dir="./"):
    """
    For each taxon, this function downloads all profiles from the JASPAR CORE
    in JASPAR format.
    """

    # Initialize
    version = 2022
    cwd = os.getcwd()
    taxons = ["fungi", "insects", "nematodes", "plants", "urochordates",
        "vertebrates"]

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
    perl_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        "jasparconvert.pl")
    taxons = ["fungi", "insects", "nematodes", "plants", "urochordates",
        "vertebrates"]

    # For each taxon...
    for taxon in taxons:

        # Initialize
        taxon_dir = os.path.join(os.path.abspath(output_dir), taxon)

        # For each profile...
        for f in os.listdir(taxon_dir):

            # JASPAR to PWMScan
            matrix_file = os.path.join(taxon_dir, f)
            pwm_file = os.path.join(taxon_dir, f"{f[:8]}.pwm")
            options = f"-noheader -o {pwm_file} -w l"
            cmd = f"{perl_script} {options} {matrix_file}"
            _ = sp.run([cmd], stdout=sp.DEVNULL, stderr=sp.DEVNULL, shell=True)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()
