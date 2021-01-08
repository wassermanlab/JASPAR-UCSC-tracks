#!/usr/bin/env python

import argparse
import os

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--devel", action="store_true", help="development mode (uses hfaistos; default = False)")
    parser.add_argument("-o", default="./", help="output directory (default = ./)", metavar="DIR")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Get profiles
    get_profiles(args.devel, args.o)

def get_profiles(devel=False, output_dir="./"):
    """
    For each taxon (except for urochordates), this function gets all profiles (in jaspar format) from the JASPAR CORE.
    """

    # Initialize
    version = 2020
    cwd = os.getcwd()
    taxons = ["fungi", "insects", "nematodes", "plants", "vertebrates"]

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
            url = "http://jaspar.genereg.net/download/CORE/JASPAR%s_CORE_%s_redundant_pfms_jaspar.zip" % (version, taxon)
            if devel:
                url = "http://hfaistos.uio.no:8002/download/CORE/JASPAR%s_CORE_%s_redundant_pfms_jaspar.zip" % (version, taxon)

            os.system("curl --silent -O %s" % url)

            # Unzip
            file_name = "JASPAR%s_CORE_%s_redundant_pfms_jaspar.zip" % (version, taxon)
            if devel:
                file_name = "JASPAR%s_CORE_%s_redundant_pfms_jaspar.zip" % (version, taxon)
            os.system("unzip -qq %s" % file_name)

            # Remove zip files
            os.remove("%s" % file_name)

            # Return to original directory
            os.chdir(cwd)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()
