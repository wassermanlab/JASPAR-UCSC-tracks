#!/usr/bin/env python2.7
import os
import optparse

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("./%prog [-o <output_dir>]")

    parser.add_option("-o", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="<output_dir>")

    (options, args) = parser.parse_args()

    return options

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
    cwd = os.getcwd()
    taxons = ["fungi", "insects", "nematodes", "plants", "vertebrates"]

    # For each taxon... #
    for taxon in taxons:
        # Skip if taxon directory already exists #
        taxon_dir = os.path.join(os.path.abspath(options.output_dir), taxon)
        if not os.path.exists(taxon_dir):
            # Create taxon dir #
            os.makedirs(taxon_dir)
            # Change directory #
            os.chdir(taxon_dir)
            # Download JASPAR profiles #
            url = "http://jaspar.genereg.net/download/CORE/JASPAR2018_CORE_%s_redundant_pfms_jaspar.zip" % taxon
            os.system("curl --silent -O %s" % url)
            # Unzip #
            file_name = "JASPAR2018_CORE_%s_redundant_pfms_jaspar.zip" % taxon
            os.system("unzip -qq %s" % file_name)
            # Remove SQL files #
            os.remove("%s" % file_name)
            # Return to original directory #
            os.chdir(cwd)
