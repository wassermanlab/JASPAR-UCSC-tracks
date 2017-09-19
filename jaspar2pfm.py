#!/usr/bin/env python2.7
import os, re
import optparse
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

    #parser = optparse.OptionParser("./%prog [-b <jaspar_bundle> -o <output_dir> -t <tax_groups>]")
    parser = optparse.OptionParser("./%prog -b <jaspar_bundle> [-o <output_dir> -t <tax_groups>]")

    parser.add_option("-b", action="store", type="string", dest="jaspar_bundle", help="JASPAR bundle (file containing many JASPAR profiles)", metavar="<jaspar_bundle>")
    parser.add_option("-o", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="<output_dir>")
    #parser.add_option("-t", default="vertebrates", action="store", type="string", dest="tax_groups", help="Taxonomic groups (comma separated; e.g. \"vertebrates,nematodes,insects\"; default = vertebrates)", metavar="<tax_groups>")

    (options, args) = parser.parse_args()

    if options.jaspar_bundle is None:
        parser.error("missing arguments: type option \"-h\" for help")

#    tax_groups = set(["fungi", "insects", "nematodes", "plants", "urochordates", "vertebrates"])
#    for tax_group in options.tax_groups.split(","):
#        if tax_group not in tax_groups:
#            parser.error("invalid taxonomic group: %s\n\tallowed taxonomic groups include: \"fungi\", \"insects\", \"nematodes\", \"plants\", \"urochordates\" and \"vertebrates\"\n\tmultiple taxonomic groups can be provided as comma-separated values: e.g. \"vertebrata,nematoda\"" % tax_group)

    return options

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Create output directory #
    if not os.path.exists(os.path.abspath(options.output_dir)):
        os.makedirs(os.path.abspath(options.output_dir))

#    # If JASPAR bundle... #
#    if options.jaspar_bundle is not None:
    # Initialize #
    profiles = []
    # For each line... #
    for line in functions.parse_file(os.path.join(options.jaspar_bundle)):
        # If header... #
        if line.startswith(">"):
            # Add profile #
            profiles.append([line])
        # ... Else... #
        else:
            m = re.search("^\s*\w\s*\[(.+)\]\s*$", line)
            # Add counts to profile #
            if m: profiles[-1].append(" ".join([i for i in re.findall("[+-]?[0-9]*[.]?[0-9]+", m.group(1))]))
    # For each profile... #
    for profile in profiles:
        # Initialize #
        matrix_id = profile[0][1:9]
        # Skip if JASPAR file already exists #
        jaspar_file = os.path.join(os.path.abspath(options.output_dir), "%s.pfm" % matrix_id)
        if not os.path.exists(jaspar_file):
            # Write #
            functions.write(jaspar_file, "\n".join(profile))
#    # ... Else... #
#    else:
#        # For each profile... #
#        for profile in functions.get_jaspar_profiles():
#            # Skip if invalid taxa #
#            if profile.tax_group not in options.tax_groups: continue
#            # Skip if JASPAR file already exists #
#            jaspar_file = os.path.join(os.path.abspath(options.output_dir), "%s.pfm" % profile.matrix_id)
#            if not os.path.exists(jaspar_file):
#                # Write #
#                functions.write(jaspar_file, ">%s %s" % (profile.matrix_id, profile.name))
#                functions.write(jaspar_file, profile.format("pfm"))
