#!/usr/bin/env python2.
import os, re
from Bio import motifs
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

    #parser = optparse.OptionParser("./%prog -m <meme_dir> [-b <jaspar_bundle> --dummy=<dummy_dir> -o <output_dir> -t <tax_groups>]")
    parser = optparse.OptionParser("./%prog -b <jaspar_bundle> -m <meme_dir> [--dummy=<dummy_dir> -o <output_dir>]")

    parser.add_option("-b", action="store", type="string", dest="jaspar_bundle", help="JASPAR bundle (file containing many JASPAR profiles)", metavar="<jaspar_bundle>")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="<dummy_dir>")
    parser.add_option("-m", action="store", type="string", dest="meme_dir", help="Full path to MEME bin directory (i.e. where all MEME executables are located; e.g. $MEME_PATH/bin)", metavar="<meme_dir>")
    parser.add_option("-o", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="<output_dir>")
    #parser.add_option("-t", default="vertebrates", action="store", type="string", dest="tax_groups", help="Taxonomic groups (comma separated; e.g. \"vertebrates,nematodes,insects\"; default = vertebrates)", metavar="<tax_groups>")

    (options, args) = parser.parse_args()

    if options.jaspar_bundle is None or options.meme_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    if not os.path.exists(os.path.join(os.path.abspath(options.meme_dir), "fimo")):
        parser.error("incorrect MEME bin directory: \"fimo\" was not found!")
        
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
    dummy_file = os.path.join(os.path.abspath(options.dummy_dir), "bundle.txt")
    profiles = {}
    header = []
    profile = None
    # If dummy file exists #
    if os.path.exists(dummy_file): os.remove(dummy_file)
    # For each line... #
    for line in functions.parse_file(os.path.join(options.jaspar_bundle)):
        # If header... #
        if line.startswith(">"):
            # Write #
            functions.write(dummy_file, line)
        # ... Else... #
        else:
            m = re.search("^\s*(\w)\s*\[(.+)\]\s*$", line)
            # Write #
            if m: functions.write(dummy_file, "%s" % " ".join(map(str, [int(float(i)) for i in re.findall("[+-]?[0-9]*[.]?[0-9]+", m.group(2))])))
    # Reformat JASPAR profiles to MEME profiles #
    process =  subprocess.check_output([os.path.join(os.path.abspath(options.meme_dir), "jaspar2meme"), "-bundle", dummy_file], stderr=subprocess.STDOUT)
    # For each line... #
    for line in process.split("\n"):
        if line.startswith("MOTIF"):
            # Initialize #
            profile = re.search("^MOTIF (MA\d{4}.\d{1}) \S+$", line)
            profiles.setdefault(profile.group(1), [])
            profiles[profile.group(1)] += header
        try: profiles[profile.group(1)].append(line)
        except: header.append(line) # no profile
    # Write JASPAR profiles in MEME format #
    for profile in sorted(profiles):
        # Skip if MEME file already exists #
        meme_file = os.path.join(os.path.abspath(options.output_dir), "%s.meme" % profile)
        if not os.path.exists(meme_file):
            # Write #
            functions.write(meme_file, "\n".join(profiles[profile]))
#    # ... Else... #
#    else:
#        # For each profile... #
#        for profile in functions.get_jaspar_profiles():
#            # Skip if invalid taxa #
#            if profile.tax_group not in options.tax_groups: continue
#            # Skip if JASPAR file already exists #
#            jaspar_file = os.path.join(os.path.abspath(options.output_dir), "%s.transfac" % profile.matrix_id)
#            if not os.path.exists(jaspar_file):
#                # Initialize #
#                functions.write(jaspar_file, "ID %s" % profile.matrix_id)
#                functions.write(jaspar_file, "BF %s" % profile.tax_group)
#                # Save JASPAR profile in TRANSFAC format #
#                functions.write(jaspar_file, motifs.write([profile], format="transfac"))
#            # Skip if MEME file already exists #
#            meme_file = os.path.join(os.path.abspath(options.output_dir), "%s.meme" % profile.matrix_id)
#            if not os.path.exists(meme_file):
#                # Reformat JASPAR profile to MEME profile #
#                process =  subprocess.check_output([os.path.join(os.path.abspath(options.meme_dir), "transfac2meme"), jaspar_file], stderr=subprocess.STDOUT)
#                # For each line... #
#                for line in process.split("\n"):
#                    # Write #
#                    functions.write(meme_file, line)
