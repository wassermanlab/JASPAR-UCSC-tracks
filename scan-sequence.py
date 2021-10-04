#!/usr/bin/env python

import click
from click_option_group import optgroup
from functools import partial
from itertools import chain
import json
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
__version__ = "2021.9.1"
__maintainer__ = "Oriol Fornes"
__email__ = "oriol@cmmt.ubc.ca"
__status__ = "Production"

# Globals
pid = os.getpid()
taxons = ["fungi", "insects", "nematodes", "plants", "urochordates",
    "vertebrates"]

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "fasta_file",
    type=click.Path(exists=True, resolve_path=True),
)
@click.argument(
    "profiles_dir",
    type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-d", "--dummy-dir",
    help="Dummy directory.",
    type=click.Path(exists=True, resolve_path=True),
    default="/tmp/",
    show_default=True,
)
@click.option(
    "-o", "--output-dir",
    help="Output directory.",
    type=click.Path(resolve_path=True),
    default="./",
    show_default=True,
)
@click.option(
    "-t", "--threads",
    help="Number of CPU threads to use.",
    type=int,
    default=1,
    show_default=True,
)
@optgroup.group("Search arguments")
@optgroup.option(
    "-b", "--background",
    help="A, C, G, T background probabilities.",
    type=float,
    nargs=4,
    default=[.25, .25, .25, .25],
    show_default=True,
)
@optgroup.option(
    "-l", "--latest",
    help="Use the latest version of each profile.",
    is_flag=True,
)
@optgroup.option(
    "--profile",
    help="Profile ID(s) to use.  [default: all]",
    multiple=True,
)
@optgroup.option(
    "--pthresh",
    help="P-value threshold.",
    type=float,
    default=.05,
    show_default=True,
)
@optgroup.option(
    "--rthresh",
    help="Relative score threshold.",
    type=float,
    default=.8,
    show_default=True,
)
@optgroup.option(
    "--taxon",
    help="Taxon(s) to use.  [default: all]",
    multiple=True,
    default=taxons,
)

def main(**params):

    # Scan sequence
    scan_sequence(params["fasta_file"], params["profiles_dir"],
        params["dummy_dir"], params["output_dir"], params["threads"],
        params["background"], params["latest"], set(params["profile"]),
        params["pthresh"], params["rthresh"], params["taxon"])

def scan_sequence(fasta_file, profiles_dir, dummy_dir="/tmp/", output_dir="./",
    threads=1, background=(.25, .25, .25, .25), latest=False, profile=set(), 
    pthresh=.05, rthresh=.8, taxon=taxons):

    # Initialize
    A, C, G, T = background

    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get profiles with which to scan sequence
    profiles = _get_profiles(profiles_dir, latest, profile, taxon)

    # Get profile names
    with open(os.path.join(profiles_dir, "names.json")) as handle:
        names = json.load(handle)

    # Scan profiles against sequence
    _scan_profiles(profiles, fasta_file, names, dummy_dir, output_dir, threads,
        A, C, G, T, pthresh, rthresh)

def _get_profiles(profiles_dir, latest=False, profile=set(), taxon=taxons):

    # Initialize
    profiles = []
    profiles_dict = {}

    # For each taxon...
    for t in taxon:

        # Initialize
        taxon_dir = os.path.join(os.path.abspath(profiles_dir), t)

        # For each profile...
        for profile_file in sorted(os.listdir(taxon_dir), reverse=True):

            # Ignore profiles
            if len(profile) > 0:
                if profile_file[:8] not in profile:
                    continue

            # Initialize key
            key = profile_file[:6]
            profiles_dict.setdefault(key, [])

            # Skip profile if only using the latest version of each profile
            if latest:
                if len(profiles_dict[key]) == 1:
                    continue

            # Add profile
            profiles_dict[key].append(os.path.join(taxon_dir, profile_file))

    # Create list of profiles
    for value_list in profiles_dict.values():
        for p in value_list:
            profiles.append(p)

    return(profiles)

def _scan_profiles(profiles, fasta_file, names, dummy_dir="/tmp/",
    output_dir="./", threads=1, A=.25, C=.25, G=.25, T=.25, pthresh=.05,
    rthresh=.8):

    # Parallelize scanning
    kwargs = {"total": len(profiles), "ncols": 100}
    pool = Pool(threads)
    p = partial(_scan_profile, fasta_file=fasta_file, names=names,
        dummy_dir=dummy_dir, output_dir=output_dir, threads=threads,
        A=A, C=C, G=G, T=T, pthresh=pthresh, rthresh=rthresh)
    for _ in tqdm(pool.imap(p, profiles), **kwargs):
        pass
    pool.close()
    pool.join()

def _scan_profile(profile_file, fasta_file, names, dummy_dir="/tmp/",
    output_dir="./", threads=1, A=.25, C=.25, G=.25, T=.25, pthresh=.05,
    rthresh=.8):

    # Initialize
    cutoff = None
    matrix_id = os.path.basename(profile_file)[:8]
    bin_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bin")
    tsv_file = os.path.join(dummy_dir, "%s.%s.%s" % \
        (os.path.basename(__file__), pid, matrix_id))
    gzipped_file = "%s.gz" % tsv_file
    output_file = os.path.join(output_dir, "%s.tsv.gz" % matrix_id)

    # Skip if profile already scanned (i.e. for speed)
    output_file = os.path.join(output_dir, "%s.tsv.gz" % matrix_id)
    if not os.path.exists(output_file):

        # Calculate distribution of PWM scores
        cmd = "%s %s" % (os.path.join(bin_dir, "matrix_prob"), profile_file)
        process = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)

        with open(tsv_file, "w") as f:

            for line in process.stdout.decode("utf-8").split("\n"):

                matches = re.findall("(\S+)", line)

                if len(matches) == 3:

                    score = matches[0]
                    p_value = float(matches[1])
                    perc = float(matches[2][:-1])

                    f.write(f"{score}\t{int(perc * 10)}\t" +\
                            f"{int(log(p_value) * 1000 / -10)}\n")

                    # Get PWM score cutoff
                    if cutoff is None:
                        cutoff = score
                    elif p_value < pthresh and perc >= rthresh * 100:
                        cutoff = score

        # Scan FASTA file (ugly code but very efficient)
        cmd_1 = "%s -m %s -c %s %s" % (os.path.join(bin_dir, "matrix_scan"),
            profile_file, cutoff, fasta_file)
        cmd_2 = "gzip > %s" % gzipped_file
        cmd = '''%s | awk -v score_tab="%s" -v name="%s" 'BEGIN { while((getline line < score_tab) > 0 ) {split(line,f," "); scores[f[1]]=f[2]; pvalues[f[1]]=f[3]} close(score_tab) } {print $1"\t"$2"\t"$3"\t"name"\t"scores[$5]"\t"pvalues[$5]"\t"$6}' | %s''' % \
            (cmd_1, tsv_file, names[matrix_id], cmd_2)
        subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT)

        # Write output
        shutil.copy(gzipped_file, output_file)

        # Remove dummy files
        os.remove(tsv_file)
        os.remove(gzipped_file)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()