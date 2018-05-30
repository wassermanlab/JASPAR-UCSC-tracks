#!/usr/bin/env python2.7
import os, re
from Bio import motifs
from numpy import log10 as log
import optparse
import shutil
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

    parser = optparse.OptionParser("./%prog -f <fasta_file> -j <matrix_id> -p <profiles_dir> -s <pwmscan_dir> [-b <background> --dummy=<dummy_dir> -o <output_dir> --pv-thresh=<p_value_thresh> --rs-thresh=<rel_score_thresh>]")

    parser.add_option("-b", "--background", default="0.25,0.25,0.25,0.25", action="store", dest="background", help="Background frequency for A, C, G, and T in the genome (e.g. for hg38 use \"0.29,0.21,0.21,0.29\"; default = \"0.25,0.25,0.25,0.25\")", metavar="<background>")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="<dummy_dir>")
    parser.add_option("-f", action="store", type="string", dest="fasta_file", help="FASTA file (e.g. hg38.fa)", metavar="<fasta_file>")
    parser.add_option("-j", action="store", type="string", dest="matrix_id", help="JASPAR matrix ID (e.g. MA0002.2)", metavar="<matrix_id>")
    parser.add_option("-o", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="<output_dir>")
    parser.add_option("-p", action="store", type="string", dest="profiles_dir", help="Profiles directory (from jaspar2meme.py and jaspar2pfm.py)", metavar="<profiles_dir>")
    parser.add_option("-s", action="store", type="string", dest="pwmscan_dir", help="Full path to PWMScan bin directory (i.e. where all PWMScan executables are located; e.g. $PWMScan/bin)", metavar="<pwmscan_dir>")
    parser.add_option("--pv-thresh", default=0.05, action="store", type="float", dest="p_value_thresh", help="P-value threshold (default = 0.05)", metavar="<p_value_thresh>")
    parser.add_option("--rs-thresh", default=0.8, action="store", type="float", dest="rel_score_thresh", help="Relative score threshold (default = 0.8)", metavar="<rel_score_thresh>")

    (options, args) = parser.parse_args()

    if options.fasta_file is None or options.matrix_id is None or options.profiles_dir is None or options.pwmscan_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")
    
    for i in options.background.split(","):
        try:
            float(i)
        except ValueError:
            parser.error("incorrect value for background frequency: \"%s\"\n\tPlease specify A, C, G and T background frequencies as comma-separated values" % i)

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

    # Skip if output file already exists #
    output_file = os.path.join(os.path.abspath(options.output_dir), "%s.bed.gz" % (options.matrix_id))
    if not os.path.exists(output_file):
        # Initialize #
        cutoff = None
        background = options.background.split(",")
        dummy_bed = os.path.join(os.path.abspath(options.dummy_dir), "%s.%s.bed" % (os.path.basename(__file__), os.getpid()))
        dummy_pwm = os.path.join(os.path.abspath(options.dummy_dir), "%s.%s.pwm" % (os.path.basename(__file__), os.getpid()))
        dummy_tsv = os.path.join(os.path.abspath(options.dummy_dir), "%s.%s.tsv" % (os.path.basename(__file__), os.getpid()))
        # Load JASPAR profile #
        if os.path.exists(os.path.join(os.path.abspath(options.profiles_dir), "%s.pfm" % options.matrix_id)):
            profile_file = os.path.join(os.path.abspath(options.profiles_dir), "%s.pfm" % options.matrix_id)
        elif os.path.exists(os.path.join(os.path.abspath(options.profiles_dir), "%s.jaspar" % options.matrix_id)):
            profile_file = os.path.join(os.path.abspath(options.profiles_dir), "%s.jaspar" % options.matrix_id)
        else: raise ValueError("PFM in JASPAR format does not exist!")
        with open(profile_file) as f:
            profile = motifs.read(f, "jaspar")
        profile.background = {"A": float(background[0]), "C": float(background[1]), "G": float(background[2]), "T": float(background[3])}
        profile.pseudocounts = motifs.jaspar.calculate_pseudocounts(profile)
        # Convert to PWMScan format #
        for i in range(len(profile.pssm["A"])):
            functions.write(dummy_pwm, "\t".join([str(int(profile.pssm[j][i] * 100)) for j in "ACGT"]))
        # Calculate distribution of matrix scores #
        try:
            process = subprocess.check_output([os.path.join(os.path.abspath(options.pwmscan_dir), "matrix_prob"), dummy_pwm], stderr=subprocess.STDOUT)
            for line in process.split("\n"):
                m = re.search("(\S+)\s+(\S+)\s+(\S+)%", line)
                if m:
                    score = m.group(1)
                    p_value = float(m.group(2))
                    perc = float(m.group(3))
                    functions.write(dummy_tsv, "%s\t%s\t%s" % (score, int(perc * 10), int(log(p_value) * 1000 / -10)))
                    if p_value < options.p_value_thresh and perc >= options.rel_score_thresh * 100:
                        cutoff = score
        except:
            raise ValueError("Could not calculate distribution of matrix scores!")
        # Scan DNA sequence for TFBS matches #
        try:
            bash_command = ''' %s -m %s -c %s %s | awk -v score_tab="%s" -v name="%s" 'BEGIN { while((getline line < score_tab) > 0 ) {split(line,f," "); scores[f[1]]=f[2]; pvalues[f[1]]=f[3]} close(score_tab) } {print $1"\t"$2"\t"$3"\t"name"\t"scores[$5]"\t"pvalues[$5]"\t"$6}' | gzip > %s ''' % (os.path.join(os.path.abspath(options.pwmscan_dir), "matrix_scan"), dummy_pwm, cutoff, os.path.abspath(options.fasta_file), dummy_tsv, profile.name, dummy_bed)
            process = subprocess.call(bash_command, shell=True, stderr=subprocess.STDOUT)
        except:
            raise ValueError("Could not scan DNA sequence file for TFBS matches!")
        # Copy dummy file #
        shutil.copy(dummy_bed, output_file)
        # Remove dummy files #
#        if os.path.exists(dummy_bed): os.remove(dummy_bed)
#        if os.path.exists(dummy_pwm): os.remove(dummy_pwm)
#        if os.path.exists(dummy_tsv): os.remove(dummy_tsv)
