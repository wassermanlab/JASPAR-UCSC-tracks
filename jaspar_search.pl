#!/usr/bin/env perl

=head1 NAME

jaspar_search.pl

=head1 SYNOPSIS

  jaspar_search.pl
            -f fasta_file
            (-id matrix_id | -m matrix_file | -n matrix_name)
            [-db jaspar_db]
            [-th threshold]
            [-o out_file]

=head1 ARGUMENTS
 
  -f fasta_file     = FASTA file containing the sequence to be scanned
  -id matrix_id     = JASPAR ID of TFBS profile with which to search the
                      sequence(s)
  -m matrix_file    = JASPAR TFBS profile with which to searc the sequence(s)
  -n matrix_name    = Name of the JASPAR TFBS profile with which to search
                      the sequence(s)
  -db               = JASPAR database name. DEFAULT = JASPAR_2016
  -th threshold     = Minimum score of putative TFBS hit to report.
                      Specify as an absolute score, i.e. 14.1 or as
                      a relative score, i.e. 75%.
                      DEFAULT = 75%
  -o out_file       = Output file to which to write TFBS hits. If not
                      provided, write to stdout.

=head1 DESCRIPTION

Retrieve the specified JASPAR TFBS profile from the database and scan the
given sequence in FASTA format. The entire sequence is scanned in chunks
with the given profile. Write out putative TFBS in BED format. The results
are written to the output file if provided or to standard output otherwise.

=head1 AUTHOR

  David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia
  
  E-mail: dave@cmmt.ubc.ca
  
  Edited by Oriol Fornes on July 5, 2017

=cut

use strict;
use warnings;

use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;
use TFBS::Matrix::PFM;
use TFBS::DB::JASPAR;

use constant JASPAR_DB_NAME => "JASPAR_2016";
use constant JASPAR_DB_HOST => "vm5.cmmt.ubc.ca";
use constant JASPAR_DB_USER => "jaspar_r";
use constant JASPAR_DB_PASS => "";

use constant TFBS_THRESHOLD => "75%";

my $fasta_file;
my $jaspar_db;
my $matrix_id;
my $matrix_file;
my $matrix_name;
my $threshold;
my $out_file;
GetOptions(
    'f=s'     => \$fasta_file,
    'db=s'    => \$jaspar_db,
    'id=s'    => \$matrix_id,
    'm=s'     => \$matrix_file,
    'n=s'     => \$matrix_name,
    'th=s'    => \$threshold,
    'o=s'     => \$out_file
);

$jaspar_db  = JASPAR_DB_NAME if !$jaspar_db;
$threshold  = TFBS_THRESHOLD if !defined $threshold;

unless ($matrix_id || $matrix_file || $matrix_name) {
    pod2usage(
        -verbose => 1,
        -msg => "Please provide either a matrix name or ID"
    );
}

my $pfm;

if ($matrix_file) {
    $pfm = TFBS::Matrix::PFM->new('-matrixfile', $matrix_file);
} else {
    my $db = TFBS::DB::JASPAR7->connect(
        "dbi:mysql:" . $jaspar_db . ":" . JASPAR_DB_HOST,
        JASPAR_DB_USER, JASPAR_DB_PASS);
    die "Error connecting to JASPAR database\n" if !$db;

    if ($matrix_id) {
        $pfm = $db->get_Matrix_by_ID($matrix_id, 'PFM');
        die "Error retrieving TFBS matrix from JASPAR\n" if !$pfm;
    } elsif ($matrix_name) {
        $pfm = $db->get_Matrix_by_name($matrix_name, 'PFM');
        die "Error retrieving TFBS matrix from JASPAR\n" if !$pfm;
    }
}

my $pwm = $pfm->to_PWM;
die "Error converting PFM to PWM\n" if !$pwm;

# Amount to overlap successive chromsome segments. Overlap by 1 nt less than
# the width of the motif to avoid duplicate hits in the overlapping regions.
my $motif_len = $pwm->length;
my $chrom_seg_overlap = $motif_len - 1;

if ($threshold !~ /%$/ && $threshold == 0) {
    my $min_score = $pwm->min_score;
    my $max_score = $pwm->max_score;

    #
    # Compute relative pct. threshold based on a raw score of 0.
    # This is neccesary as the search function appears to consider a
    # threshold of 0 as a percent threshold even if no "%" is given.
    #
    my $rel_threshold = -($min_score / ($max_score - $min_score)) * 100;

    $threshold = "$rel_threshold%";
}

if ($out_file) {
    open(OFH, ">$out_file")
        || die "Error opening output file $out_file\n";
} else {
    open(OFH, ">-")
        || die "Error opening standard output\n";
}

my $seqio = Bio::SeqIO->new(-file => $fasta_file, '-format' => 'Fasta');
# While sequences... #
while(my $seq = $seqio->next_seq) {
    # do stuff with $string
    my $chrom = $seq->display_id;
    my $chrom_size = $seq->length;

    my $siteset = $pwm->search_seq(-seq => $seq->seq, -threshold => "$threshold");

    if ($siteset && $siteset->size > 0) {
        my $it = $siteset->Iterator(-sort_by => 'start');
        while (my $site = $it->next) {
            my $strand = $site->strand ? ($site->strand > 0 ? '+' : '-') : '?';

            printf OFH "$chrom\t%d\t%d\t%s\t%.3f\n",
            $site->start,
            $site->end,
            $strand,
            $site->rel_score;
        }
    }
}

close(OFH);

exit;
