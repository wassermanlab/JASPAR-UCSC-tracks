#!/usr/bin/env perl

# FILE: jasparconvert
# CREATE DATE: 8/27/2013
# AUTHOR: Giovanna Ambrosini 
#
# Part of the code is based on an implemetation by
# William Stafford Noble and Timothy L. Bailey
# Created in 1999
# ORIG: transfac2meme.pl

# Giovanna Ambrosini 18/10/2017 
# Add pseudo weight fraction to correct frequencies
# 

use Math::Round;

# Set up global variables. Assume uniform.
my @bases = ("A", "C", "G", "T");
my $num_bases = 4;
my %bg = ();
$bg{"A"} = 0.25;
$bg{"C"} = 0.25;
$bg{"G"} = 0.25;
$bg{"T"} = 0.25;
my $c = 0;				# default pseudocount fraction
my $logscl = 100;
my $minscore = -10000;

my $usage = "USAGE: jasparconvert.pl [options] <matrix file>

  Options: -species <name>              not used yet
	   -bg <background file>	set of f_a
	   -c <pseudo weight>	        add an arbitrary pseudo weight fraction <c> to each freq
					default: $c
           -m <low value score>         set low value score
                                        default: $minscore
           -n <log scaling factor>      set log scaling factor (int)
                                        default: 100
           -noheader                    write raw matrix (without header)
           -o <outfile>                 output file name
                                        default: no output file
	   -w <l>|<p>                   <l> select log-odds matrix format for written output
	                                <p> select letter-probability matrix format for written output
                                        default format: original JASPAR frequency matrix

  Convert a JASPAR matrix file to MEME format (i.e. integer log likelihoods or letter-probability matrix).
  \n";

my $logodds = 0;
my $letprob = 0;
my $defout  = 0;
my $out_file = "";
my $ofile = 0;

# Process command line arguments.
if (scalar(@ARGV) == 0) {
  printf(STDERR $usage);
  exit(1);
}
my $species = "";
my $header = 1;
while (scalar(@ARGV) > 1) {
  $next_arg = shift(@ARGV);
  if ($next_arg eq "-species") {
    $species = shift(@ARGV);
  } elsif ($next_arg eq "-bg") {
    $bg_file = shift(@ARGV);
  } elsif ($next_arg eq "-c") {
    $c = shift(@ARGV);
  } elsif ($next_arg eq "-m") {
    $minscore = shift(@ARGV);
  } elsif ($next_arg eq "-n") {
    $logscl = shift(@ARGV);
  } elsif ($next_arg eq "-noheader") {
    $header = 0;
  } elsif ($next_arg eq "-o") {
    $out_file = shift(@ARGV);
    $ofile = 1;
  } elsif ($next_arg eq "-w") {
    my $f = shift(@ARGV);
    if ($f eq "l") {
      $logodds = 1;
    } elsif ($f eq "p") {
      $letprob = 1;
    } else {
      print(STDERR "Illegal argument ($next_arg) $f\n");
      exit(1);
    }
  } else {
    print(STDERR "Illegal argument ($next_arg)\n");
    exit(1);
  }
}
($matrix_file) = @ARGV;

if (($logodds or $letprob) and $ofile == 0) {
  print(STDERR "Please, give a filename for your output file (-o <outfile>)\n");
  exit(1);
}
if ($logodds == 0 and $letprob == 0 and $ofile == 1) {
  $defout = 1;
}
# read the background file
if (defined($bg_file)) {
  open($bg_file, "<$bg_file") || die("Can't open $bg_file.\n");
  $total_bg = 0;
  while (<$bg_file>) {
    next if (/^#/);			# skip comments
    ($a, $f) = split;
    if ($a eq "A" || $a eq "a") {
      $bg{"A"} = $f; 
      $total_bg += $f;
    } elsif ($a eq "C" || $a eq "c") {
      $bg{"C"} = $f; 
      $total_bg += $f;
    } elsif ($a eq "G" || $a eq "g") {
      $bg{"G"} = $f; 
      $total_bg += $f;
    } elsif ($a eq "T" || $a eq "t") {
      $bg{"T"} = $f; 
      $total_bg += $f;
    }
  }
  # make sure they sum to 1
  foreach $key (keys %bg) {
    $bg{$key} /= $total_bg;
    #printf STDERR "$key $bg{$key}\n";
  }
}  # background file

# Open the matrix file for reading.
open(MF, "<$matrix_file") || die("Can't open $matrix_file.\n");
# Open the output file
if ($ofile) {
  open(OF, ">$out_file") || die("Can't open $out_file.\n");
}

# Print the MEME header.
print("ALPHABET= ACGT\n");
print("strands: + -\n");
print("Background letter frequencies (from dataset with add-one prior applied):\n");
printf("A %f C %f G %f T %f\n\n",  $bg{"A"}, $bg{"C"}, $bg{"G"}, $bg{"T"});

# Read the input file.
my $num_motifs = 0;
my $num_skipped = 0;
my $i_motif = 0;
my $width = 0;
my $i_base = 0;
my $num_seqs = 0;
my $curspos = 0;
my $prev_curspos = 0;
while ($line = <MF>) {
  $i_motif = 0;
  $i_base = 0;
  
  #print $line;

  # Check header (new matrix)
  if ($line =~ /^>\s*(\S+)(\s+(\S+))?/ or $line =~ /\s*A\s*\[\s*(.*)\s*\]/ or $line =~ /\s*A\s+(.*)/) {
    if (!defined $2) {
      print $line;
      $matrix_name = "Unknown";
      $matrix_id = "";
      @counts = split(' ', $1);
      $width = scalar(@counts);
      for ($i_motif = 0; $i_motif < $width;  $i_motif++) {
        $motif{$i_base, $i_motif} = shift(@counts);
        #print "  motif{$i_base, $i_motif} : $motif{$i_base, $i_motif}";
      }
      $i_base++;
    } else {
      $matrix_id = $1;
      $matrix_name = $2;
    }
    $matrix_id =~ s/\s+//g;
    $matrix_name =~ s/\s+//g;
    # Read the motif.
    print (STDERR "Read the motif for $matrix_id $matrix_name\n\n");
    while () {
      $line = <MF>;
      $curspos = tell(MF);
      #print "CURSOR POS : $curspos\n";
      
      if (! defined $line ) {
        print "\n";
	last;
      }
      print $line;

      if ($line =~ /\s*([ACGT])\s*\[\s*(.*)\s*\]/ or $line =~ /\s*([ACGT])\s+(.*)/) {
        @counts = split(' ', $2);
      }

      if ($line =~ />/) {
        seek(MF, $prev_curspos, 0);  
        #print("Moving cursor to POS $prev_curspos\n");
	last;
      }
      $width = scalar(@counts);

      # Store the contents of this row.
      for ($i_motif = 0; $i_motif < $width;  $i_motif++) {
	$motif{$i_base, $i_motif} = shift(@counts);
        #print "  motif{$i_base, $i_motif} : $motif{$i_base, $i_motif}";
      }
      #print "\n";
      $i_base++;
      if ($curspos ne $prev_curspos) {
       $prev_curspos = $curspos;
       #print "PREV CURSOR POS : $prev_curspos\n";
      }
    } # END OF WHILE ON SINGLE MATRIX
    if ($defout) { 
      print OF ">MOTIF $matrix_id $matrix_name len=$width\n";
    }
    print "MOTIF $matrix_id $matrix_name len=$width\n";
    for ($i_motif = 0; $i_motif < $width; $i_motif++) {
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
        if ($defout) { 
          printf (OF "%7d ",  $motif{$i_base, $i_motif});
        }
        printf ("%7d ",  $motif{$i_base, $i_motif});
      }
      if ($defout) {
        print OF "\n";
      }
      print "\n";
    }
    print "\n";
    # Convert the motif to frequencies.
    # If $c != 0, add pseudocount fraction $c
    for ($i_motif = 0; $i_motif < $width; $i_motif++) {
      # motif columns may have different counts
      $num_seqs = 0;
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
        $num_seqs += $motif{$i_base, $i_motif};
      }
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
         $motif{$i_base, $i_motif} =
         ($motif{$i_base, $i_motif} + ($bg{$bases[$i_base]} * $num_seqs * $c) ) /
         ($num_seqs * (1 + $c));
      }
    }
#   Old stuff
#   for ($i_motif = 0; $i_motif < $width; $i_motif++) {
#     # motif columns may have different counts
#     $num_seqs = 0;
#     for ($i_base = 0; $i_base < $num_bases; $i_base++) {
#       $num_seqs += $motif{$i_base, $i_motif};
#     }
#     for ($i_base = 0; $i_base < $num_bases; $i_base++) {
#	$motif{$i_base, $i_motif} = 
#          ($motif{$i_base, $i_motif} + ($b * $bg{$bases[$i_base]}) ) / 
#          ($num_seqs + $b);
#     }
#   }
    ###### Decide whether to print the motif.

    # If no criteria are given, then print it.
    $print_it = 1;

    # If we were given a species.
    if ($species ne "") {
      # is this the right species?
      $print_it = ($this_species =~ m/$species/);
      if ($this_species eq "") {
	print(STDERR "Warning: No species given for $matrix_name.\n");
      }
    }

    # Print the motif.
    if ($print_it) {
      $num_motifs++;
      print(STDERR "Printing motif $matrix_name.\n\n");

      # PSSM
      if ($logodds) {
        if ($header) {
          print(OF ">log-odds matrix $matrix_id $matrix_name: alength= $num_bases w= $width n= 0 bayes= 0 E= 0\n");
        }
      }
      print("log-odds matrix $matrix_id $matrix_name: alength= $num_bases w= $width n= 0 bayes= 0 E= 0\n");
      for ($i_motif = 0; $i_motif < $width; $i_motif++) {
	for ($i_base = 0; $i_base < $num_bases; $i_base++) {
          if ($logodds) {
            if ($motif{$i_base, $i_motif}) {
	      printf(OF "%7d ", round((log( $motif{$i_base, $i_motif} / $bg{$bases[$i_base]} )/log(2.0))*$logscl) );
            } else {
              printf(OF "%7d ", $minscore);
            }
          }
          if ($motif{$i_base, $i_motif}) { #log2
	    printf("%7d ", round((log( $motif{$i_base, $i_motif} / $bg{$bases[$i_base]} )/log(2.0))*$logscl) );
          } else {
            printf("%7d ", $minscore);
          }
	}
        if ($logodds) {
	  print(OF "\n");
        }
	print("\n");
      }
      print("\n");

      # Letter-probability matrix (PSFM)
      if ($letprob) {
        if ($header) {
          print(OF ">letter-probability matrix $matrix_id $matrix_name: ");
          print(OF "alength= $num_bases w= $width nsites= $num_seqs E= 0\n");
        }
      }
      print("letter-probability matrix $matrix_id $matrix_name: ");
      print("alength= $num_bases w= $width nsites= $num_seqs E= 0\n");
      for ($i_motif = 0; $i_motif < $width; $i_motif++) {
	for ($i_base = 0; $i_base < $num_bases; $i_base++) {
          if ($letprob) {
	    printf(OF "  %8.6f\t", $motif{$i_base, $i_motif});
          }
	  printf("  %8.6f\t", $motif{$i_base, $i_motif});
	}
        if ($letprob) {
	  print(OF "\n");
        }
	print("\n");
      }
      print("\n");
    } else {
      $num_skipped++;
      #print(STDERR "Skipping motif $matrix_name.\n");
    }
  }
}
print(STDERR "Converted $num_motifs motifs.\n");
print(STDERR "Skipped $num_skipped motifs.\n");

close(MF);
if ($ofile) {
  close(OF);
}
