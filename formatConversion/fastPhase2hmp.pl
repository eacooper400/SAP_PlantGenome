#!/usr/bin/perl -w
#
# fastPhase2hmp.pl
#
# Convert fastPHASE output into hapmap format
#
# August 15, 2014
#
# Elizabeth A. Cooper

use strict;
use Data::Dumper;

# Take as input:
# The original fastPHASE input file (need this to get the positions),
# The genotypes out file from PHASE
# The chromosome
# A name for the new hapmap file to create
my ($USAGE) = "\n$0 <phase.inp> <genotypes.out> <chr> <Remove Low Confidence? (Y/N)> <new.hmp.txt>\n
\tphase.inp = The  original fastPHASE input file (need this to get the positions),
\tgenotypes.out = The genotypes out file from fastPHASE,
\tchr = Chromosome Number (only 1 chromosome at a time)
\tRemove Low Confidence = Enter Y if you want to replace bracketed (low confidence) genotypes with missing values; enter N to keep all imputed genotype calls regardless of confidence.
\tnew.hmp.txt = A name for the new hapmap file to create\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}

my ($phasein, $genin, $chr, $remove, $output) = @ARGV;

# First, open the original phase input file and get all of the positions in a list
print "Gathering Position Info. from the Original Input File...\n";

my @positions = ();
open (PHASEIN, $phasein) || die "\nUnable to open the file $phasein!\n";
while (<PHASEIN>) {
  chomp $_;
  if ($. == 3) {
    $_ =~ s/^P\s{1,}//;
    @positions = split(/\s{1,}/, $_);
    last;
  } else {
    next;
  }
}
close(PHASEIN);
print "Getting all of the Genotypes at Every Position...\n";

# Open the output file and print the first part of the header
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";
print OUT "rs\t", "alleles\t", "chrom\t", "pos\t", "strand\t", "assembly\t", "center\t", "protLSID\t", "assayLSID\t", "panelLSID\t", "QCcode";

# Create a hash to save the genotypes for each position
# Print out PI numbers to the header of the hapmap file as they are read in
# Also keep track of the number of PIs
my $num_ind = 0;
my %genotypes = ();

open (IN, $genin) || die "\nUnable to open the file $genin!\n";
my $gflag = 0;
my $counter = 0;

while (<IN>) {
  chomp $_;

  # Skip all of the lines before the genotypes start
  # Skip the lines that say "BEGIN GENOTYPES" and "END GENOTYPES"
  if ($_ =~ /BEGIN\sGENOTYPES/) {
    $gflag = $.;
    next;
  } elsif ($_ =~ /END GENOTYPES/) {
    $gflag=0;
    last;
  } elsif ($gflag) {
    my $check = $. - $gflag;
    if (($counter == 0) || ($check == ($counter+3))) {
      print OUT "\t", $_;
      $num_ind++;
      $counter = $check;
    } else {
      # Save the genotypes in a hash
      my $string = $_;
      $string =~ s/\s//ig;
      if ($remove =~ /Y/) {
	$string =~ s/\[[A-Z]\]/N/ig;
      } else {
	$string =~ s/(\[|\])//ig;
      }
      for (my $i = 0; $i < (length($string)); $i++) {
	my $base = substr($string, $i, 1);
	if (exists $genotypes{$positions[$i]}) {
	  $genotypes{$positions[$i]} .= $base;
	} else {
	  $genotypes{$positions[$i]} = $base;
	}
      }
    } 
  } else {
    next;
  }
}
close(IN);
print OUT "\n";
# Now, for each position, print out the necessary info to the output file,
# line by line
print "Converting Haplotype Calls to Diploid Genotypes...\n";
LOOP:for (my $p = 0; $p < scalar @positions; $p++) {
  
  # Check that the number of genotypes matches the number of individuals
  my $string = $genotypes{$positions[$p]};
  if (((length $string)/2) > $num_ind) {
    print length $string, "\t", $num_ind, "\n";
    next LOOP;
  }
  
  print OUT "S", $chr, "_", $positions[$p], "\t";
  
  # Figure out the two alleles at the site
  my %alleles = ();
  for (my $g = 0; $g < length $string; $g++) {
    my $tmp = substr($string, $g, 1);
    unless ($tmp =~ /N/) {
      if (exists $alleles{$tmp}) {
	$alleles{$tmp} += 1;
      } else {
	$alleles{$tmp} = 1;
      }
    }
  }
  my @alt = keys %alleles;
  my $astring = join("/", @alt);
  print OUT $astring, "\t", $chr, "\t", $positions[$p], "\t", "+", "\t", "NA\t", "NA\t", "NA\t", "NA\t","NA\t","NA\t";

  # Now print out all of the genotypes
  my @gens = ();
  for (my $s = 0; $s < ((length $string) - 1); $s+=2) {
    my $diploid = substr($string, $s, 2);
    my $single = get_iupac($diploid);
    push (@gens, $single);
  }
  my $outline = join("\t", @gens);
  print OUT $outline, "\n";
}
print "Done\n";
close(OUT);
exit;

###########################################################################
# Subroutine(s)
###########################################################################
sub get_iupac {
  my ($code) = @_;

  my %iupac = (
	       'AA' => 'A',
	       'AC' => 'M',
	       'AG' => 'R',
	       'AT' => 'W',
	       'AN' => 'N',
	       'CA' => 'M',
	       'CC' => 'C',
	       'CG' => 'S',
	       'CT' => 'Y',
	       'CN' => 'N',
	       'GA' => 'R',
	       'GC' => 'S',
	       'GG' => 'G',
	       'GT' => 'K',
	       'GN' => 'N',
	       'TA' => 'W',
	       'TC' => 'Y',
	       'TG' => 'K',
	       'TT' => 'T',
	       'TN' => 'N',
	       'NA' => 'N',
	       'NG' => 'N',
	       'NT' => 'N',
	       'NN' => 'N'
	       );

  my $new = $iupac{$code};
  return($new);
}
