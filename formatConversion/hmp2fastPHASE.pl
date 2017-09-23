#!/usr/bin/perl -w
#
# hmp2fastPHASE.pl
#
# Convert sites from TASSEL hapmap format into fastPHASE input format
#
# August 5, 2014
#
# Elizabeth A.Cooper

use strict;
use Data::Dumper;

# From the command line arguments, collect the following information:
# The input and output file names
my ($USAGE) = "$0 <input.hmp.txt> <output.phase.inp>\n
\tinput.hmp.txt = An input file in the TASSEL hapmap format
\toutput.phase.inp = The name of the fastPHASE input file to be created\n";
unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output) = @ARGV;

# Create a hash to save het code information
my %hets = (
	    'M' => ['A', 'C'],
	    'R' => ['A', 'G'],
	    'W' => ['A', 'T'],
	    'S' => ['G', 'C'],
	    'Y' => ['T', 'C'],
	    'K' => ['G', 'T']
	    );

# Create a hash keyed by individual ID to save all of the genotypes at every saved position
# Create another hash keyed the individual order/index to keep track of PI number orders
# Create an array to save all of the position IDs
my %genotypes = ();
my %pi_index = ();
my @positions = ();

# Open the input hapmap file, and start processing each line
open (IN, $input) || die "\nUnable to open the file $input!\n";

FILELOOP: while (<IN>) {
  chomp $_;
  my @info = split(/\s{1,}/, $_);

  # Get all of the individual names from the first line of the file
  if ($_ =~ /^rs/) {
    my @names = @info[11..(scalar @info - 1)];
    for (my $n = 0; $n < scalar @names; $n++) {
      my $pi = (split/\:/, $names[$n])[0];
      $genotypes{$pi} = '';
      $pi_index{$n} = $pi;
    }
    next FILELOOP;
  }

  # For the remaining lines in the file, get all of the genotypes at that position
  # Store the genotypes for each inidividual, and save the position in a sep. array
  my @bases = @info[11..(scalar @info - 1)];
  my $longstring = join('', @bases);

  # Check that the site has only 2 alleles
  my %temp_alleles = ();
  for (my $pos = 0; $pos < length $longstring; $pos++) {
    my $base = substr($longstring, $pos, 1);
    if ($base =~ /[ACGT]/) {
      if (exists $temp_alleles{$base}) {
	$temp_alleles{$base} += 1;
      } else {
	$temp_alleles{$base} = 1;
      }
    } elsif ($base =~ /[MRWSYK]/) {
      my @codes = @{$hets{$base}};
      foreach my $code (@codes) {
	if (exists $temp_alleles{$code}) {
	$temp_alleles{$code} += 1;
      } else {
	$temp_alleles{$code} = 1;
      }
      }
    }
  }
  unless ((scalar (keys %temp_alleles)) > 2) {
    
    # If the site meets all of the requirements, add the genotype to the individual string in the hash
    # Add the position to the positions list
    push (@positions, $info[3]);
    for (my $i = 0; $i < length $longstring; $i++) {
      my $individual = $pi_index{$i};
      $genotypes{$individual} .= substr($longstring, $i, 1);
    }
  }
}
close(IN);

# Open the output file for printing
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

# First, count the number of individuals (in the PI index) and print that to the first line of the file (make sure to skip the reference genome)
my $numIND = 0;
my @inds = sort {$a <=> $b} (keys %pi_index);
foreach my $i (@inds) {
  my $pi = $pi_index{$i};
  unless ($pi =~ /REFERENCE/) {
    $numIND++;
  }
}
print OUT $numIND, "\n";

# Count the number of SNP positions, and print that to the output file
my $numSNPs = scalar @positions;
print OUT $numSNPs, "\n";

# Print out the list of positions (in numerical order along the chromosomes) preceded by the letter P
print OUT "P";
foreach my $p (@positions) {
  print OUT " ", $p;
}
print OUT "\n";

# Now go through the individual PIs, get the genotypes, convert het. calls and print 2 lines to the output file
foreach my $i (@inds) {
  my $pi = $pi_index{$i};
  unless ($pi =~ /REFERENCE/) {
    my $genstring = $genotypes{$pi};
    my $hstring1 = '';
    my $hstring2 = '';

    for (my $g = 0; $g < length $genstring; $g++) {
      my $basecall = substr($genstring, $g, 1);
      if ($basecall =~ /[MRWSYK]/) {
	my @pair = @{$hets{$basecall}};
	$hstring1 .= $pair[0];
	$hstring2 .= $pair[1];
      } elsif ($basecall =~ /N/) {
	$hstring1 .= '?';
	$hstring2 .= '?';
      } else {
	$hstring1 .= $basecall;
	$hstring2 .= $basecall;
      }
    }
    print OUT $pi, "\n";
    print OUT $hstring1, "\n";
    print OUT $hstring2, "\n";
  }
}
close(OUT);
exit;
