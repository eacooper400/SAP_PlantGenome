#!/usr/bin/perl -w
#
# hmp2linkage.pl
#
# Convert hapmap format to the linkage format used by haploview
#
# September 30, 2014
# Liz Cooper

use strict;

# From the command line arguments, collect the following information:
# The input and output file names
my ($USAGE) = "$0 <input.hmp.txt> <output.linkage> <output.markerInfo>
\tinput.hmp.txt = An input file in the TASSEL hapmap format
\toutput.linkage.txt = The name of the linkage format file to be created
\toutput.markerInfo = The name of a locus info. file to be created\n\n";
unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output, $locusfile) = @ARGV;

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
my @pos_IDs = ();

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
      $pi =~ s/\.//g;
      $pi =~ s/_//g;
      $genotypes{$pi} = '';
      $pi_index{$n} = $pi;
    }
    next FILELOOP;
  }
  # For the remaining lines in the file, get all of the genotypes at that position
  # Store the genotypes for each inidividual, and save the position in a sep. array
  my @bases = @info[11..(scalar @info - 1)];
  my $longstring = join('', @bases);

  # If the site meets all of the requirements, add the genotype to the individual string in the hash
  # Add the position to the positions list
  push (@positions, $info[3]);
  push (@pos_IDs, $info[0]);
  for (my $i = 0; $i < length $longstring; $i++) {
    my $individual = $pi_index{$i};
    $genotypes{$individual} .= substr($longstring, $i, 1);
  }
}
close(IN);

# Open the output file for printing
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

# For each individual, convert the genotype to 2 alleles, and print all positions to 1 line 
# Make sure to skip the reference
my @inds = sort {$a <=> $b} (keys %pi_index);
foreach my $i (@inds) {
  my $pi = $pi_index{$i};
  unless ($pi =~ /REFERENCE/) {
    my $genstring = $genotypes{$pi};
    my @new_gens = ();

    for (my $g = 0; $g < length $genstring; $g++) {
      my $basecall = substr($genstring, $g, 1);
      my $alleles = '';
      if ($basecall =~ /[MRWSYK]/) {
	my @pair = @{$hets{$basecall}};
	$alleles = $pair[0] . " " . $pair[1];
	push (@new_gens, $alleles);
      } elsif ($basecall =~ /[ACGT]/) {
	$alleles = $basecall . " " . $basecall;
	push (@new_gens, $alleles);
      } elsif ($basecall =~ /N/) {
	$alleles = "0 0";
	push (@new_gens, $alleles);
      } else {
	print $basecall;
      }
    }
    my $newstring = join("\t", @new_gens);
    
    my $ped = $i + 1;
    print OUT $ped, "\t", $pi, "\t", "0\t", "0\t", "1\t", "0\t", $newstring, "\n";
  }
}
close(OUT);

# Open the marker locus file for printing
open (LOCUS, ">$locusfile") || die "\nUnable to open the file $locusfile!\n";
for (my $p = 0; $p < scalar @positions; $p++) {
  print LOCUS $pos_IDs[$p], "\t", $positions[$p], "\n";
}
close(LOCUS);

exit;
