#!/usr/bin/perl -w
#
# missing_hmp.pl
#
# Calculate percent missingness for every individual in a hapmap format file
#
# July 14, 2014
# Liz Cooper

use strict;

# From the command line collect:
# The input file name in hapmap format
# The output file name in text format
my ($USAGE) = "\n$0 <input.hmp> <output.txt>
\tinput.hmp = The input file in hapmap format
\toutput.txt = A tab-delimited file with 2 columns:
\t\tColumn1: Individual/Sample Name
\t\tColumn2: The fraction of missing sites for that individual\n\n";
unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output) = @ARGV;

# Create 2 arrays to save the individual names and the number of Ns for the output file
# Save the total count of SNP sites outside of the file loop
my @individuals = ();
my @num_N = ();
my $total_Sites = 0;

# Read through the input file, store the individual names from the first line
# Then tally up the Ns for each individual at all subsequent lines
open (IN, $input) || die "\nUnable to open the file $input!\n";
while (<IN>) {
  chomp $_;
  my @info = split(/\t/, $_);
  
  # Check if this is the first line
  # If so, get the list of PI numbers, and initialize the numN array with zeros
  if ($_ =~ /^rs/) {
    my @names = @info[12..(scalar @info - 1)];
    for (my $n = 0; $n < scalar @names; $n++) {
      my $pi = (split/\:/, $names[$n])[0];
      push (@individuals, $pi);
      $num_N[$n] = 0;
    }
    next;
  }

  # Otherwise, get the genotype for each individual, 
  # and see if it is an N
  else {
    my @genotypes = @info[12..(scalar @info - 1)];
    $total_Sites += 1;
    for (my $g = 0; $g < scalar @genotypes; $g++) {
      if ($genotypes[$g] =~ /N/) {
	$num_N[$g] += 1;
      }
    }
  }
}
close(IN);

# Open the outfile to print the results for each individual
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

# For each individual, calculate the fraction missing, and print to the out file with the indiviudal name
for (my $i = 0; $i < scalar @individuals; $i++) {
  my $frac = $num_N[$i]/$total_Sites;
  print OUT $individuals[$i], "\t", $frac, "\n";
}
close(OUT);
exit;
