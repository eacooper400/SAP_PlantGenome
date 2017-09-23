#!/usr/bin/perl -w
#
# hmpSubset.pl
#
# Extract a specific subset of individuals from a hapmap formatted SNP data file
#
# July 10, 2014
# Liz Cooper

use strict;

# From the command line arguments, collect the following information:
# The input and output file names,
# A text file containing a list of individuals to extract, in a single column
my ($USAGE) = "\n$0 <input.hmp.txt> <output.hmp.txt> <list.txt>
\tinput.hmp.txt = An input file in the TASSEL hapmap format
\toutput.hmp.txt = The name of a hapmap format output file to create
\tlist.txt = A text file containing a list of individuals to extract, in a single column\n\n";
unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output, $list) = @ARGV;

# Create a hash keyed by individual order/index to keep track of PI number orders
my %pi_index = ();

# Open and read in the complete list file
open (LIST, $list) || die "\nUnable to open the file $list!\n";
my @subset_list = <LIST>;
close(LIST);

# Open the output file so that it can be printed to as each line is processed
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

# Open the input file, and get the subset of genotypes one line at a time
open (IN, $input) || die "\nUnable to open the file $input!\n";

FILELOOP: while (<IN>) {
  chomp $_;
  my @info = split(/\s{1,}/, $_);

  # For the first line of the file, find the index position of the subset of desired individuals
  if ($_ =~ /^rs/) {
    my @headers = @info[0..10];
    my $header_line = join("\t", @headers);
    print OUT $header_line;
    my @names = @info[11..(scalar @info - 1)];

    # Check if each name is in the list
    for (my $n = 0; $n < scalar @names; $n++) {
      my $pi = (split/\:/, $names[$n])[0];
      foreach my $search (@subset_list) {
	if ($search =~ /$pi/) {
	  $pi_index{$n} = $pi;
	  print OUT "\t", $names[$n];
	}
      }
    }
    print OUT "\n";
  } else {

    # For the SNP rows, get the subset of genotypes for the saved PI numbers
    my @misc = @info[0..10];
    my $misc_line = join("\t", @misc);
    print OUT $misc_line;
    my @gens = @info[11..(scalar @info - 1)];

    my @saved = sort {$a <=> $b} (keys %pi_index);

    foreach my $s (@saved) {
      print OUT "\t", $gens[$s];
    }
    print OUT "\n";
  }
}
close(IN);
close(OUT);
exit;
