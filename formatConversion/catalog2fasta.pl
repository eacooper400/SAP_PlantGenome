#!/usr/bin/perl -w
#
# catalog2fasta.pl
#
# Convert catalog list of tags (from Stacks) to a fasta file
#
# August 13, 2013
# Liz Cooper

use strict;

# Get the name of the input and output files from the command line
my ($USAGE) = "\n$0 <input.tsv> <output.fasta>
\tinput.tsv = A tsv format file create by Stacks
\toutput.fasta = The name of the fasta format file to create\n\n";
unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output) = @ARGV;

# Open the input and output files, and start processing the input
open (IN, $input) || die "\nUnable to open the file $input!\n";
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

while (<IN>) {
  chomp $_;
  my @info = split(/\t/, $_);
  print OUT ">", $info[2], "\n";
  print OUT $info[9], "\n";
}
close(IN);
close(OUT);
exit;
