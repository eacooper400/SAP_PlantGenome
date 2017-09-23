#!/usr/bin/perl -w
#
# sam2fastq.pl
#
# Convert SAM format files to fastq files
#
# January 3, 2014
# Liz Cooper

use strict;

# Get the input and output file names from the command line
my ($USAGE) = "\n$0 <input.sam> <output.fastq>
\tinput.sam = The input file in .sam format
\toutput.fastq = The name of the output file to create\n\n";
unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output) = @ARGV;

open (IN, $input) || die "\nUnable to open the file $input!\n";
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

while (<IN>) {
  chomp $_;
  my @info = split(/\s{1,}/, $_);
  print OUT "@", $info[0], "\n";
  print OUT $info[9], "\n";
  print OUT "+", "\n";
  print OUT $info[10], "\n";
}
close(IN);
close(OUT);
exit;
