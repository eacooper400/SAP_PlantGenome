#!/usr/bin/perl -w
#
# SNP_distance.pl
#
# Calculate the distance between consecutive SNPs in a hapmap file
#
# January 20, 2015
# Liz Cooper

use strict;

# Get the name of the input hapmap file and the output distance file from the command line
my ($USAGE) = "\n$0 <input.hmp.txt> <output.txt>
\tinput.hmp.txt = An input file in hapmap format
\toutput.txt = An output file with a single column of distances\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output) = @ARGV;

open (IN, $input) || die "\nUnable to open the file $input!\n";
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

my $prev_chr = 0;
my $prev_pos = 0;

while (<IN>) {
  chomp $_;
  
  if ($_ =~ /^rs/) {
    next;
  }

  my @info = split(/\s{1,}/, $_);
  my $chr = $info[2];
  my $pos = $info[3];

  if ($chr == $prev_chr) {
    my $distance = $pos - $prev_pos;
    print OUT $distance, "\n";
    $prev_pos = $pos;
  } else {
    $prev_chr = $chr;
    $prev_pos = $pos;
  }
}
close(IN);
close(OUT);
exit;
