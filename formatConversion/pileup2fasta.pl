#!/usr/bin/perl -w
#
# pileup2fasta.pl
#
# Convert pileup results to a fasta file to use in building a reference
#
# December 11, 2012
# Liz Cooper

use strict;

# Get the names of the input and output files from the command line
my ($USAGE) = "\n$0 <input.pileup> <output.fasta>
\tinput.pileup = An input file in the pileup format created by Samtools
\toutput.fasta = The name of the output file in .fasta format to create\n\n";
unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output) = @ARGV;

open (IN, $input) || die "\nUnable to open the file $input!\n";
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

my $chr = '';
my $string = '';
my $count = 0;
my $prev_pos = 0;

while (<IN>) {
  chomp $_;
  my @info = split(/\t/, $_);
  my @temp = split("_", $info[0]);
  my $name = pop @temp;
  if ($count == 0) {
    $chr = $name;
    $string .= $info[2];
    $prev_pos = $info[1];
    $count++;
  } elsif ($chr eq $name) {
    if (($info[1] == ($prev_pos +1)) && (length($string) < 94)) {
      $string .= $info[2];
      $prev_pos = $info[1];
    } else {
      print OUT ">", $chr, "_", $count, "\n";
      print OUT $string, "\n";
      $count++;
      $string = $info[2];
      $prev_pos = $info[1];
    }
  } else {
    print OUT ">", $chr, "_", $count, "\n";
    print OUT $string, "\n";
    $chr = $name;
    $string = $info[2];
    $prev_pos = $info[1];
    $count = 1;
  }
}
print OUT  ">", $chr, "_", $count, "\n";
print OUT $string, "\n";

close(IN);
close(OUT);
exit;
