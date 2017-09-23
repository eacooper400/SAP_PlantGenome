#!/usr/bin/perl -w
#
# fastq2fasta.pl
#
# Convert fastq file to fasta format
#
# December 6, 2012
# Liz Cooper

use strict;

my ($USAGE) = "\n$0 input.fastq
\tinput.fastq = The name of the input .fastq file.  The output .fasta file will automatically have the same prefix as this input file\n\n";
unless(@ARGV) {
  print $USAGE;
  exit;
}
my $input = $ARGV[0];
my ($prefix, $ext) = split(/\./, $input);
my $output = $prefix . ".fasta";

open (IN, $input) || die "\nUnable to open the file $input!\n";
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

my $flag = 0;

while (<IN>) {
  chomp $_;
  if ($_ =~ /^\@/) {
    $_ =~ s/^\@/>/;
    print OUT $_, "\n";
    $flag = 1;
  } elsif ($_ =~ /\+/) {
    $flag = 0;
    next;
  } elsif ($flag) {
    print OUT $_, "\n";
  } else {
    next;
  }
}
close(IN);
close(OUT);
exit;
