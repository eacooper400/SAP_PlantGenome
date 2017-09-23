#!/usr/bin/perl -w
#
# hmp2numeric.pl
#
# Convert SNPs in IUPAC format (A,C,G,T + het codes) to numerical format
# Where 2 = homozygous minor allele, 0 = homozygous major allele, and 1 = heterozygous
#
# March 5, 2015
# Liz Cooper

use strict;

# Take as input the text file with the SNPs, and the name of an output file
my ($USAGE) = "\n$0 <input.hmp.txt> <output.txt>
\tinput.hmp.txt = The name of the input file in hapmap format
\toutput.txt = The name of the output file to create with numeric SNP codes\n\n";
unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output) = @ARGV;

# Open the output file for printing
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

# Open the input file, and start reading through one line at a time
open (IN, $input) || die "\nUnable to open the file $input!\n";

while (<IN>) {
  chomp $_;

  # Print the first line of the input file to the output, without changing anything
  if ($. == 1) {
    print OUT $_, "\n";
    next;
  }

  # If the line contains the allele possibilities, store the information in an array,
  # and then print the same line back to the output file
  else {
    my @temp = split(/\s{1,}/, $_);
    my $minor_allele = (split(/\//, $temp[1]))[1];
    my @genotypes = @temp[11..(scalar @temp - 1)];
    for (my $g = 0; $g < scalar @genotypes; $g++) {
      if ($genotypes[$g] =~ /N/) {
	$genotypes[$g] = '-9';
      } elsif ($genotypes[$g] =~ /[MRWSYK]/) {
	$genotypes[$g] = '1';
      } elsif ($genotypes[$g] =~ /$minor_allele/) {
	$genotypes[$g] = '2';
      } else {
	$genotypes[$g] = '0';
      }
    }
    my $new1 = join("\t", @temp[0..10]);
    my $new2 = join("\t", @genotypes);
    print OUT $new1, "\t", $new2, "\n";
  }
}
close(IN);
close(OUT);
exit;
