#!/usr/bin/perl -w
#
# filter_sites_hmp.pl
#
# Filter the sites based on missingness and allele frequency in a hapmap file
#
# October 31, 2014
# Liz Cooper

use strict;

# Get the name of the input file, the name of the output file, 
# the min. percent of individuals and the min MAF from the command line
my ($USAGE) = "\n$0 <input.hmp> <output.hmp> <min_percent_ind> <MAF>
\tinput.hmp = The input file in hapmap format
\toutput.hmp = The name of the output file, which will also be hapmap format
\tmin_percent_ind = The minimum fraction of individuals needed to be non-missing to keep a given site (b/t 0 and 1)
\tMAF = The minimum minor allele frequency to retain a site (b/t 0 and 1)\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output, $minCount, $minMAF) = @ARGV;

# Open the input and output files; start processing the input one line at a time
open (IN, $input) || die "\nUnable to open the file $input!\n";
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

while (<IN>) {
  chomp $_;
  if ($_ =~ /^rs/) {
    print OUT $_, "\n";
    next;
  } else {
    my @info = split(/\s{1,}/, $_);
    my @snps = @info[11..(scalar @info - 1)];
    my $numInd = scalar @snps;

    # Check if the number of non-missing individuals is above the threshold
    my $string = join('', @snps);
    $string =~ s/N//g;
    my $non_missing = length $string;
    my $perc_nonMissing = $non_missing/$numInd;
    if ($perc_nonMissing >= $minCount) {
      
      # Check if the minor allele frequency is above the threshold
      $string =~ s/A/AA/g;
      $string =~ s/C/CC/g;
      $string =~ s/G/GG/g;
      $string =~ s/T/TT/g;
      $string =~ s/M/AC/g;
      $string =~ s/R/AG/g;
      $string =~ s/W/AT/g;
      $string =~ s/S/GC/g;
      $string =~ s/Y/CT/g;
      $string =~ s/K/GT/g;
      my @alleles = split(/\//, $info[1]);
      $string =~ s/$alleles[0]//ig;
      my $maf = (length $string)/($numInd * 2);
      if ($maf >= $minMAF) {
	print OUT $_, "\n";
      } else {
	next;
      }
    } else {
      next;
    }
  }
}
close(IN);
close(OUT);
exit;

