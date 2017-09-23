#!/usr/bin/perl -w
#
# plink2hmp_numeric.pl
#
# Convert plink format (ped and map) files to a numeric hapmap format file
#
# December 9, 2015
# Liz Cooper

use strict;

# Take as input the name of the ped, map, and output hmp files
my ($USAGE) = "\n$0 <input.ped> <input.map> <output.hmp.txt>
\tinput.ped = The input .ped file from Plink
\tinput.map = The input .map file from Plink
\toutput.hmp.txt = The name of the hapmap file to create\n\n";
unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($pedfile, $mapfile, $output) = @ARGV;

# Read through the ped file first, and save 2 hashes
# First, for each position, save the alleles
# Next, for each individual, save the genotype string
my %alleles = ();
my %genotypes = ();

open (PED, $pedfile) || die "\nUnable to open the file $pedfile!\n";
while (<PED>) {
  chomp $_;
  my @info = split(/\s{1,}/, $_);
  my $name = $info[1];
  my @markers = @info[6..(scalar @info - 1)];
  my $pair = '';
  for (my $m = 0; $m < scalar @markers; $m++) {
    if (0 == $m % 2) {
      unless (length $pair == 0) {
	my $gen = genotype_lookup($pair);
	$genotypes{$name} .= $gen;
      }
      $pair = '';
      $pair .= $markers[$m];
    } else {
      $pair .= $markers[$m];
    }
  }
}
close(PED);

# Now get the alleles for every possible position
my @genstrings = values %genotypes;
my $first = $genstrings[0];
for (my $p = 0; $p < length $first; $p++) {
  my %temp = ();
  foreach my $string (@genstrings) {
    my $gcode = substr($string, $p, 1);
    my @acodes = ();
    allele_lookup(\$gcode, \@acodes);
    foreach my $a (@acodes) {
      unless ($a =~ /N/) {
	if (exists $temp{$a}) {
	  $temp{$a} += 1;
	} else {
	  $temp{$a} = 1;
	}
      }
    }
  }
  @{$alleles{$p}} = (keys %temp);
}

# Open the output file for printing
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

# Print a header line to the output file
print OUT "rs\t", "alleles\t", "chrom\t", "pos\t", "strand\t", "assembly\t", "center\t", "protLSID\t", "assayLSID\t", "panelLSID\t", "Qcode\t";
my @names = keys %genotypes;
my $namestring = join("\t", @names);
print OUT $namestring, "\n";

# Open the map file, then read and process one position at a time
# Filter out positions with more or fewer than 2 alleles
open (MAP, $mapfile) || die "\nUnable to open the file $mapfile!\n";
my $linecounter = 0;

while (<MAP>) {
  chomp $_;
  my ($chrom, $rs, $sex, $pos) = split(/\s{1,}/, $_);
  $rs =~ s/chr/S/ig;
  $rs =~ s/\:/_/ig;
  my @pos_alleles = @{$alleles{$linecounter}};
  unless ((scalar @pos_alleles < 2) || (scalar @pos_alleles > 2)) {
    print OUT $rs, "\t";
    my $allelecol = join("/", @pos_alleles);
    print OUT $allelecol, "\t", $chrom, "\t", $pos, "\t", "NA\t", "NA\t", "NA\t", "NA\t", "NA\t", "NA\t", "NA"; 
    foreach my $name (@names) {
      my $string1 = substr($genotypes{$name}, $linecounter, 1);
      $string1 =~ s/$pos_alleles[0]/0/ig;
      $string1 =~ s/$pos_alleles[1]/2/ig;
      $string1 =~ s/N/9/ig;
      $string1 =~ s/[MRWSYK]/1/ig;
      print OUT "\t", $string1;
    }
    print OUT "\n";
  }
  $linecounter++;
}
close(MAP);
close(OUT);
exit;
###################################################
# Subroutines
###################################################
sub genotype_lookup {
  my ($string) = @_;
  my %codes = (
	       'AA' => 'A',
	       'AC' => 'M',
	       'AG' => 'R',
	       'AT' => 'W',
	       'AN' => 'A',
	       'CA' => 'M',
	       'CC' => 'C',
	       'CG' => 'S',
	       'CT' => 'Y',
	       'CN' => 'C',
	       'GA' => 'R',
	       'GC' => 'S',
	       'GG' => 'G',
	       'GT' => 'K',
	       'GN' => 'G',
	       'TA' => 'W',
	       'TC' => 'Y',
	       'TG' => 'K',
	       'TT' => 'T',
	       'TN' => 'T',
	       'NA' => 'A',
	       'NC' => 'C',
	       'NG' => 'G',
	       'NT' => 'T',
	       'NN' => 'N',
	       '00' => '0',
	       '11' => '1',
	       '01' => 'N',
	       '10' => 'N');

  my $single = $codes{$string};
  return($single);
}

####################################################
sub allele_lookup {
  my ($string, $array) = @_;
  my %codes = (
	       'A' => 'A,A',
	       'C' => 'C,C',
	       'G' => 'G,G',
	       'T' => 'T,T',
	       'M' => 'A,C',
	       'R' => 'A,G',
	       'W' => 'A,T',
	       'S' => 'C,G',
	       'Y' => 'C,T',
	       'K' => 'G,T',
	       'N' => 'N,N',
	       '0' => '0,0',
	       '1' => '1,1');
  my $list = $codes{$$string};
  @$array = split(/,/, $list);
}
