#!/usr/bin/perl -w
#
# SNP_density.pl
#
# Find the number of SNPs in windows of pre-specified size, and return values for genome-wide plots
#
# June 25, 2014
# Edited: September 18, 2014
# Liz Cooper

use strict;

# Take as input arguments the:
# Name of the hapmap file with the SNPs
# Name of the output file with windows and numbers of SNPs
# Window size for counting SNPs (in base pairs)
# Amount of overlap between windows (in base pairs)

my ($USAGE) = "$0 <input.hmp.txt> <output.txt> <window_size> <overlap>
\t\tinput.hmp.txt = A SNP file in tab-delimited format (like the hapmap format used by TASSEL)
\t\toutput.txt = The name of the output file to create
\t\twindow_size = The size of the window (in base pairs)
\t\toverlap = The amount of overlap between windows (in base pairs)\n";
unless (@ARGV) {
  print $USAGE;
  exit;
}

my ($input, $output, $winsize, $overlap) = @ARGV;

# For each window, save the count of SNPs and the window positions outside of the loop
# Also, save the count of SNPs in the region of overlap between the current window and the next
my $win_start = 0;
my $win_end = $winsize;
my $snp_count = 0;
my $overlap_start = $win_end - $overlap;
my $overlap_snp_count = 0;

my $current_chr = 1;

# Open the output file, to print the results from each window
open (OUT, ">$output") || die"\nUnable to open the file $output!\n";

# Open the input file, and start reading through one line at a time
open (IN, $input) || die "\nUnable to open the file $input!\n";

while (<IN>) {
  chomp $_;
  
  # Skip the header line
  if ($. < 2) {
    next;
  }

  # Split the line into an array, then check if the chromosome and position are within the same window
  my @info = split(/\s{1,}/, $_);
  $info[2] =~ s/Chr0//g;
  $info[2] =~ s/Chr//g;
  $info[2] =~ s/super_/1/g;
  if (($info[2] == $current_chr) && ($info[3] <= $win_end)) {
        $snp_count++;

    # Check if this position also falls into the region of overlap with the next window
	if ($info[3] >= $overlap_start) {
	  $overlap_snp_count++;
	}

  } elsif (($info[2] == $current_chr) && ($info[3] > $win_end)) {
        
    # If the chromosome is the same and the position belongs in the next window,
    # then print out the information for the current window, and re-set the position counts for the new window
    print OUT $current_chr, "\t", $win_start, "\t", $win_end, "\t", $snp_count, "\n";

    my $flag = 0;
    do {
      $win_start = $overlap_start;
      $win_end = $win_start + $winsize;
      $snp_count = $overlap_snp_count;
      $current_chr = $info[2];
      $overlap_start = $win_end - $overlap;
      $overlap_snp_count = 0;
      if ($info[3] > $win_end) {
	 print OUT $current_chr, "\t", $win_start, "\t", $win_end, "\t", $snp_count, "\n";
	 $flag = 0;
       } else {
	 $snp_count += 1;
	 $flag = 1;
       }
    } until ($flag > 0);

  } elsif ($info[2] > $current_chr) {

    # If this is the next chromosome, 
    # then print the information for the current window, and re-set all of the positions and counts
    print OUT $current_chr, "\t", $win_start, "\t", $win_end, "\t", $snp_count, "\n";
    $win_start = 0;
    $win_end = $winsize;
    $snp_count = 1;
    $current_chr = $info[2];
    $overlap_start = $win_end - $overlap;
    $overlap_snp_count = 0;
  }
}
close (IN);

# Print the information from the last window, and then close the output file
print OUT $current_chr, "\t", $win_start, "\t", $win_end, "\t", $snp_count, "\n";

close(OUT);
exit;
