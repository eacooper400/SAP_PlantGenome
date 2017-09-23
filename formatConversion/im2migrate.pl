#!/usr/bin/perl -w
#
# im2migrate.pl
#
# Convert IMa2 formatted input file to Migrate formatted input file
#
# October 3, 2013
# Liz Cooper

use strict;

# Get the IMa2 input file and the name of the migrate output file from the command line
my ($USAGE) = "\n$0 <input.ima2> <output.migrate>
\tinput.ima2 = The input file in IMa2 format
\toutput.migrate = The output file to create, in MIGRATE required input format\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output) = @ARGV;

# Save the number of loci, population list, and number of populations outside of the loop
my $num_loci = 0;
my $num_pops = 0;
my @pop_names = ();
my %pop_info = ();
my $header = '';
my @lengths = ();
my @pop_sizes = ();


open (IN, $input) || die "\nUnable to open the file $input!\n";

# Start reading through the first few lines of the IMa2 files
while (<IN>) {
  chomp $_;
  if ($_ =~ /^\#/) {
    $_ =~ s/^\#\s{1,}//;
    $header = $_;
  } elsif ($. == 2) {
    $num_pops = $_;
  } elsif ($. == 3) {
    @pop_names = split(/\s/, $_);
  } elsif ($_ =~ /\(/) {
    next;
  } elsif ($. == 5) {
    $num_loci = $_;
  } elsif ($_ =~ /^Locus/) {
    my @info = split(/\s/, $_);
    my $len = $info[$num_pops+1];
    push (@lengths, $len); 
    @pop_sizes = @info[1..$num_pops];
  } else {
    foreach my $popname (@pop_names) {
      if ($_ =~ /$popname/) {
	if (exists $pop_info{$popname}) {
	  push (@{$pop_info{$popname}}, $_);
	} else {
	  @{$pop_info{$popname}} = ($_);
	}
      }
    }
  }
}
close(IN);

# Print the lines to the Migrate output file
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";
print OUT $num_pops, " ", $num_loci, " ", $header, "\n";
print OUT "@lengths\n";
for (my $p = 0; $p < scalar @pop_names; $p++) {
  my @lines = @{$pop_info{$pop_names[$p]}};
  my $numIn = $pop_sizes[$p];
  print OUT $numIn, "\t", $pop_names[$p], "\n";
  foreach my $line (@lines) {
    print OUT $line, "\n";
  }
}
close(OUT);
exit;
