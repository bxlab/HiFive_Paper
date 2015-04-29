#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use Data::Dumper;
use List::Util qw(first);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <input file> <output file> <output bin index file> <binsize>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $ofn = $ARGV[1];
my $ofn_bins = $ARGV[2];
my $binsize = $ARGV[3];

# hash table with all fends
my %fends;
my $index = 1;
my %bins;

##########################################################################################
# traverse fends file
##########################################################################################

open(IN, $ifn) || die;
print STDERR "Reading input file $ifn into hash...\n";
my $header = <IN>;
chomp($header);
my %h = parse_header($header);

while (my $line = <IN>) {
	chomp $line;
	my @f = split("\t", $line);

	my $fend = $f[$h{fend}];
	my $chr = $f[$h{chr}];
	my $coord = $f[$h{coord}];

	# get bin
	my $bcoord = int($coord / $binsize) * $binsize;

	# add to hash
	$fends{$chr} = {} if !defined($fends{$chr});
	if (!defined($fends{$chr}->{$bcoord})) {
		$fends{$chr}->{$bcoord} = {};
		$fends{$chr}->{$bcoord}->{index} = $index;
		$fends{$chr}->{$bcoord}->{lines} = {};
		$bins{$index} = {};
		$bins{$index}->{chr} = $chr;
		$bins{$index}->{from} = $bcoord;
		$bins{$index}->{to} = $bcoord + $binsize;
		$bins{$index}->{count} = 0;
		$index++;
	}
	# update bin counter
	$bins{$fends{$chr}->{$bcoord}->{index}}->{count}++;

	# add line
	$fends{$chr}->{$bcoord}->{lines}->{$fend} = $line;
};
close(IN);

open(OUT, ">", $ofn) || die;
print OUT $header, "\tcoord_bin\n";
print STDERR "generating $ofn\n";

foreach my $chr (sort keys %fends) {
	foreach my $coord (sort {$a <=> $b} keys %{$fends{$chr}}) {
		foreach my $fend (sort {$a <=> $b} keys %{$fends{$chr}->{$coord}->{lines}}) {
			print OUT $fends{$chr}->{$coord}->{lines}->{$fend}, "\t", $fends{$chr}->{$coord}->{index}, "\n";
		}
	}
}
close(OUT);


open(OUT, ">", $ofn_bins) || die;
print OUT "cbin\tchr\tfrom.coord\tto.coord\tcount\n";
print STDERR "generating $ofn_bins\n";
foreach my $bin (sort {$a <=> $b} keys %bins) {
	print OUT $bin, "\t", $bins{$bin}->{chr}, "\t", 
	          $bins{$bin}->{from}, "\t", $bins{$bin}->{to}, "\t", $bins{$bin}->{count}, "\n";
}
close(OUT);

##########################################################################################

######################################################################################################
# Subroutines
######################################################################################################

sub parse_header
{
	my ($header) = @_;
	chomp($header);
	my @f = split("\t", $header);
	my %result;
	for (my $i = 0; $i <= $#f; $i++) {
		$result{$f[$i]} = $i;
	}
	return %result;
}
