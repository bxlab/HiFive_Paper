#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <input prefix> <max distance> <output file>\n";
	exit 1;
}

my $prefix = $ARGV[0];
my $max_distance = $ARGV[1];
my $out_fn = $ARGV[2];

my $site_fn = $prefix.".cbinned";
my $cbins_fn = $prefix.".cbins";
my $mat_fn = $prefix.".mat";

if ($max_distance > 0) {
	print STDERR "Reporting cis up to ".$max_distance."bp\n";
} else {
	print STDERR "Reporting all contacts\n";
}

#############################################################################################
# cbin file
#############################################################################################


our %cbins;

open(IN, $cbins_fn) || die;
print STDERR "Reading file $cbins_fn into hash...\n";
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
	chomp $line;
    # cbin    chr     from.coord      to.coord
	my @f = split("\t", $line);
	my $cbin = $f[$h{cbin}];

	$cbins{$cbin} = {};
	$cbins{$cbin}->{chr} = $f[$h{chr}];
	$cbins{$cbin}->{from} = $f[$h{"from.coord"}];
	$cbins{$cbin}->{to} = $f[$h{"to.coord"}];
}
close(IN);

#############################################################################################
# fend file
#############################################################################################

our %fends;

open(IN, $site_fn) || die;
print STDERR "Reading file $site_fn into hash...\n";
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) {
	chomp $line;
    # fend    frag    strand  chr     coord   frag_len_bin    fragend_len_bin coord_bin
	my @f = split("\t", $line);
	my $fend = $f[$h{fend}];

	!defined($fends{$fend}) or die "non-unique fend";
	$fends{$fend} = {};
	$fends{$fend}->{fend} = $fend;
	$fends{$fend}->{chr} = $f[$h{chr}];
	$fends{$fend}->{coord} = $f[$h{coord}];
	$fends{$fend}->{cbin} = $f[$h{coord_bin}];
}
close(IN);

#############################################################################################
# mat file
#############################################################################################

my $appr_lines = apprx_lines($mat_fn);
print STDERR "traversing file $mat_fn, with about ".int($appr_lines/1000000)."M lines\n";

my %bin_counters;

my $count = 0;
open(IN, $mat_fn) || die;
$header = <IN>;
%h = parse_header($header);
while (my $line = <IN>) 
{
	$count++;
	print STDERR "line: $count\n" if ($count % 1000000 == 0);

	chomp $line;
	my @f = split("\t", $line);
	# fend1 chr1 coord1 fend2 chr2 coord2 count
	my $fend1 = $f[$h{fend1}];
	my $fend2 = $f[$h{fend2}];
	my $count = $f[$h{count}];

	next if (!defined ($fends{$fend1}) or !defined ($fends{$fend2}));

	my $coord1 = $fends{$fend1}->{coord};
	my $coord2 = $fends{$fend2}->{coord};
	my $chr1 = $fends{$fend1}->{chr};
	my $chr2 = $fends{$fend2}->{chr};

	# limit to of max_distance > 0
	next if ($max_distance > 0 && (($chr1 ne $chr2) || (abs($coord1 - $coord2) > $max_distance)));

	my $cbin1 = $fends{$fend1}->{cbin};
	my $cbin2 = $fends{$fend2}->{cbin};

	my $bin = make_bin($cbin1, $cbin2);
	if (!defined($bin_counters{$bin}))
	{
		$bin_counters{$bin} = {};
		$bin_counters{$bin}->{count} = 0; # real count
		$bin_counters{$bin}->{unique_count} = 0; # count each pair once
	}

	$bin_counters{$bin}->{count} += $count;
	$bin_counters{$bin}->{unique_count} += 1;
}
close(IN);

print STDERR "Writing output file: $out_fn\n";

open(OUT, ">", $out_fn) || die;
#print OUT "coord_bin1\tchr1\tfrom1\tto1\tcoord_bin2\tchr2\tfrom2\tto2\tobserved_count\tunique_observed_count\n";
print OUT "coord_bin1\tchr1\tfrom1\tto1\tcoord_bin2\tchr2\tfrom2\tto2\tobserved_count\n";
foreach my $bin (sort sort_bins keys %bin_counters)
{
	my ($cbin1, $cbin2) = split(",", $bin);
	my $ucount = $bin_counters{$bin}->{unique_count};

	print OUT $cbin1, "\t", $cbins{$cbin1}->{chr}, "\t", $cbins{$cbin1}->{from}, "\t", $cbins{$cbin1}->{to}, "\t";
	print OUT $cbin2, "\t", $cbins{$cbin2}->{chr}, "\t", $cbins{$cbin2}->{from}, "\t", $cbins{$cbin2}->{to}, "\t";
	print OUT $ucount, "\n";
}
close(OUT);

######################################################################################################
# Subroutines
######################################################################################################

sub sort_bins
{
	my @as = split(",", $a);
	my @bs = split(",", $b);
	return ( ($as[0] <=> $bs[0]) or ($as[1] <=> $bs[1]) );
}

sub make_bin
{
	my ($bin1, $bin2) = @_;
	return join(",", sort {$a <=> $b} $bin1, $bin2);
}

sub apprx_lines
{
	my ($fn) = @_;
	my $tmp = "/tmp/".$$."_apprx_lines.tmp";
	system("head -n 100000 $fn > $tmp");
	my $size_head = -s $tmp;
	my $size_all = -s $fn;
	return (int($size_all/$size_head*100000));
}

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
