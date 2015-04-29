#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use Data::Dumper;
use List::Util qw(first);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <input prefix> <output prefix> <min segment length> <max segment length> <rejected file> <pair file1, pair file2, ...>\n";
	print STDERR "   fend output indices start from one\n";
	print STDERR "   set segment max length to 0 to disable\n";
	exit 1;
}

my $in_fends_fn = $ARGV[0];
my $oprefix = $ARGV[1];
my $min_len = $ARGV[2];
my $max_len = $ARGV[3];
my $rejected_fn = $ARGV[4];

shift; shift; shift; shift; shift;
my @pair_files = @ARGV;

print "reading paired read files:\n", join("\n", @pair_files), "\n";

my $out_stats_fn = $oprefix.".mat_stats";
my $out_matrix_fn = $oprefix.".mat";

# remove stats file
system("rm -rf $out_stats_fn") == 0 || die;
print STDERR "log file: $out_stats_fn\n";

# hash table with all fends
our %fends;
our %fend_matrix;

# to quickly get from approximate coord to fend
our %coord2index;

##########################################################################################
# read fends file
##########################################################################################

open(IN, $in_fends_fn) || die "cannot open $in_fends_fn";
print STDERR "Reading input file $in_fends_fn into hash...\n";
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
	chomp $line;
	my @f = split("\t", $line);
	my $fend = $f[$h{fend}];
	my $chr = $f[$h{chr}];
	my $coord = $f[$h{coord}];

	!defined($fends{$fend}) or die "non-unique fend";
	$fends{$fend} = {};
	$fends{$fend}->{fend} = $fend;
	$fends{$fend}->{frag} = $f[$h{frag}];
	$fends{$fend}->{chr} = $chr;
	$fends{$fend}->{coord} = $coord;
	$fends{$fend}->{frag_len} = $f[$h{frag_len}];

	if (!defined($coord2index{$chr}))
	{
		$coord2index{$chr} = {};
		$coord2index{$chr}->{coords} = {};
	}
	$coord2index{$chr}->{coords}->{$coord} = $fend;
}
close(IN);

#########################################################################################
# compute sorted coords per chrom
##########################################################################################

for my $chr (%coord2index)
{
	my @sorted = sort {$a <=> $b} keys %{$coord2index{$chr}->{coords}};
	$coord2index{$chr}->{sorted_coords} = \@sorted;
}

##########################################################################################
# parse pair files
##########################################################################################


open(REJECT, ">", $rejected_fn) || die;
print REJECT "chr1\tcoord1\tdist1\tchr2\tcoord2\tdist2\n";

my $read_count = 0;
my %stats;

foreach my $fn (@pair_files)
{
	read_pair_file($fn);
}
close(REJECT);

sub read_pair_file
{
	my ($fn) = @_;

	my $appr_lines = apprx_lines($fn);
	print STDERR "traversing file $fn, with about ".int($appr_lines/1000000)."M lines\n";


	open(IN, $fn) || die;

	# !!! Chan files do not have a header
	# <IN>;

	while (my $line = <IN>) 
	{
		$read_count++;
		print STDERR "line: $read_count\n" if ($read_count % 1000000 == 0);

		chomp $line;
		my @f = split("\t", $line);
		my $chr1 = $f[0];
		my $coord1 = $f[1];
		my $strand1 = $f[2];
		my $chr2 = $f[3];
		my $coord2 = $f[4];
		my $strand2 = $f[5];

		my ($fend1, $fend_coord1, $index1) = find_closest_fend($chr1, $coord1, $strand1);
		my ($fend2, $fend_coord2, $index2) = find_closest_fend($chr2, $coord2, $strand2);

		if ($fend1 == -1 or $fend2 == -1)
		{
			$stats{not_found}++;
			next;
		}

		my $dist1 = abs($fend_coord1 - $coord1);
		my $dist2 = abs($fend_coord2 - $coord2);
		my $dist = $dist1 + $dist2;

		# check length threshold
		if ($dist < $min_len || ($max_len > 0 && $dist > $max_len))
		{
			$stats{len_outofrange}++;
			print REJECT "$chr1\t$coord1\t$dist1\t$chr2\t$coord2\t$dist2\n";
			next;
		}
		my $frag1 = $fends{$fend1}->{frag};
		my $frag2 = $fends{$fend2}->{frag};

		my $not_ligated = ($chr1 eq $chr2) && 
			              (between($coord1, $fend_coord1, $coord2) ||
						   between($coord2, $fend_coord2, $coord1));

		# site got broken and got fused back
		my $re_ligated = ($chr1 eq $chr2) && (abs($frag1-$frag2) == 1) && (abs($index1-$index2) == 1);

		# two ends of one fragment got fused
		my $self_ligated = ($chr1 eq $chr2) && ($frag1 eq $frag2);

		# no ligation
		if ($not_ligated)
		{
			$stats{no_ligation}++;
			($frag1 == $frag2) or die "fend1=$fend1, fend2=$fend2. $frag1 != $frag2, line:$line\n";
			(abs($index1-$index2) <= 1) or die;
			next;
		}

		if ($re_ligated)
		{
			$stats{re_ligation}++;
		} elsif ($self_ligated)
		{
			$stats{self_ligation}++;
		} else {
			$stats{mapped}++;
		}

		my $fend_small = ($fend1 < $fend2) ? $fend1 : $fend2;
		my $fend_large = ($fend1 < $fend2) ? $fend2 : $fend1;
		$fend_matrix{$fend_small} = {} if !defined($fend_matrix{$fend_small});
		$fend_matrix{$fend_small}->{$fend_large} = 0 if !defined($fend_matrix{$fend_small}->{$fend_large});
		$fend_matrix{$fend_small}->{$fend_large}++;
	}
	close(IN);

}

######################################################################################################
# write mat file
######################################################################################################

my %fend_stats;
my $fend_count = 0;
open(OUT, ">", $out_matrix_fn) || die;
print OUT "fend1\tfend2\tcount\n";
foreach my $fend1 (sort { $a <=> $b } keys %fend_matrix)
{
	foreach my $fend2 (sort { $a <=> $b } keys %{$fend_matrix{$fend1}})
	{
		$fend_count++;
		my $count = $fend_matrix{$fend1}->{$fend2};
		print OUT $fend1, "\t" , $fend2, "\t" , $count, "\n";

		# fend to itself (small chance in single cell)
		if ($fend1 == $fend2) {
			$fend_stats{same_fend}++;
			next;
		}

		# trans
		if ($fends{$fend1}->{chr} ne $fends{$fend2}->{chr}) {
			$fend_stats{trans}++;
			next;
		}

		# cis
		if ($fends{$fend1}->{frag} == $fends{$fend2}->{frag}) {
			$fend_stats{self_ligation}++;
		} elsif (abs($fends{$fend1}->{coord} - $fends{$fend2}->{coord}) < 10) {
			$fend_stats{re_ligation}++;
		} elsif (abs($fends{$fend1}->{coord} - $fends{$fend2}->{coord}) < 10000) {
			$fend_stats{close_cis}++;
		} else {
			$fend_stats{far_cis}++;
		}

	}
}
close(OUT);

#open(OUT, ">", $out_stats_fn) || die;
#print OUT "Total reads: $read_count\n";
#print OUT 
#	"no site found (probably chrom end): ", perc_str($stats{not_found}, $read_count), "\n",
#	"segment length out of range: ", perc_str($stats{len_outofrange}, $read_count), "\n",
#	"no ligation: ", perc_str($stats{no_ligation}, $read_count), "\n",
#	"self ligation: ", perc_str($stats{self_ligation}, $read_count), "\n",
#	"re ligation: ", perc_str($stats{re_ligation}, $read_count), "\n";
#print OUT "normal reads: ", perc_str($stats{mapped}, $read_count), "\n"; 
#print OUT "------------------------------------------\n";
#print OUT "unique fend pairs: $fend_count\n";
#print OUT "fend to itself: ", perc_str($fend_stats{same_fend}, $fend_count), "\n";
#print OUT "frag self ligation: ", perc_str($fend_stats{self_ligation}, $fend_count), "\n";
#print OUT "re ligation: ", perc_str($fend_stats{re_ligation}, $fend_count), "\n";
#print OUT "close cis pairs (<10k): ", perc_str($fend_stats{close_cis}, $fend_count), "\n";
#print OUT "far cis pairs (>10k): ", perc_str($fend_stats{far_cis}, $fend_count), "\n";
#print OUT "trans pairs: ", perc_str($fend_stats{trans}, $fend_count), "\n";
#print OUT "-------------------------------------\n";
#close(OUT);

######################################################################################################
# Subroutines
######################################################################################################

# check if C is between and A,B
sub between
{
	my ($A, $B, $C) = @_;
	$B = $B - $A;
	$C = $C - $A;
	return (($B>$C && $C>0) || ($B<$C && $C<0));
}

# coord should be rounded to binsize
sub add_fend
{
	my ($chr, $bin, $fend, $strand) = @_;
	$coord2index{$chr} = {} if !defined($coord2index{$chr});
	$coord2index{$chr}->{$bin} = {} if !defined($coord2index{$chr}->{$bin});

	# mark bin if multiple fends map to it
	if (!defined($coord2index{$chr}->{$bin}->{$strand}))
	{	
		$coord2index{$chr}->{$bin}->{$strand} = $fend;
	}
	else 
	{
		$coord2index{$chr}->{$bin}->{$strand} = -1;
	}
}

sub apprx_lines
{
	my ($fn) = @_;
	my $tmp = "/tmp/".$$."_apprx_lines.tmp";
	system("head -n 1000 $fn > $tmp");
	my $size_head = -s $tmp;
	my $size_all = -s $fn;
	$size_head > 0 or die;
	return (int($size_all/$size_head*1000));
}

sub perc_str
{
	my ($n, $total) = @_;
    return ($n." (".(int(1000 * $n / $total)/10)."%)");
}

sub find_closest_fend
{
	my ($chr, $coord, $strand) = @_;
	return ((-1,-1)) if (!defined($coord2index{$chr}));
	my $index = ($strand eq "+") ?
		binary_search($coord2index{$chr}->{sorted_coords}, $coord, 1) :
		binary_search($coord2index{$chr}->{sorted_coords}, $coord, 0);
	return ((-1,-1)) if ($index == -1);
	my $fend_coord = $coord2index{$chr}->{sorted_coords}[$index];
	defined ($coord2index{$chr}->{coords}->{$fend_coord}) and defined($fend_coord) or die;
	return (($coord2index{$chr}->{coords}->{$fend_coord}, $fend_coord, $index));
}

# returns first element above/below value in sorted array
sub binary_search 
{
	my $arr = shift;
	my $value = shift;
	my $above = shift; 

	my $left = 0;
	my $right = $#$arr;
	
	while ($left <= $right) {
		my $mid = ($right + $left) >> 1;
		my $c = $arr->[$mid] <=> $value;
		return $mid if ($c == 0);
		if ($c > 0) {
			$right = $mid - 1;
		} else {
			$left  = $mid + 1;
		}
	}
	$left = -1 if ($left > $#$arr);
	$right = -1 if ($right < 0);
	return ($above ? $left : $right);
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
