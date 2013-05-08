#!/usr/bin/perl

use warnings;
use strict;
use AT::Tools::RangeHandler;

my ($FT_FN) = @ARGV;
unless($FT_FN) {
    die "usage: perl $0 <bed-file>\n";
}

# read features to calculate density for
print STDERR "Reading features...\n";
my %fts_by_chr;
open IN, $FT_FN or die "could not open $FT_FN";
while(my $line = <IN>) {
    if($line =~ /^browser/ or $line =~ /^track/) {
	print $line;
	next;
    }
    chomp $line;
    my ($chr, $start, $end) = split /\t/, $line;
    my $strand = '.';  #$strand = '.' unless($SEP_STRANDS);
    push @{$fts_by_chr{$chr}}, [$start+1,$end];
}
close IN;

# process each chromosome
foreach my $chr (sort keys %fts_by_chr) {
    print STDERR "Processing $chr...\n";
    my $fts = $fts_by_chr{$chr};
    my $ft_ivs = calc_expressed_intervals($fts);
    foreach my $iv (@$ft_ivs) {
	print join("\t", $chr, $iv->[0]-1, $iv->[1]), "\n";
    }
}

print STDERR "Done!\n";


sub calc_expressed_intervals
{
    my ($exons_ref) = @_;

    return [] unless(@$exons_ref);

    my @hsps = sort {$a->[0] <=> $b->[0]} @$exons_ref;

    my @iv;
    my ($start, $end) = @{$hsps[0]};
    for my $i (1..@hsps-1) {
	my ($my_start, $my_end) = @{$hsps[$i]};
	if($my_start > $end) {
	    push @iv, [$start,$end];
	    ($start, $end) = ($my_start, $my_end); 
	}
	else {
	    $end = $my_end if ($end < $my_end);
	}
    }
    push @iv, [$start,$end];
    print STDERR "collapsed ", scalar(@$exons_ref)," into ",scalar(@iv)," intervals\n";
    return \@iv;
}
