#!/usr/bin/perl -w
use warnings;
use strict;

while(my $line  = <STDIN>) {

    # Parse line
    chomp $line;
    my ($gene_id, $type, $symbol, undef, undef, $chr, $gene_start, $gene_end, $strand) = split /\t/, $line;

    # Check that we have all fields we need
    next if($type ne 'Gene' or $symbol eq 'null' or $chr eq 'null' or
	    $gene_start eq 'null' or $gene_end eq 'null' or $strand eq 'null');

    # Change chr name to UCSC-style and check that it is valid
    $chr = 'M' if($chr eq 'MT');
#    unless($valid_chr{"chr$chr"}) {
#	$skipped_chr{$chr} = 1;
#	next;
#    }
    $chr = "chr$chr";

    print join("\t", 
	       $chr,
	       $gene_start-1,
	       $gene_end,
	       $gene_id,
	       0,
	       $strand), "\n";

    
}
