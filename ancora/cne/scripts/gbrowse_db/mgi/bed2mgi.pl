#!/usr/bin/perl -w
use warnings;
use strict;

my ($BED_FN, $MGI_FN) = @ARGV; 
unless($BED_FN and $MGI_FN) {
    die "Usage: perl $0 <bed-in-file> <mgi-in-file>\nOutput on stdout.\n";
}

# Read locations from bed file and index by id
my %id_to_loc;
open BED, $BED_FN or die "could not open file $BED_FN\n";
while(my $line  = <BED>) {
    chomp $line;
    my ($chr, $start, $end, $id, undef, $strand) = split /\t/, $line;
    $id_to_loc{$id} = [$chr,$start+1,$end,$strand];
}
close BED;

# Process MGI file
open MGI, $MGI_FN or die "could not open file $MGI_FN\n";
# Read and print header
my $header = <MGI>;
die "Header not found in MGI file" unless($header =~ /^MGI accession id/);
print $header;
# Process data in file
while(my $line  = <MGI>) {

    # Parse line
    chomp $line;
    my ($gene_id, $type, $symbol, $name, $other_id) = split /\t/, $line;

    # Check that we have all fields we need
    next if($type ne 'Gene' or $symbol eq 'null');

    # Lookup location in bed file
    my $loc = $id_to_loc{$gene_id};
    next unless($loc);
    my ($chr, $start, $end, $strand) = @$loc;

    $chr =~ s/^chr//;
    $chr = 'MT' if($chr eq 'M');
    print join("\t", $gene_id, $type, $symbol, $name, $other_id, $chr, $start, $end, $strand), "\n";
    
}
close MGI;
