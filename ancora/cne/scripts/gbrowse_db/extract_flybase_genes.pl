#!/usr/bin/perl -w
use strict;
use warnings;
use Class::Struct;
use AT::DB::Binner;
use AT::Tools::RangeHandler;

my $SOURCE = "FlyBase";
my ($FLYBASE_FN) = @ARGV;

unless($FLYBASE_FN) {

    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Usage: 
    
perl $cmd <flybase file>

The input file should be a GFF annotation file from FlyBase.
For FlyBase Release 2007_1 the file is:
ftp://ftp.flybase.net/releases/FB2007_01/dmel_r5.2_FB2007_01/gff/dmel-all-r5.2.gff.gz

This script will parse the FlyBase GFF file and create the following files:
flyBaseGene.txt
flyBaseNoncoding.txt
flyBase2004Xref.txt
flybase-gbrowse.gff

The first three follow the UCSC database schema for FlyBase genes,
while the last is a GFF file that can be imported into a GBrowse Bio::DB::GFF database.

WARNING: these files will be overwritten if they exist!

EndOfUsage

    exit;
}
    

struct(transcript => [gene_id => '$',
		      type => '$',
		      transcript_symbol => '$',
		      chr => '$',
		      start => '$',
		      end => '$',
		      strand => '$',
		      cds_start => '$',
		      cds_end => '$',
		      exons => '@']);

struct(gene => [gene_symbol => '$',
		chr => '$',
		start => '$',
		end => '$',
		strand => '$']);

my $BINNER = AT::DB::Binner->new();

my $transcripts = {};
my $genes = {};
my %gene_id_to_symbol;
my %transcript_symbols;

open IN, $FLYBASE_FN or die "could not open $FLYBASE_FN";

print STDERR "Reading $FLYBASE_FN...\n";

while(my $line = <IN>) {
    next if($line =~ /^\#/);
    chomp $line;
    my ($chr, $source, $type, $start, $end, undef, $strand, undef, $extras) = 
	split /\t/, $line;
    next unless(defined($source) and $source eq 'FlyBase');
    if($type eq 'gene') {
	my ($gene_id) = $extras =~ /ID=(\w+);/;
	my ($gene_symbol) = $extras =~ /Name=([^;]+);/;
	die "duplicate gene " if($genes->{$gene_id});
	my $gene = $genes->{$gene_id} = gene->new();
	$gene->chr($chr);
	$gene->start($start);
	$gene->end($end);
	$gene->strand($strand);
	$gene->gene_symbol($gene_symbol);
    }
    elsif($type eq 'CDS') {
	foreach my $transcript_id (get_parents($extras)) {
	    my $transcript = get_transcript($transcripts, $transcript_id);
	    $transcript->cds_start($start) if(!defined($transcript->cds_start) or $start < $transcript->cds_start);
	    $transcript->cds_end($end) if(!defined($transcript->cds_end) or $end > $transcript->cds_end);
	}
    }
    elsif($type eq 'exon') {
	foreach my $transcript_id (get_parents($extras)) {
	    my $transcript = get_transcript($transcripts, $transcript_id);
	    push @{$transcript->exons}, [$start, $end];
	}
    }
    elsif($type =~ /RNA$/ or $type eq 'pseudogene') {
	my ($transcript_id) = $extras =~ /ID=(\w+);/;
	my ($transcript_symbol) = $extras =~ /Name=([^;]+);/;
	my @parents = get_parents($extras);
	my $transcript = get_transcript($transcripts, $transcript_id);
	die "duplicate transcript $transcript_id" if($transcript->start);
	die "duplicate transcript symbol $transcript_symbol" if($transcript_symbols{$transcript_symbol});
	die "multiple genes for transcript $transcript_id" if(@parents > 1);
	$transcript_symbols{$transcript_symbol} = 1;
	$transcript->type($type);
	$transcript->chr($chr);
	$transcript->start($start);
	$transcript->end($end);
	$transcript->strand($strand);
	$transcript->transcript_symbol($transcript_symbol);
	$transcript->gene_id($parents[0]);
    }

# noncoding types: tRNA, ncRNA, rRNA, miRNA, snoRNA, snRNA, snRNA
# pseudogene: pseudogene
}

close IN;

validate_transcripts($transcripts);

if(%$transcripts) {
    write_ucsc_tables($transcripts, $genes);
    write_gbrowse_gff($transcripts, $genes, $SOURCE);
}
else {
    print STDERR "No transcript data found in the input files.\nNo files created.\n";
}

print STDERR "Done.\n";
exit;


sub validate_transcripts
{
    my ($transcripts) = @_;
    while(my ($transcript_id, $t) = each %$transcripts) {
	if($t->cds_start) {
	    warn "Warning: transcript $transcript_id has CDS coordinates but is not classified as mRNA.\n" if($t->type ne 'mRNA');
	}
	else {
	    warn "Warning: mRNA $transcript_id does not have CDS coordinates.\n" if($t->type eq 'mRNA');
	}
	my @exons = sort { $a->[0] <=> $b->[0] } @{$t->exons};
	if($exons[0][0] != $t->start or $exons[-1][1] != $t->end) {
	    warn "Warning: exon and transcript boundaries do not match for $transcript_id.\n";
	    if($exons[0][0] >= $t->start and $exons[-1][1] <= $t->end and $exons[0][0] <= $exons[-1][1]) {
		$t->start($exons[0][0]);
		$t->end($exons[-1][1]);
		warn "  Transcript truncated to exon coordinates.\n";
	    }
	    else {
		warn "  Exon coordinates are inconsistent or exceed transcript boundaries.\n";
	    }
	}
    }
}


sub write_ucsc_tables
{
    my ($transcripts, $genes) = @_;

    print STDERR "Creating files: flyBaseGene.txt, flyBaseNoncoding.txt and flyBase2004Xref.txt...\n";
    open GENE, ">flyBaseGene.txt";
    open NC, ">flyBaseNoncoding.txt";
    open XREF, ">flyBase2004Xref.txt";
    while(my ($transcript_id, $t) = each %$transcripts) {

	# check that we have all data
	my $gene_id = $t->gene_id or die "no gene_id for $transcript_id";
	my $transcript_symbol = $t->transcript_symbol or die "no symbol for $transcript_id";
	my $type = $t->type or die "no type for $transcript_id";
	my $exons = $t->exons or die "no exons for $transcript_id";
	my $chr = $t->chr or die "no chr for $transcript_id";
	my $start = $t->start or die "no start for $transcript_id";
	my $end = $t->end or die "no end for $transcript_id";
	my $strand = $t->strand or die "no strand for $transcript_id";
	my $gene_symbol = $genes->{$t->gene_id}->gene_symbol or die "no symbol for $gene_id";
	$start--;
	
	if($t->cds_start) {
	    # get cds start, end
	    my $cds_start = $t->cds_start;
	    my $cds_end = $t->cds_end or die "no cds end for $transcript_id";
	    $cds_start--;
	    
	    # sort exons and check start, end
	    my (@exon_starts, @exon_ends);
	    foreach my $e (@$exons) {
		push @exon_starts, $e->[0] - 1;
		push @exon_ends, $e->[1];
	    }
	    @exon_starts = sort {$a <=> $b} @exon_starts;
	    @exon_ends = sort {$a <=> $b} @exon_ends;

	    # output
	    print GENE join("\t",
			    $transcript_id,
			    "chr$chr",
			    $strand,
			    $start,
			    $end,
			    $cds_start,
			    $cds_end,
			    scalar @$exons,
			    join(',',@exon_starts,''),
			    join(',',@exon_ends,'')
		), "\n";	
	}
	else {
	    # sort exons and check start, end
	    my (@rel_exon_starts, @exon_sizes);
	    foreach my $e (sort { $a->[0] <=> $b->[0] } @$exons) {
		push @rel_exon_starts, $e->[0] - $start - 1;
		push @exon_sizes, $e->[1] - $e->[0] + 1;
	    }
	    #print STDERR $rel_exon_starts[0], " ", $start+$rel_exon_starts[-1]+$exon_sizes[-1], " ", $end, "\n";
	    #print STDERR join("\t", $start, $end, join(",", @rel_exon_starts), join(",", @exon_sizes)), "\n";

	    print NC join("\t",
			  $BINNER->bin_from_coord_range($start+1, $end),
			  "chr$chr",
			  $start,
			  $end,
			  $transcript_id,
			  0,
			  $strand,
			  $start,
			  $start,
			  0,
			  scalar @$exons,
			  join(',',@exon_sizes,''),
			  join(',',@rel_exon_starts,'')
		), "\n";
	    
	    
	}
	
	print XREF join("\t",
			$transcript_id,
			# ^ can't use transcript symbol here because it is not unique in case-insensitive comparisons
			$gene_symbol,
			$transcript_symbol, # this column is really meant for the synonyms
			$transcript_id,
			$gene_id,
			'.', #fbpp
			'.', #fban
			$type), "\n";
	
    }
}


sub write_gbrowse_gff
{
    my ($transcripts, $genes, $source) = @_;

    print STDERR "Creating files: flybase-gbrowse.gff...\n";
    open GBROWSE, ">flybase-gbrowse.gff";

    while(my ($gene_id, $g) = each %$genes) {
	# check that we have all data
	my $gene_symbol = $g->gene_symbol or die "no symbol for $gene_id";
	my $chr = $g->chr or die "no chr for $gene_id";
	my $start = $g->start or die "no start for $gene_id";
	my $end = $g->end or die "no end for $gene_id";
	my $strand = $g->strand or die "no strand for $gene_id";
	my $group = "Gene $gene_id ; Alias \"$gene_symbol\"";
	print GBROWSE join("\t", "chr$chr", $source, 'gene', $start, $end,'.', $strand, '.',$group), "\n";
    }

    while(my ($transcript_id, $t) = each %$transcripts) {

	# check that we have all data
	my $gene_id = $t->gene_id or die "no gene_id for $transcript_id";
	my $transcript_symbol = $t->transcript_symbol or die "no symbol for $transcript_id";
	my $type = $t->type or die "no type for $transcript_id";
	my $exons_ref = $t->exons or die "no exons for $transcript_id";
	my $chr = $t->chr or die "no chr for $transcript_id";
	my $start = $t->start or die "no start for $transcript_id";
	my $end = $t->end or die "no end for $transcript_id";
	my $strand = $t->strand or die "no strand for $transcript_id";
	
	$chr = "chr$chr";
	my $group = "Transcript $transcript_id";
	my $group_w_symbol = "$group ; Symbol \"$transcript_symbol\"";
	
	# Get sorted exons
	my @exons = sort { $a->[0] <=> $b->[0] } @$exons_ref;
	
	if($t->cds_start) {
	    warn "Warning: transcript $transcript_id has CDS coordinates but is not classified as mRNA." if($type ne 'mRNA');
	    
	    # Get CDS start, end
	    my $cds_start = $t->cds_start;
	    my $cds_end = $t->cds_end or die "no cds end for $transcript_id";

	    # Compute CDS and UTR coords
	    my $fputr = AT::Tools::RangeHandler->compute_intersection(\@exons, [[$start, $cds_start-1]]);
	    my $cds = AT::Tools::RangeHandler->compute_intersection(\@exons, [[$cds_start, $cds_end]]);
	    my $tputr = AT::Tools::RangeHandler->compute_intersection(\@exons, [[$cds_end+1, $end]]);
	    ($fputr, $tputr) = ($tputr, $fputr) if ($strand eq '-');

	    # Print feature spanning entire transcript
	    print GBROWSE join("\t", $chr, $source, 'mRNA', $start, $end,'.', $strand, '.',$group_w_symbol), "\n";

	    # Print CDS and UTR features
	    print_range_set_as_gff(\*GBROWSE, $fputr, $chr, $source, "5'-UTR", '.', $strand, '.', $group);
	    print_range_set_as_gff(\*GBROWSE, $cds, $chr, $source, "CDS", '.', $strand, '.', $group);
	    print_range_set_as_gff(\*GBROWSE, $tputr, $chr, $source, "3'-UTR", '.', $strand, '.', $group);
	}
	else {
	    print GBROWSE join("\t", $chr, $source, 'transcript', $start, $end,'.', $strand, '.',$group_w_symbol), "\n";
	    print_range_set_as_gff(\*GBROWSE, \@exons, $chr, $source, "exon", '.', $strand, '.', $group);
	}
    }
}

    
sub get_parents
{
    my $entire_str = shift;
    my ($parent_str) = $entire_str =~ /Parent=([^;]+);/;
    my @parents = split ',', $parent_str;
    return @parents;
}


sub get_transcript
{
    my ($index, $id) = @_;
    return($index->{$id}) if($index->{$id});
    return($index->{$id} = transcript->new());
}

sub print_range_set_as_gff {
    my ($fh, $ranges, $chr, $source, $type, $score, $strand, $phase, $group) = @_;
    foreach my $range (@$ranges) {
	print $fh join("\t",
		       $chr, $source, $type, $range->[0], $range->[1],
		       $score, $strand, $phase, $group), "\n";
    }
    # note: phase will not be correct unless all ranges divisible by 3
}
