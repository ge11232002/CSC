#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use Class::Struct;
use AT::DB::Binner;
use AT::Tools::RangeHandler;

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

# Parse commandline args
my %args;
getopts(\%args, 'u:');
my $CREATE_UCSC_FILES = $args{'u'};
my ($MODE, $GFF_FN, $ID_FN) = @ARGV;
$MODE = lc $MODE;

# Check args and print usage if incorrect
unless(($MODE eq 'flybase' and $GFF_FN) or
       ($MODE eq 'wormbase' and $GFF_FN and $ID_FN)) {

    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Usage: 
    
perl $cmd [-u] flybase <flybase-gff3-file> 
perl $cmd wormbase <wormbase-gff3-file> <wormbase-id-file>

The first input file should be a GFF3 annotation file from FlyBase or WormBase.
For WormBase, an additional file that maps gene and transcript IDs to symbols is
required.

For FlyBase Release 2007_1 the file is:
ftp://ftp.flybase.net/releases/FB2007_01/dmel_r5.2_FB2007_01/gff/dmel-all-r5.2.gff.gz

For WormBase Release WS176 the files are:
ftp://ftp.wormbase.org/pub/wormbase/genomes/elegans/genome_feature_tables/GFF3/elegansWS176.gff3.gz
ftp://ftp.wormbase.org/pub/wormbase/genomes/elegans/annotations/gene_IDs/geneIDs.WS176.gz

This script will parse the input files and print, on standard output, GFF data
suitable for import into a GBrowse Bio:DB:GFF database.

If the option -u is given in FlyBase mode, the following files will also be
generated:
flyBaseGene.txt
flyBaseNoncoding.txt
flyBase2004Xref.txt
The files follow the UCSC database schema for FlyBase genes.
WARNING: these files will be overwritten if they exist!

EndOfUsage

    exit;
}

# Do the work and exit
main($MODE, $GFF_FN, $ID_FN);
exit;


sub main {
    my ($mode, $gff_fn, $id_fn) = @_;

# We need to treat FlyBase and WormBase files separately because they differ in 
# several respects. For example:
#
# FlyBase has gene symbols in comment field: Name=symbol
# WormBase has gene symbols in a separate file 
#
# FlyBase always has source FlyBase
# WormBase has sources: Coding_transcript, Non_coding_transcript, Pseudogene, *RNA
#   (where *RNA = miRNA, ncRNA, rRNA, scRNA, snlRNA, snoRNA, snRNA, tRNA)
#
# FlyBase pseudogenes as transcripts undger gene features
# WormBase has pseudogenes on top level

    my ($transcripts, $genes, $source);

    print STDERR "Reading input file(s)...\n";
    if($mode eq 'flybase') {
	($transcripts, $genes) = read_flybase_gff($gff_fn);
	$source = "FlyBase";
    }
    elsif($MODE eq 'wormbase') {
	($transcripts, $genes) = read_wormbase_gff($gff_fn, $id_fn);
	$source = "WormBase";
    }
    else {
	die "Illegal mode $mode.\n";
    }
    
    print STDERR "Validating transcripts...\n";
    validate_transcripts($transcripts);
    
    if(%$transcripts) {
	print STDERR "Writing output...\n";
	write_gbrowse_gff($transcripts, $genes, $source);
	write_ucsc_tables($transcripts, $genes) if($MODE eq 'flybase' and $CREATE_UCSC_FILES);
    }
    else {
	print STDERR "No transcript data found in the input files.\nNo files created.\n";
    }
    
    print STDERR "Done.\n";
}


sub read_flybase_gff {

    my ($gff_fn) = @_;

    my $transcripts = {};
    my $genes = {};
    my $transcript_symbols = {};

    open IN, $gff_fn or die "could not open $gff_fn";

    while(my $line = <IN>) {
	    next if($line =~ /^\#/);
	    chomp $line;
	    my ($chr, $source, $type, $start, $end, undef, $strand, undef, $extras) = 
	    split /\t/, $line;
	    next unless(defined($source) and $source eq 'FlyBase');
	    if($type eq 'gene') {
	      # Handle gene feature
	      my ($gene_id) = $extras =~ /ID=(\w+);/;
	      my ($gene_symbol) = $extras =~ /Name=([^;]+);/;
	      add_gene($genes, $gene_id, $gene_symbol, $chr, $start, $end, $strand);
	    }
  	  elsif($type eq 'CDS') {
	      # Handle CDS feature
	      my $parents = get_parents($extras);
	      add_cds($transcripts, $parents, $start, $end);
	    }
	    elsif($type eq 'exon') {
	      # Handle exon feature
	      my $parents = get_parents($extras);
	      add_exon($transcripts, $parents, $start, $end);
	    }
	    elsif($type =~ /RNA$/ or $type eq 'pseudogene') {
	      # Handle transcript / pseudogene feature (coding or noncoding)
	      # The regex is supposed to capture: mRNA, tRNA, ncRNA, rRNA, miRNA, snoRNA, scRNA, snoRNA, snRNA
	      my ($transcript_id) = $extras =~ /ID=(\w+);/;
	      my ($transcript_symbol) = $extras =~ /Name=([^;]+);/;
	      my $parents = get_parents($extras);
	      die "no gene for transcript $transcript_id" if(@$parents == 0);
	      die "multiple genes for transcript $transcript_id" if(@$parents > 1);
	      add_transcript($transcripts, $transcript_symbols, $transcript_id, $transcript_symbol,
			   $type, $chr, $start, $end, $strand, $extras, $parents->[0]);
	    }
	    # we could also look for: nc_primary_transcript, which is used in the wormbase gff3 files
    }

    close IN;

    return ($transcripts, $genes);
}


sub read_wormbase_gff {

    my ($gff_fn, $id_fn) = @_;

    my $transcripts = {};
    my $genes = {};
    my $transcript_symbols = {};

    # Read gene ID -> symbol and gene ID -> alternative gene ID assocations
    my (%gene2symbol, %gene2alias);
    open IN, $id_fn or die "could not open $id_fn";
    while(my $line = <IN>) {
	chomp $line;
	my ($gene_id, $symbol, $alias) = split ',', $line;
	$gene2symbol{$gene_id} = $symbol if($gene_id and $symbol);
	$gene2alias{$gene_id} = $alias if($gene_id and $alias);
    }
    close IN;

    open IN, $gff_fn or die "could not open $gff_fn";
    
    while(my $line = <IN>) {
	next if($line =~ /^\#/);
	chomp $line;
	my ($chr, $source, $type, $start, $end, undef, $strand, undef, $extras) = 
	    split /\t/, $line;
	next unless(defined($source) and ($source eq 'Coding_transcript' or $source eq 'Non_coding_transcript' or
					  $source eq 'Pseudogene' or $source =~ /RNA$/));
	if($type eq 'gene' or $type eq 'Pseudogene') {
	    # Handle gene / pseudogene feature
	    my ($gene_id) = $extras =~ /ID=Gene:([^;]+)/;
	    my $gene_symbol = $gene2symbol{$gene_id} || $gene2alias{$gene_id} || '';
	    add_gene($genes, $gene_id, $gene_symbol, $chr, $start, $end, $strand);
	}
	elsif($type eq 'CDS') {
	    # Handle CDS feature
	    my $parents = get_parents($extras, 'Transcript');
	    add_cds($transcripts, $parents, $start, $end);
	}
	elsif($type eq 'exon' or $type eq 'five_prime_UTR' or $type eq 'three_prime_UTR') {
	    # Handle exon feature. Note: in wormbase UTRs do not belong to exons.
	    my $parents = get_parents($extras, 'Transcript');
	    add_exon($transcripts, $parents, $start, $end);
	}
	elsif($type eq 'mRNA' or $type eq 'ncRNA') {
	    # Handle transcript feature (coding or noncoding)
	    my ($transcript_id) = $extras =~ /ID=Transcript:([^;]+)/;
	    my $parents = get_parents($extras, 'Gene');
	    die "no gene for transcript $transcript_id" if(@$parents == 0);
	    die "multiple genes for transcript $transcript_id" if(@$parents > 1);
	    my $gene_id = $parents->[0];
	    my $transcript_symbol = $gene2symbol{$gene_id} || '';
	    add_transcript($transcripts, $transcript_symbols, $transcript_id, $transcript_symbol,
			   $type, $chr, $start, $end, $strand, $gene_id);
	}
    }

    close IN;

    return ($transcripts, $genes);
}


sub add_gene
{
    my ($genes, $id, $symbol, $chr, $start, $end, $strand) = @_;
    die "duplicate gene " if($genes->{$id});
    my $gene = $genes->{$id} = gene->new();
    $gene->chr($chr);
    $gene->start($start);
    $gene->end($end);
    $gene->strand($strand);
    $gene->gene_symbol($symbol);
    return $gene;
}


sub add_cds
{
    my ($transcripts, $transcript_ids, $start, $end) = @_;
    foreach my $transcript_id (@$transcript_ids) {
	my $transcript = get_transcript($transcripts, $transcript_id);
	$transcript->cds_start($start) if(!defined($transcript->cds_start) or $start < $transcript->cds_start);
	$transcript->cds_end($end) if(!defined($transcript->cds_end) or $end > $transcript->cds_end);
    }
}


sub add_exon
{
    my ($transcripts, $transcript_ids, $start, $end) = @_;
    foreach my $transcript_id (@$transcript_ids) {
	my $transcript = get_transcript($transcripts, $transcript_id);
	push @{$transcript->exons}, [$start, $end];
    }
}


sub add_transcript 
{
    my ($transcripts, $symbols, $id, $symbol, $type, $chr, $start, $end, $strand, $gene_id) = @_;
    my $transcript = get_transcript($transcripts, $id);
    die "duplicate transcript $id" if($transcript->start);
    #die "duplicate transcript symbol $symbol" if($symbol and $symbols->{$symbol});
    $symbols->{$symbol} = 1;
    $transcript->type($type);
    $transcript->chr($chr);
    $transcript->start($start);
    $transcript->end($end);
    $transcript->strand($strand);
    $transcript->transcript_symbol($symbol);
    $transcript->gene_id($gene_id);
    return $transcript;
}


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
	my $exons = collapse_regions($t->exons);
	$t->exons($exons);
	if($exons->[0][0] != $t->start or $exons->[-1][1] != $t->end) {
	    warn "Warning: exon and transcript boundaries do not match for $transcript_id.\n";
	    if($exons->[0][0] >= $t->start and $exons->[-1][1] <= $t->end and $exons->[0][0] <= $exons->[-1][1]) {
		$t->start($exons->[0][0]);
		$t->end($exons->[-1][1]);
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
    my $binner = AT::DB::Binner->new();
    open GENE, ">flyBaseGene.txt";
    open NC, ">flyBaseNoncoding.txt";
    open XREF, ">flyBase2004Xref.txt";
    while(my ($transcript_id, $t) = each %$transcripts) {

	# check that we have all data
	my $gene_id = $t->gene_id or die "no gene_id for $transcript_id";
	my $transcript_symbol = $t->transcript_symbol || ''; # or die "no symbol for $transcript_id";
	my $type = $t->type or die "no type for $transcript_id";
	my $exons = $t->exons or die "no exons for $transcript_id";
	my $chr = $t->chr or die "no chr for $transcript_id";
	my $start = $t->start or die "no start for $transcript_id";
	my $end = $t->end or die "no end for $transcript_id";
	my $strand = $t->strand or die "no strand for $transcript_id";
	my $gene_symbol = $genes->{$t->gene_id}->gene_symbol || ''; # or die "no symbol for $gene_id";
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
			  $binner->bin_from_coord_range($start+1, $end),
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

    while(my ($gene_id, $g) = each %$genes) {
	# check that we have all data
	my $gene_symbol = $g->gene_symbol || ''; # or die "no symbol for $gene_id";
	my $chr = $g->chr or die "no chr for $gene_id";
	my $start = $g->start or die "no start for $gene_id";
	my $end = $g->end or die "no end for $gene_id";
	my $strand = $g->strand or die "no strand for $gene_id";
	my $group = "Gene $gene_id";
	$group = "$group ; Alias \"$gene_symbol\"" if($gene_symbol);
	$group = "$group ; Note \"$source Gene\"";
	print join("\t", "chr$chr", $source, 'gene', $start, $end,'.', $strand, '.',$group), "\n";
    }

    while(my ($transcript_id, $t) = each %$transcripts) {

	# check that we have all data
	my $gene_id = $t->gene_id or die "no gene_id for $transcript_id";
	my $transcript_symbol = $t->transcript_symbol || ''; # or die "no symbol for $transcript_id";
	my $type = $t->type or die "no type for $transcript_id";
	my $exons_ref = $t->exons or die "no exons for $transcript_id";
	my $chr = $t->chr or die "no chr for $transcript_id";
	my $start = $t->start or die "no start for $transcript_id";
	my $end = $t->end or die "no end for $transcript_id";
	my $strand = $t->strand or die "no strand for $transcript_id";
	
	$chr = "chr$chr";
	my $group = "Transcript $transcript_id";
	my $group_w_symbol = $transcript_symbol ? "$group ; Symbol \"$transcript_symbol\"" : $group;
	
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
	    print join("\t", $chr, $source, 'mRNA', $start, $end,'.', $strand, '.',$group_w_symbol), "\n";

	    # Print CDS and UTR features
	    print_range_set_as_gff(\*STDOUT, $fputr, $chr, $source, "5'-UTR", '.', $strand, '.', $group);
	    print_range_set_as_gff(\*STDOUT, $cds, $chr, $source, "CDS", '.', $strand, '.', $group);
	    print_range_set_as_gff(\*STDOUT, $tputr, $chr, $source, "3'-UTR", '.', $strand, '.', $group);
	}
	else {
	    print join("\t", $chr, $source, 'transcript', $start, $end,'.', $strand, '.',$group_w_symbol), "\n";
	    print_range_set_as_gff(\*STDOUT, \@exons, $chr, $source, "exon", '.', $strand, '.', $group);
	}
    }
}

    
sub get_parents
{
    my ($entire_str, $type) = @_;
    my ($parent_str) = $entire_str =~ /Parent=([^;]+)/;
    my @parents = split ',', $parent_str;
    if($type) {
	@parents = grep { ($_ =~ s/^$type://) ? $_ : '' } @parents;
    }
    return \@parents;
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

sub collapse_regions
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
    return \@iv;
}
