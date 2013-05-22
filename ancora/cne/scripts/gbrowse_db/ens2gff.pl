use warnings;
use strict;
use Class::Struct;
use AT::DB::GenomeAssemblyTwoBit;
use AT::Tools::RangeHandler;
use Getopt::Std;

my $SOURCE = "Ensembl";
my @FIELD_NAMES = (
    'Ensembl Gene ID', 'Chromosome Name', 'Ensembl Transcript ID',
    'Gene Start (bp)', 'Gene End (bp)', 'Transcript Start (bp)', 'Transcript End (bp)',
    'Strand', #'External Gene ID',
    'Associated Gene Name',
    #'Exon Start (bp)', 'Exon End (bp)', 'Coding Start (bp)', 'Coding End (bp)');
    'Exon Chr Start (bp)', 'Exon Chr End (bp)', 'Genomic coding start', 'Genomic coding end');

my %args;
getopts('h',\%args);
my $ALLOW_HEADER_MISMATCH = $args{h};

my ($TWOBIT_FN, $MART_FN) = @ARGV;

unless($MART_FN) {
    
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    my $fields = join("\n ",@FIELD_NAMES);
    print STDERR <<EndOfUsage;
    
$cmd

Usage: 
    
$cmd [options] <2bit file> <biomart file>

Required arguments:

<2bit file>     Name of .2bit file for assembly
<biomart file>  Name of result file from biomart query

To generate the biomart file, go to www.ensembl.org and follow the biomart link.
Choose the most recent Ensembl version for your assembly of interest.
(Use the "View previous release of page..." link if you need an earlier version of 
Ensembl than the current.)
Under Dataset, choose genes for the organism of interest (e.g. "Homo sapiens genes").
Do not apply any filters.
Choose to get the following attributes (categorized under "Structures"):
 $fields

Options:

-h   Proceed even if biomart file header does not match the expected


EndOfUsage

    exit;
}


struct ENSG => [ensg_id => '$',
    symbol => '$',
    chr => '$',
    start => '$',
    end => '$',
    strand => '$',
    transcripts => '%',
    ];

struct ENST => [enst_id => '$',
    symbol => '$',
    # ^ Although it appears that all transcripts for the same gene
    #   have the same symbol, we allow different symbols in case this changes.
    start => '$',
    end => '$',
    exons => '@',
    cds => '@',
    ];

my %genes;

# Open input file
open IN, $MART_FN or die "could not open $MART_FN";

# Find which order the columns are in (this can vary between ensembl versions & species)
my $header = <IN>;
die "Input file is empty\n" unless(defined $header);
my $column_order = parse_header($header, \@FIELD_NAMES);


# Get list of valid chromsomes from 2bit file
my $asm = AT::DB::GenomeAssemblyTwoBit->new(file => $TWOBIT_FN)
    or die "could not open assembly file $TWOBIT_FN";
my %valid_chr = map { $_ => 1 } $asm->get_chr_names;
my %skipped_chr;


# Read all genes & transcripts
print STDERR "Reading input file...\n";
my ($nr_genes, $nr_tx)  = (0,0);  # just to have something to say at the end
my %total_ensg_counts; # count nr of genes for each gene id
my %total_enst_counts; # count nr of transcripts of each transcript id
while(my $line  = <IN>) {

    chomp $line;
    my @fields = split /\t/, $line;
    
    # Parse line
    my ($ensg_id, $chr, $enst_id, $gene_start, $gene_end, $tx_start, $tx_end, $strand, $symbol,
        $exon_start, $exon_end, $cds_start, $cds_end) = @fields[@$column_order];

    # Change chr name to UCSC-style and check that it is valid
    $chr = 'M' if($chr eq 'MT');
    if($valid_chr{"chr$chr"}) {
	$chr = "chr$chr";
    }
    elsif(!$valid_chr{$chr}) {
	$skipped_chr{$chr} = 1;
	next;
    }

    # Create gene if we have not seen this gene before
    # Because Ensembl genes in different locations (e.g. pseudoautosomal regions)
    # may have the same id, we handle these cases by indexing by chr and gene_start
    my $gene = $genes{$ensg_id}{$chr}{$gene_start};
    unless($gene) {
	$gene = ENSG->new(ensg_id => $ensg_id,
			  symbol => $symbol, 
                          # ^ it should be ok to set this to the first symbol we see for this gene
			  #   because each gene seems to have a unique symbol
			  chr => $chr,
			  start => $gene_start,
			  end => $gene_end,
			  strand => $strand,
			  transcripts => {});
	$genes{$ensg_id}{$chr}{$gene_start} = $gene;
	$total_ensg_counts{$ensg_id}++;
	$nr_genes++;
    }

    # Create transcript if we have not seen this transcript before
    my $tx = $gene->transcripts($enst_id);
    unless($tx) {
	$tx = ENST->new(enst_id => $enst_id,
			symbol => $symbol, # just in case some gene has multiple symbols
			start => $tx_start,
			end => $tx_end,
			exons => [],
		        cds => []);
	$gene->transcripts($enst_id, $tx);
	$total_enst_counts{$enst_id}++;
	$nr_tx++;
    }
    
    # Keep exon and cds positions
    push @{$tx->exons()}, [$exon_start,$exon_end];
    push @{$tx->cds()}, [$cds_start,$cds_end] if($cds_start);
}
print STDERR "Read $nr_genes genes and $nr_tx transcripts for ".scalar(keys %genes)." unique gene ids.\n";
print STDERR "Skipped genes on sequences: ".join(", ", keys %skipped_chr)." beacuse they are not in the 2bit file.\n"
  if(scalar keys %skipped_chr);


# Print genes and transcripts
foreach my $genes_by_chr (values %genes) {
    my $ensg_count = 0;
    my %enst_count;
    foreach my $genes_by_start (values %$genes_by_chr) {
	foreach my $gene (values %$genes_by_start) {
	    $ensg_count++;
	    my $ensg_id = $gene->ensg_id;

	    # Convert strand from number to character
	    my $strand = $gene->strand==1 ? '+' : $gene->strand==-1 ? '-' : die "Invalid strand ".$gene->strand;

	    # Make GFF group string for gene
	    my $gene_group = "Gene $ensg_id ; Note \"Ensembl Gene\"";
	    $gene_group .= ' ; Alias "'.$gene->symbol.'"' if($gene->symbol);

	    # Print gene as GFF
	    print join("\t", $gene->chr, $SOURCE, 'gene', $gene->start, $gene->end, '.', $strand, '.', $gene_group), "\n";

	    # Process each transcript for the gene
	    my $transcripts = $gene->transcripts;
	    foreach my $tx (values %$transcripts) {
		my $enst_id = $tx->enst_id;
		$enst_count{$enst_id}++;

		# Make GFF group string for transcript
		my $tx_group;
		if($total_enst_counts{$enst_id} > 1) {
		    $tx_group = "Transcript \"$enst_id mapping $enst_count{$enst_id}\"";
		}
		else {
		    $tx_group = "Transcript $enst_id";
		}
		my $tx_group_w_symbol = $tx->symbol ? "$tx_group ; Symbol ".$tx->symbol : $tx_group;

		# Get sorted exon and cds coordinates
		my @exons = sort { $a->[0] <=> $b->[0] } @{$tx->exons};
		my @cds = sort { $a->[0] <=> $b->[0] } @{$tx->cds};

		# Does the transcript have a CDS?
		if(@cds) {
		    # For coding transcripts: compute UTR locations, then print UTRs and CDS
		    my $fputr = AT::Tools::RangeHandler->compute_intersection(\@exons, [[$tx->start, $cds[0][0]-1]]);
		    my $tputr = AT::Tools::RangeHandler->compute_intersection(\@exons, [[$cds[-1][1]+1, $tx->end]]);
		    ($fputr, $tputr) = ($tputr, $fputr) if ($strand eq '-');
		    print join("\t", $gene->chr, $SOURCE, 'mRNA', $tx->start, $tx->end,'.', $strand, '.',$tx_group_w_symbol), "\n";
		    print_range_set_as_gff($fputr, $gene->chr, $SOURCE, "5'-UTR", '.', $strand, '.', $tx_group);
		    print_range_set_as_gff(\@cds, $gene->chr, $SOURCE, "CDS", '.', $strand, '.', $tx_group);
		    print_range_set_as_gff($tputr, $gene->chr, $SOURCE, "3'-UTR", '.', $strand, '.', $tx_group);
		}
		else {
		    # For nonocding transcripts: just print the exons
		    print join("\t", $gene->chr, $SOURCE, 'transcript', $tx->start, $tx->end,'.', $strand, '.',$tx_group_w_symbol), "\n";
		    print_range_set_as_gff(\@exons, $gene->chr, $SOURCE, "exon", '.', $strand, '.', $tx_group);
		}
	    }
	}
    }
}

print STDERR "Done.\n";


sub print_range_set_as_gff {
    my ($ranges, $chr, $source, $type, $score, $strand, $phase, $group) = @_;
    foreach my $range (@$ranges) {
	print join("\t",
		   $chr, $source, $type, $range->[0], $range->[1],
		   $score, $strand, $phase, $group), "\n";
    }
    # note: phase will not be correct unless all ranges divisible by 3
}


sub parse_header {
    my ($header, $ordered_field_list) = @_;
    chomp $header;
    my @actual_field_list = split /\t/, $header;
    my %field_index;
    my $i = 0;
    foreach my $f (@actual_field_list) {
	$field_index{$f} = $i++;
    }
    my @ordered_indices;
    foreach my $f (@$ordered_field_list) {
	$i = $field_index{$f};
	die "Header is missing field '$f'\n" unless(defined $i);
	push @ordered_indices, $i;
    }
    return \@ordered_indices;
}
