#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Std;
use DBI;
use Sys::Hostname;
use AT::DB::GenomeAssemblyTwoBit;
use AT::Tools::RangeHandler;

my %args;
getopts('a:cd:h:', \%args);
my $TWOBIT_FN = $args{'a'};
my $STRIP_CHR = $args{'c'};
my $UCSC_DB_NAME = $args{'d'};
my $DB_HOST = $args{'h'} || 'localhost';

my (@FEATURES) = @ARGV;

unless(@FEATURES) {

    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Usage: 
    
$cmd [options] <features>

Required arguments:

<features>  White-space separated list of feature types to read.
  The following feature types are recognized:
  assembly - Assembly sequences (typically chromosomes or scaffolds)
  refGene - RefSeq genes (refGene+refLink tables)
  knownGene - UCSC Known Genes (knownGene+knownIsoforms+kgXref tables)          
  wormBaseGene - WormBase Genes (sangerGene+sangerIsoforms+orfToGene tables)
  blast<asm>KG:<db2>:<desc> - Blasted UCSC Known Gene peptides from assembly <asm>
                 (tables: blast<asm>KG from primary db + kgXref from <db2>)             
  zfin:<file> - ZFIN Genes. ZFIN genes from tables vegaGene + vegaInfoZfish + 
                refGene + refLink and <file> which should map ZFIN ids to Entrez
                gene ids. Presently the file is at:
                http://zfin.org/data_transfer/Downloads/entrezgene.txt
  miRna - miRNAs (miRNA table)
  wgRna - snoRNAs and miRNAs (wgRna table)
  cpgIsland - CpG islands (cpgIslandExt table)
  rmsk - RepeatMasker repeats (*_rmsk tables)
  gap - Assembly gaps (*_gap table)
  oreganno - Regulatory annotation (oreganno+oregannoAttr tables)

If assembly features are requested a twobit file must be specified in options.
For the other features, a UCSC database name must be specified in options.
Note that any gbrowse database must contain assembly features.

Options:

-a <twobit file>     Name of .2bit file for assembly.
-c                   Strip the prefix "chr" from chromosome names
-d <db name>         Name of UCSC database.
-h <database host>   Database host. Default: localhost.

EndOfUsage

    exit;
}

# Connect to UCSC db if specified
my $ucsc_dbh;
if($UCSC_DB_NAME) {
    $ucsc_dbh = DBI->connect("dbi:mysql:host=$DB_HOST;database=$UCSC_DB_NAME","nobody","")
	or die "could not connect to db $UCSC_DB_NAME @ $DB_HOST";
}

# Open twobit file if specified
my $asm;
if($TWOBIT_FN) {
    $asm = AT::DB::GenomeAssemblyTwoBit->new(file => $TWOBIT_FN)
	or die "could not open assembly file $TWOBIT_FN";
}

foreach my $f (@FEATURES) {
    if($f eq 'assembly') {
       # Print out reference (assembly) sequences in GFF
	die "Need twobit file to print assembly sequences" unless($asm);
	print_assembly_sequences_as_gff($asm);
    }
    else {
	die "Need UCSC db name to get features" unless($ucsc_dbh);
	if($f eq 'refGene') {
	    # Print out Refseq gene features in GFF
	    print_refGene_features_as_gff($ucsc_dbh, 'RefSeq');
	}
	elsif($f eq 'knownGene') {
	    # Print out UCSC gene features in GFF
	    print_ucscGene_features_as_gff($ucsc_dbh, 'UCSC');
	}
	elsif($f eq 'wormBaseGene') {
	    # Print out WormBase gene features in GFF
	    print_wormBaseGene_features_as_gff($ucsc_dbh, 'WormBase');
	}
	elsif($f =~ /^(blast\w+KG):(\w+):(.+)$/) {
	    # Print out xeno-blasted UCSC peptides in GFF
	    my $table = $1;
	    my $ucsc_db2_name = $2;
	    my $desc = $3;
	    #my $ucsc_dbh2 = DBI->connect("dbi:mysql:host=$DB_HOST;database=$ucsc_db2_name","at_read","tittut")
		#or die "could not connect to db $ucsc_db2_name @ $DB_HOST";
	    print_blastKG_features_as_gff($ucsc_dbh, $ucsc_db2_name, $table, $desc);
	}
	elsif($f =~ /^zfin:(.+)$/) {
	    # Print out zfin gene bounds in GFF, mapped through Entrez gene ids
	    my $zfin2entrez_fn = $1;
	    print_zfin_features_as_gff($ucsc_dbh, $zfin2entrez_fn);
	}
	elsif($f eq 'wgRna') {
	    # Print out UCSC gene features in GFF
	    print_wgRna_features_as_gff($ucsc_dbh);
	}
	elsif($f eq 'miRna') {
	    # Print out UCSC gene features in GFF
	    print_miRna_features_as_gff($ucsc_dbh);
	}
	elsif($f eq 'cpgIsland') {
	    # Print out CpG island features in GFF
	    print_cpgIsland_features_as_gff($ucsc_dbh);
	}
	elsif($f eq 'rmsk') {
	    # Print out RepeatMasker features in GFF
	    print_rmsk_features_as_gff($ucsc_dbh);
	}
	elsif($f eq 'gap') {
	    # Print assembly gap features in GFF
	    print_gap_features_as_gff($ucsc_dbh);
	}
	elsif($f eq 'oreganno') {
	    # Print oreganno features in GFF
	    print_oreganno_features_as_gff($ucsc_dbh);
	}
	else {
	    warn "Skipping unrecognized feature type $f.\n";
	}
    }
}

# Done.
print STDERR "Done.\n";
exit;

#
# SUBROUTINES
#

sub format_seq_id
{
    my $id = shift;
    $id =~ s/^chr// if($STRIP_CHR);
    return $id;
}


sub print_assembly_sequences_as_gff
{
    my $asm = shift;
    my @seq_names = $asm->get_chr_names();
    foreach my $seq_name (@seq_names) {
	my $size = $asm->get_chr_size($seq_name);
	my $seq_type;
	my $group;
	if($seq_name =~ /^chr/) {
	    $seq_type = "chromosome";
	    my $seq_alias = $seq_name;
	    $seq_alias =~ s/^chr//;
	    ($seq_name, $seq_alias) = ($seq_alias, $seq_name) if($STRIP_CHR);
	    $group = "Sequence $seq_name ; Alias $seq_alias";
	}
	else {
	    $seq_type = "supercontig";
	    $group = "Sequence $seq_name";
	}
	print join("\t", $seq_name, 'assembly', $seq_type, 1, $size, '.', '+', '.', $group), "\n";
    }
}


sub print_refGene_features_as_gff {
    my ($dbh, $source) = @_;

    my $query = "SELECT distinct locusLinkId, refGene.name, refGene.name2, 'RefSeq Gene', chrom, strand, exonStarts, exonEnds, cdsStart, cdsEnd
                 FROM refGene, refLink WHERE refGene.name = refLink.mrnaAcc
                 ORDER BY locusLinkId, name2, name";
    my $sth = $dbh->prepare($query) or die "could not prepare query [$query]";
    $sth->execute();
    print_gene_features_as_gff($sth, $source, 'gene', 'rg');
}


sub print_zfin_features_as_gff {
    my ($dbh, $fn) = @_;

    # Read entrez to zfin and zfin to symbol associations from file
    my %zfin_to_symbol;
    my %entrez_to_zfin;
    open IN, $fn or die "could not open file $fn for reading";
    while(my $line = <IN>) {
	chomp $line;
	my ($zfin_id, $symbol, $entrez_id) = split /\t/, $line;
	if($entrez_to_zfin{$entrez_id}) {
	    warn "Multiple zfin ids for entrez id $entrez_id; will only keep first."
	}
	else {
	    $entrez_to_zfin{$entrez_id} = $zfin_id;
	    $zfin_to_symbol{$zfin_id} = $symbol;
	}
    }
    close IN;

    # Hash to store transcript bounds for each zfin id
    my %zfin_transcripts;

    # Read transcript bounds from Entrez RefSeq in UCSC db
    my $entrez_query = "SELECT locusLinkId, chrom, strand, txStart, txEnd FROM refGene, refLink
                        WHERE refGene.name = refLink.mrnaAcc";
    my $entrez_sth = $dbh->prepare($entrez_query) or die "could not prepare query [$entrez_query]";
    $entrez_sth->execute();
    while(my ($entrez_id, $chr, $strand, $tx_start, $tx_end) = $entrez_sth->fetchrow_array()) {
	my $zfin_id = $entrez_to_zfin{$entrez_id} or next;
	push @{$zfin_transcripts{$zfin_id}{$chr}{$strand}}, [$tx_start+1, $tx_end];
    }

    # Read transcript bounds from Vega table in UCSC db
    my $vega_query = "SELECT zfinId, zfinSymbol, chrom, strand, txStart, txEnd FROM vegaGene, vegaInfoZfish
                        WHERE vegaGene.name = vegaInfoZfish.transcriptId and zfinId != ''";
    my $vega_sth = $dbh->prepare($vega_query) or die "could not prepare query [$vega_query]";
    $vega_sth->execute();
    while(my ($zfin_id, $symbol, $chr, $strand, $tx_start, $tx_end) = $vega_sth->fetchrow_array()) {
	push @{$zfin_transcripts{$zfin_id}{$chr}{$strand}}, [$tx_start+1, $tx_end];
	$zfin_to_symbol{$zfin_id} = $symbol unless($zfin_to_symbol{$zfin_id});
    }

    # Collapse transcript bounds into gene bounds and output
    while(my ($zfin_id, $bounds_by_chr_strand) = each %zfin_transcripts) {
	while(my ($chr, $bounds_by_strand) = each %$bounds_by_chr_strand) {
	    while(my ($strand, $bound_list) = each %$bounds_by_strand) {
		my $collapsed_bounds = collapse_regions($bound_list);
		foreach my $bounds (@$collapsed_bounds) {
		    my $gene_group =
			'Gene '.$zfin_id.
			' ; Alias "'.($zfin_to_symbol{$zfin_id}||$zfin_id).'"'.
			' ; Note "ZFIN Gene"';
			print join("\t", $chr, 'ZFIN', 'gene', @$bounds,'.', $strand, '.',$gene_group), "\n";
		}
	    }
	}
    }

}


sub print_ucscGene_features_as_gff {
    my ($dbh, $source) = @_;

    my $query = "SELECT distinct clusterId, kgID, geneSymbol, 'UCSC Gene', chrom, strand, exonStarts, exonEnds, cdsStart, cdsEnd
                 FROM knownGene, kgXref, knownIsoforms WHERE knownGene.name = kgXref.kgID and knownGene.name=knownIsoforms.transcript
                 ORDER BY clusterId, geneSymbol, kgID";
    my $sth = $dbh->prepare($query) or die "could not prepare query [$query]";
    $sth->execute();
    print_gene_features_as_gff($sth, $source, 'gene', 'ug');
}


sub print_wormBaseGene_features_as_gff {
    my ($dbh, $source) = @_;

    my $query = "SELECT distinct clusterId, sangerGene.name, orfToGene.value, 'WormBase Gene', chrom, strand, exonStarts, exonEnds, cdsStart, cdsEnd
                 FROM sangerGene, orfToGene, sangerIsoforms WHERE sangerGene.name = orfToGene.name and sangerGene.name=sangerIsoforms.transcript
                 ORDER BY clusterId, orfToGene.value, sangerGene.name";
    my $sth = $dbh->prepare($query) or die "could not prepare query [$query]";
    $sth->execute();
    print_gene_features_as_gff($sth, $source, 'gene', 'ug');
}



sub print_blastKG_features_as_gff {
    my ($dbh, $db2_name, $table, $desc) = @_;

# get results ordered by symbols, accn
# A few accessions have multiple symbols. Since they are few (~100) and acconting for multiple
# symbols complicate the processing, we just grab one symbol here. This should not be a serious issue.
# We could do: group_concat(distinct geneSymbol separator '__') as symbols
# But then we have to split it in print_gene_features_as_gff()

    my $query = "SELECT 0, qName, geneSymbol, '$desc', tName, strand, tSize, tStarts, blockSizes
                 FROM $table LEFT JOIN $db2_name.kgXref ON qName = mRNA
                 ORDER BY geneSymbol, qName";
    my $sth = $dbh->prepare($query) or die "could not prepare query [$query]";
    $sth->execute();
    print_gene_features_as_gff($sth, $table, 'psl_prot', 'xg');
}


sub print_miRna_features_as_gff {
    my ($dbh) = @_;

    my $query = "SELECT distinct chrom, strand, chromStart, chromEnd, name FROM miRNA";
    my $sth = $dbh->prepare($query) or die "could not prepare query [$query]";
    $sth->execute();
    while(my ($chr, $strand, $start, $end, $name) = $sth->fetchrow_array) {
	my $group = "miRNA \"$name\"";
	$chr = format_seq_id($chr);
	print join("\t", $chr, 'UCSC', 'small_rna', $start+1, $end, '.', $strand, '.', $group), "\n";
    }
}

sub print_wgRna_features_as_gff {
    my ($dbh) = @_;

    my $query = "SELECT distinct chrom, strand, chromStart, chromEnd, name, type FROM wgRna";
    my $sth = $dbh->prepare($query) or die "could not prepare query [$query]";
    $sth->execute();
    while(my ($chr, $strand, $start, $end, $name, $type) = $sth->fetchrow_array) {
	$type =~ s/Rna$/RNA/;
	my $group = "$type \"$name\"";
	$chr = format_seq_id($chr);
	print join("\t", $chr, 'UCSC', 'small_rna', $start+1, $end, '.', $strand, '.', $group), "\n";
    }
}

sub print_cpgIsland_features_as_gff {
    my ($dbh) = @_;

    my $query = "SELECT distinct chrom, chromStart, chromEnd, cpgNum, perCpg, perGc FROM cpgIslandExt";
    my $sth = $dbh->prepare($query) or die "could not prepare query [$query]";
    $sth->execute();
    my $i = 1;
    while(my ($chr, $start, $end, $cpgNum, $perCpg, $perGc) = $sth->fetchrow_array) {
	my $group = "CpG_island $i ; perGC $perGc ; perCpG $perCpg";
	$chr = format_seq_id($chr);
	print join("\t", $chr, 'UCSC', 'CpG_island', $start+1, $end, '.', '+', '.', $group), "\n";
	$i++;
    }
}


sub print_rmsk_features_as_gff {
    my ($dbh, $table) = @_;

    my @tables = grep /[_\W]rmsk.?$/, $dbh->tables(undef, undef, undef, 'TABLE');
    print STDERR "Reading rmsk features from ".scalar(@tables)." tables.\n";
    my $i = 1;
    foreach my $table (@tables) {
	my $query = "SELECT genoName, genoStart, genoEnd, strand, milliDiv, repName, repClass from $table";
	my $sth = $dbh->prepare($query) or die "could not prepare query [$query]";
	$sth->execute();
	while(my ($chr, $start, $end, $strand, $divergence, $name, $class) = $sth->fetchrow_array) {
	    my $group = "$class \"$name\"";
	    $chr = format_seq_id($chr);
	    print join("\t", $chr, 'RepeatMasker', 'repeat', $start+1, $end, $divergence/10,  $strand, '.', $group), "\n";
	    $i++;
	}
    }
}


sub print_gap_features_as_gff {
    my ($dbh, $table) = @_;

    my @tables = grep /[_\W]gap.?$/, $dbh->tables(undef, undef, undef, 'TABLE');
    print STDERR "Reading gap features from ".scalar(@tables)." tables.\n";
    my $i = 1;
    foreach my $table (@tables) {
	my $query = "SELECT distinct chrom, chromStart, chromEnd from $table";
	my $sth = $dbh->prepare($query) or die "could not prepare query [$query]";
	$sth->execute();
	while(my ($chr, $start, $end) = $sth->fetchrow_array) {
	    my $group = "Gap $i";
	    $chr = format_seq_id($chr);
	    print join("\t", $chr, 'assembly', 'gap', $start+1, $end, '.',  '+', '.', $group), "\n";
	    $i++;
	}
    }
}


sub print_oreganno_features_as_gff {
    my ($dbh, $table) = @_;

    my $query1 = "SELECT a.id, chrom, chromStart, chromEnd, strand, b.attrVal
                   FROM oreganno a, oregannoAttr b 
                  WHERE a.id = b.id and b.attribute = 'type'";
    my $query2 = "SELECT attrVal
                    FROM oregannoAttr
                   WHERE id = ? and attribute = 'Dataset'";
    my $sth1 = $dbh->prepare($query1) or die "could not prepare query [$query1]";
    my $sth2 = $dbh->prepare($query2) or die "could not prepare query [$query2]";
    $sth1->execute();
    while(my ($id, $chr, $start, $end, $strand, $type) = $sth1->fetchrow_array) {
	if($type eq 'REGULATORY REGION') { $type = 'Regulatory region' }
	elsif($type eq 'TRANSCRIPTION FACTOR BINDING SITE') { $type = 'TFBS' }
	elsif($type eq 'REGULATORY POLYMORPHISM') { $type = 'Regulatory polymorphism' }
	$sth2->execute($id);
	my ($dataset) = $sth2->fetchrow_array();
	$sth2->finish;
	my $group = "ORegAnno $id";
	$group .= "; Dataset \"$dataset\"" if($dataset);
	$chr = format_seq_id($chr);
	print join("\t", $chr, 'ORegAnno', $type, $start+1, $end, '.',  $strand, '.', $group), "\n";
    }
}


sub print_gene_features_as_gff 
{
    my ($sth, $source, $format, $id_prefix) = @_;

    my $prev_tx_id = '';
    my $mapping_nr = 0;
    my @next_row = $sth->fetchrow_array;
    my %bounds_by_chr_strand;
    my $generated_gene_id = 1;

    while(my ($gene_id,$tx_id,$symbol,$note,$chr,@coord_fields) = @next_row) {

	# strip 'chr' from chromosome id requested
	$chr = format_seq_id($chr);

	# Number mappings when several have the same id
	@next_row = $sth->fetchrow_array;
	my $next_gene_id = $next_row[0] || '';
	my $next_tx_id = $next_row[1] || '';
	my $next_symbol = $next_row[2] || '';
	if($tx_id eq $prev_tx_id) {
	    $mapping_nr++;
	}
	elsif($tx_id eq $next_tx_id) {
	    $mapping_nr = 1;
	}
	else {
	    $mapping_nr = 0;
	}
	$prev_tx_id = $tx_id;

	# Assemble the group field
	my $group = $mapping_nr ? "Transcript \"$tx_id mapping $mapping_nr\"" : "Transcript $tx_id";
	my $group_w_symbol = $symbol ? "$group ; Symbol \"$symbol\"" : $group;

        # Get exon coords
	my ($exons, $strand, $start, $end, $cdsStart, $cdsEnd);
	if($format eq 'gene') {
	    ($exons, $strand, $start, $end, $cdsStart, $cdsEnd) = get_gene_exon_coords(@coord_fields);
	}
	elsif($format eq 'psl_prot') {
	    ($exons, $strand, $start, $end, $cdsStart, $cdsEnd) = get_psl_prot_exon_coords(@coord_fields);
	}
	else {
	    die "unknown format $format";
	}

	if($cdsEnd - $cdsStart > 0) {
	    # Compute CDS and UTR coords
	    my $fputr = AT::Tools::RangeHandler->compute_intersection($exons, [[$start, $cdsStart-1]]);
	    my $cds = AT::Tools::RangeHandler->compute_intersection($exons, [[$cdsStart, $cdsEnd]]);
	    my $tputr = AT::Tools::RangeHandler->compute_intersection($exons, [[$cdsEnd+1, $end]]);
	    ($fputr, $tputr) = ($tputr, $fputr) if ($strand eq '-');

	    # Print feature spanning entire transcript
	    print join("\t", $chr, $source, 'mRNA', $start, $end,'.', $strand, '.',$group_w_symbol), "\n";

	    # Print CDS and UTR features
	    print_range_set_as_gff($fputr, $chr, $source, "5'-UTR", '.', $strand, '.', $group);
	    print_range_set_as_gff($cds, $chr, $source, "CDS", '.', $strand, '.', $group);
	    print_range_set_as_gff($tputr, $chr, $source, "3'-UTR", '.', $strand, '.', $group);
	}
	else {
	    print join("\t", $chr, $source, 'transcript', $start, $end,'.', $strand, '.',$group_w_symbol), "\n";
	    print_range_set_as_gff($exons, $chr, $source, "exon", '.', $strand, '.', $group);
	}

	# Keep track of transcript bounds.
	# Collapse them into gene bounds and output if we have reached the last transcript for this gene.
	push @{$bounds_by_chr_strand{$chr}{$strand}}, [$start, $end];
	if(($gene_id and $gene_id ne $next_gene_id) or !$symbol or $symbol ne $next_symbol) {
	    # Collapse and count
	    my $nr_genes = 0;
	    foreach my $bounds_by_strand (values %bounds_by_chr_strand) {
		foreach my $strand ('+','-') {
		    next unless($bounds_by_strand->{$strand});
		    $bounds_by_strand->{$strand} = collapse_regions($bounds_by_strand->{$strand});
		    $nr_genes += scalar @{$bounds_by_strand->{$strand}}
		}
	    }
	    my $i = 1;
	    while(my ($chr, $bounds_by_strand) = each %bounds_by_chr_strand) {
		while(my ($strand, $bound_list) = each %$bounds_by_strand) {
		    foreach my $bounds (@$bound_list) {
			#my $gene_group = $nr_genes == 1 ? "Gene \"$symbol\"" : "Gene \"$symbol mapping $i\"";
			my $gene_group =
			    'Gene '.$id_prefix.($gene_id||$generated_gene_id).
			    ' ; Alias "'.($symbol||$tx_id).'"'.
			    ' ; Note "'.$note.'"';
			print join("\t", $chr, $source, 'gene', @$bounds,'.', $strand, '.',$gene_group), "\n";
			$i++;
			$generated_gene_id++;
		    }
		}
	    }
	    %bounds_by_chr_strand = ();
	}
    }

}


sub get_gene_exon_coords {
    my ($strand,$exonStart_str,$exonEnd_str,$cdsStart,$cdsEnd) = @_;
    my @exonStarts = map { $_ + 1 } split ',', $exonStart_str;
    my @exonEnds = split ',', $exonEnd_str;
    my $start = $exonStarts[0];
    my $end = $exonEnds[-1];
    my @exons;
    while(@exonStarts) {
	push @exons, [shift @exonStarts, shift @exonEnds];
    }
    return (\@exons, $strand, $start, $end, $cdsStart, $cdsEnd);
}


sub get_psl_prot_exon_coords {
    my ($strand, $chr_size, $blockStart_str, $blockSize_str) = @_;
    my $query_is_prot = 1;
    my @blockStarts = split ',', $blockStart_str;
    my @blockSizes = split ',', $blockSize_str;
    @blockSizes = map { $_ * 3 } @blockSizes if($query_is_prot);
    my ($qStrand, $tStrand) = split(//, $strand);
    $qStrand = $qStrand eq '-' ? -1 : 1;
    $tStrand = $tStrand eq '-' ? -1 : 1;
    # Note: in the above calculations, it is important that we default to 1, since
    # tStrand may be an empty string and should then be interpreted as '+'.
    $strand = $qStrand * $tStrand == 1 ? '+' : '-';
    if($tStrand == -1) {
	@blockStarts = reverse @blockStarts;
	@blockSizes = reverse @blockSizes;
	for my $i (0..@blockStarts-1) {
	    $blockStarts[$i] = $chr_size - $blockStarts[$i] - $blockSizes[$i]; 
	}
    }
    my $ftStart = $blockStarts[0]+1;
    my $ftEnd = $blockStarts[-1]+$blockSizes[-1];
    my @exons;
    while(@blockStarts) {
	my $blockStart = shift @blockStarts;
	my $blockSize = shift @blockSizes;
	push @exons, [$blockStart+1, $blockStart+$blockSize];
    }
    return (\@exons, $strand, $ftStart, $ftEnd, $ftStart, $ftEnd);
}


# For each gene we need:
# Symbol, start, end, strand, mapping nr (if several)
# We need to group them by symbol+strand+overlap
# So:
#  - index s,e by symbol+strand
#  - for each s,e list in the index: collapse ranges. output each resulting range as a gene.

sub print_range_set_as_gff {
    my ($ranges, $chr, $source, $type, $score, $strand, $phase, $group) = @_;
    foreach my $range (@$ranges) {
	print join("\t",
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
