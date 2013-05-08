use warnings;
use strict;
use DBI;
use Sys::Hostname;
use Getopt::Std;
use AT::DB::GenomeAssemblyTwoBit;
use MyPerlVars;

my %args;
getopts('a:b:g:p:f:c:r',\%args);

my $BUILD = $args{'b'};
my $DB_HOST = (hostname() =~ /^angband\./) ? 'localhost' : 'localhost';
my $DB_PORT = (hostname() =~ /^angband\./) ? 3306 : 3306;

my $GENE_FEATURES = $args{'g'} || "";
my $PSL_FEATURES = $args{'p'} || "";
my $EXON_FILE = $args{'f'};
my $CHROM_FILE = $args{'c'};
my $ASSEMBLY_FILE = $args{'a'};
my $REPEATS = $args{'r'};

#usage
unless($BUILD and ($GENE_FEATURES or $PSL_FEATURES or $EXON_FILE or $REPEATS)) {
    print <<EndOfUsage;

This script creates a non-redundant exon representation by calculating the union
of exons from one or more UCSC gene tables, for subsequent use in CNE filtering stages

Usage:
perl $0 [options]

Options:
 -b   Assembly name. Specifies which of the UCSC_<assembly> databases is read
 -g   UCSC gene table in standard gene format
 -p   UCSC gene table in psl format (mapped sequence). Coordinates may be protein based!
 -f   Tab delimited file listing chr, start, end for exons
 -c   Tab delimited file listing all chromosome names for the species used
 -a   Twobit assembly file (for chromosome names)
 -r   Include repeats in the filter set
      
EndOfUsage
   exit;
}

my $UCSC_DB_NAME = "UCSC_$BUILD";
my @UCSC_GENE_TABLES = split ",", $GENE_FEATURES if $GENE_FEATURES;
my @UCSC_PSL_TABLES = split ",", $PSL_FEATURES if $PSL_FEATURES;

#some gene tables give exon features in protein coordinates..
my %isProtein = ('all_est' => 1,
		 'xenoMrna' => 1,
		 'blastHg18KG' => 1
		 );

# Connect to database
my $dbh = DBI->connect("dbi:mysql:host=$DB_HOST;database=$UCSC_DB_NAME;port=$DB_PORT",$MyPerlVars::sqlUser,$MyPerlVars::sqlPass)
    or die "could not connect to db $UCSC_DB_NAME @ $DB_HOST:$DB_PORT";

#bed head
print 'track name="collExons" ';
print "description=\"Collapsed $BUILD features from ";
print "${GENE_FEATURES}," if ${GENE_FEATURES};
print "${PSL_FEATURES}," if ${PSL_FEATURES};
print "${EXON_FILE}," if ${EXON_FILE};
print "repeatMasker" if $REPEATS;
print '" ';
print "color=0,0,120\n";

my $features_by_chr = read_exons_from_bedfile($EXON_FILE) if ($EXON_FILE);

#read chromosome list
my @chrom_list;
#my $chrom_list = read_chrom_list($CHROM_FILE) if ($CHROM_FILE);

my $genome_db = AT::DB::GenomeAssemblyTwoBit->new (
						   file => $ASSEMBLY_FILE,
						   id => 'test'
						   ) if $ASSEMBLY_FILE;
@chrom_list = $genome_db->get_chr_names();

#my $chr_sth = $dbh->prepare("select distinct chrom from genscan"); #get nonred list of chromosomes
#my $chr_sth = $dbh->prepare("select distinct chrom from all_mrna"); #get nonred list of chromosomes
#$chr_sth->execute();
#while (my ($chr) = $chr_sth->fetchrow_array) {
#    push @chrom_list, $chr;
#}

#foreach my $chr (@$chrom_list){
foreach my $chr (@chrom_list){
    print STDERR "processing chr $chr...\n";
    my $features = $features_by_chr->{$chr} || [];
    foreach my $gene_table (@UCSC_GENE_TABLES){
	get_gene_features($dbh, $gene_table, $chr, $features);
    }

    foreach my $psl_table (@UCSC_PSL_TABLES){
	get_psl_features($dbh, $psl_table, $chr, $features, $isProtein{$psl_table} ? 1 : 0);
    }
    if($REPEATS){
	get_repeat_features($dbh, $chr, $features);
    }

    my $collapsed = calc_expressed_intervals($features);
    foreach my $iv (@$collapsed) {
	print join("\t", $chr, $iv->[0]-1, $iv->[1]), "\n";
    }
    $features = undef; # free memory
}

sub get_repeat_features
{
    my ($dbh, $chr, $list) = @_;

    # Determine options to select statement
#    my $table_exists = grep /^$table$/, $dbh->tables(undef, undef, undef, 'TABLE');
    my $table = $chr . "_rmsk";
    my $table_exists = grep /.*$table.*/, $dbh->tables(undef, undef, undef, 'TABLE');
    $table = "rmsk" unless ($table_exists);
    
    my $query = "SELECT genoStart + 1,genoEnd FROM $table WHERE genoName = ?";
    my $sth = $dbh->prepare($query) or die "could not prepare query [$query]";
    $sth->execute($chr);

    while(my ($repeatStart,$repeatEnd) = $sth->fetchrow_array) {
	push @$list, [$repeatStart, $repeatEnd];
    }
}


sub calc_expressed_intervals
{
    my ($features_ref) = @_;

    return [] unless(@$features_ref);

    my @hsps = sort {$a->[0] <=> $b->[0]} @$features_ref;

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
    print STDERR "collapsed ", scalar(@$features_ref)," into ",scalar(@iv)," intervals\n";
    return \@iv;
}


sub get_gene_features
{
    my ($dbh, $table, $chr, $list) = @_;

    my $query = "SELECT exonStarts, exonEnds FROM $table WHERE chrom = ?";
    my $sth = $dbh->prepare($query) || die "Could not prepare $query\n";
    $sth->execute($chr) || die "Could not execute $query\n";
    while(my ($exonStart_str,$exonEnd_str) = $sth->fetchrow_array) {
	my @exonStarts = map { $_ + 1 } split /,/, $exonStart_str;
	my @exonEnds = split /,/, $exonEnd_str;
	while(@exonStarts) {
	    push @$list, [shift @exonStarts, shift @exonEnds];
	}
    }
}


sub get_psl_features {

    my ($dbh, $table, $chr, $list, $query_is_prot) = @_;

    my $query = "SELECT strand, tSize, tStarts, blockSizes FROM $table WHERE tName = ?";
    my $sth = $dbh->prepare($query) or die "could not prepare query [$query]";
    $sth->execute($chr);

    while(my ($strand, $chr_size, $blockStart_str, $blockSize_str) = $sth->fetchrow_array) {
	my @blockStarts = split /,/, $blockStart_str;
	my @blockSizes = split /,/, $blockSize_str;
	@blockSizes = map { $_ * 3 } @blockSizes if($query_is_prot);
	my ($qStrand, $tStrand) = split(//, $strand);
	$qStrand = $qStrand eq '-' ? -1 : 1;
	$tStrand = $tStrand eq '-' ? -1 : 1;
        # Note: in the above calculations, it is important that we default to 1, since
        # tStrand may be an empty string and should then be interpreted as '+'.
	$strand = $qStrand * $tStrand;
	if($tStrand == -1) {
	    @blockStarts = reverse @blockStarts;
	    @blockSizes = reverse @blockSizes;
	    for my $i (0..@blockStarts-1) {
		$blockStarts[$i] = $chr_size - $blockStarts[$i] - $blockSizes[$i]; 
	    }
	}
	while(@blockStarts) {
	    my $blockStart = shift @blockStarts;
	    my $blockSize = shift @blockSizes;
	    push @$list, [$blockStart+1, $blockStart+$blockSize];
	}
    }
}


sub read_exons_from_bedfile
{
    my ($fn) = @_;
    open IN, $fn or die "could not open $fn\n";
    my %exons;
    while(my $line = <IN>) {
	chomp $line;
	my ($chr, $exon_start, $exon_end) = (split /\t/, $line)[0,1,2];
	push @{$exons{$chr}}, [$exon_start, $exon_end];
    }
    return \%exons;
}

sub read_chrom_list
{
    my ($fn) = @_;
    open IN, $fn or die "could not open $fn\n";
    my %chrom;
    while(my $line = <IN>) {
	chomp $line;
	my ($chr) = split /\t/, $line;
	$chrom{$chr} = 1;
    }
    my @chromlist = sort keys %chrom;
    return \@chromlist;
}


#sub pslIsProtein {
#    my ($psl) = @_;

#       boolean pslIsProtein(const struct psl *psl)
#        /* is psl a protein psl (are it's blockSizes and scores in protein space) 
#        */
#        {
#        int lastBlock = psl->blockCount - 1;
#        
#        return  (((psl->strand[1] == '+' ) &&
#             (psl->tEnd == psl->tStarts[lastBlock] + 3*psl->blockSizes[lastBlock])) 
#        ||
#            ((psl->strand[1] == '-') &&
#             (psl->tStart == (psl->tEnd-(psl->tStarts[lastBlock] + 3*psl->blockSizes[lastBlock]))))); }

#}
