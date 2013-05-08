#!/usr/bin/env perl -w
use strict;
use lib "/home/engstrom/DEVEL/AT/lib";
use AT::DB::GenomeAssembly;
use AT::DB::GenomeAssemblyNibs;
use AT::DB::GenomeMapping;
use AT::DB::TranscriptSeq;
use AT::MapAlignmentFactory;
use Bio::AlignIO;
use Getopt::Std;

use constant DEFAULT_HOST => 'nautilus.cgb.ki.se';
use constant DEFAULT_USER => 'at_read';
use constant DEFAULT_PASS => 'tittut';
use constant DEFAULT_TRSEQ_DB => 'TRANSCRIPT_SEQ';

my %args;
getopts('lc:h:u:p:n:t:2', \%args);
my $LONG_ALN = $args{'l'};
my $COMBINE_ALN = $args{'c'};
my $dbhost = $args{'h'} || DEFAULT_HOST;
my $dbuser = $args{'u'} || DEFAULT_USER;
my $dbpass = $args{'p'} || DEFAULT_PASS;
my $trseqdb_name = $args{'t'} || DEFAULT_TRSEQ_DB;
my $nibdir = $args{'n'};
my $mapdb_class;
if($args{'2'}) {  # use exerimental v.2 db schema
    require "AT/DB/GenomeMapping2.pm";
    $mapdb_class = 'AT::DB::GenomeMapping2';
}
else {
    $mapdb_class = 'AT::DB::GenomeMapping';
}

my ($database, @accs) = @ARGV;
usage() unless(@accs);

my ($assembly) = $database =~ /AT_(.._...\d\d)/;

# Connect to databases

my $query_db = $mapdb_class->connect( -dbname => $database,
			    	      -dbhost => $dbhost,
				      -dbuser => $dbuser,
				      -dbpass => $dbpass);

my $target_db;
if($nibdir) {
    $target_db = AT::DB::GenomeAssemblyNibs->new(assembly_name => $assembly,
						 dir => $nibdir);
}
else {
    $target_db = AT::DB::GenomeAssembly->connect( -dbname => $assembly,
					  -dbhost => $dbhost,
					  -dbuser => $dbuser,
					  -dbpass => $dbpass);
}

my $qseq_db = AT::DB::TranscriptSeq->connect( -dbname => $trseqdb_name,
					  -dbhost => $dbhost,
					  -dbuser => $dbuser,
					  -dbpass => $dbpass);

# Get alignments

my $af = AT::MapAlignmentFactory->new( compact => $LONG_ALN ? 0 : 1,
				       target_db => $target_db,
				       query_db => $query_db,
				       query_seq_db => $qseq_db,
				       upstr_context => 100,
				       downstr_context => 100 );

my @alignments = $af->get_alignments( accs => \@accs,
     				      combine => $COMBINE_ALN);

# Output alignments

my $outstream = Bio::AlignIO->new ( -fh => \*STDOUT,
				    -format => "clustalw" );
foreach my $aln (@alignments) {
    $outstream->write_aln ($aln);
    print "//\n";
}


sub usage
{
    warn "\n";
    warn "Usage: perl $0 [options] <database> <region | accessions>\n";
    warn "Example: perl $0 -c all AT_HS_JUL03 NM_004082 NM_023019\n";
    warn "Example: perl $0 -c all AT_HS_JUL03 chr1:128425129-128429399\n";
    warn "Options:\n";
    warn "  -l             long alignment; do not clip introns\n";
    warn "  -c all         combine all alignments (creates semi-MSA)\n";
    warn "  -c overlapping combine overlapping alignments (creates one or more semi-MSA)\n";
    warn "  -h <host>      database host (default: ",DEFAULT_HOST,")\n";
    warn "  -u <username>  database username (default: ",DEFAULT_USER,")\n";
    warn "  -p <password>  database password\n";
    warn "  -t <trseq_db>  transcript sequence db name (default: ",DEFAULT_TRSEQ_DB,")\n";
    warn "  -n <nibdir>    nibfile directory (default: use SQL assembly db)\n";
    warn "  -2             use experimental v.2 schema (this option will eventually be removed)\n";
    warn "\n";
    exit;
}


