#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Std;
use MyPerlVars;
use CNE::CNE;
use CNE::DB;

my %args;
getopts('d:h:r', \%args);
my $DB_HOST = $args{'h'} || 'localhost';
my $DB_NAME = $args{'d'};
my $REMOVE_TABLES = $args{'r'};
my @IN_FILES = @ARGV;

unless($DB_NAME and @IN_FILES) {

    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Usage: 
    
$cmd <options> <input files>

Options:
-d <database name>     Required.
-h <database host>     Optional. Default: localhost.
-r                     Remove CNE table(s) if already exist

EndOfUsage

    exit;
}


# Connect to db
my $db = CNE::DB->connect(dbhost => $DB_HOST,
			  dbname => $DB_NAME,
			  dbuser => $MyPerlVars::sqlUser,
			  dbpass => $MyPerlVars::sqlPass)
    or die "could not connect to db $DB_NAME @ $DB_HOST";

# Load each input file
foreach my $fn (@IN_FILES) {

    # Parse filename to get: asm1, asm2, length, id
    my $filename_wo_dir = $fn;
    $filename_wo_dir =~ s/.*\///;
    my ($asm1, $asm2, $min_score, $win_size) = $filename_wo_dir =~ /^cne2w[A-Za-z0-9]*_(\w+)_(\w+)_(\d+)_(\d+)$/;
    unless($asm1 and $asm2 and $min_score and $win_size) {
	print STDERR "Error parsing file name $fn.\n";
	print STDERR "Filename should follow the format cne2w[Bf]_<asm1>_<asm2>_<minScore>_<winSize>, ";
	print STDERR "e.g. cne2wBf_hg18_danRer4_27_30\n";
	exit;
    }
    my $identity = int(1000 * $min_score / $win_size) / 1000;
    # Truncate identity to three decimals. Note that we truncate rather than round
    # so that the stated threshold never is higher than the actual threshold.

    # Create a database table
    my $table_name = $db->get_cne_table_name(assembly1 => $asm1,
					     assembly2 => $asm2,
					     min_length => $win_size,
					     min_identity => $identity,
					     version => 1);
    # ^version (here 1) should also be configurable

    if($db->cne_table_exists($table_name)) {
	if($REMOVE_TABLES) {
	    $db->drop_cne_table($table_name);
	}
	else {
	    die "Table $table_name for file $fn already exists. Use option -r to remove.\n";
	}
    }

    $db->create_cne_table($table_name);

    # Load data into table
    load_cne_data($db, $fn, $table_name, $asm1, $asm2);

    # Index table
    $db->index_cne_table($table_name);
}


sub load_cne_data
{
    my ($db, $fn, $table_name, $asm1, $asm2) = @_;

    open IN, $fn or die "could not open $fn";

    my $importer = $db->create_cne_importer($table_name);
    
    while(my $line = <IN>) {
	next if($line =~/^browser/ or $line =~ /^track/ or $line =~ /^\s*$/);
	chomp $line;
	my ($chr1, $start1, $end1, $chr2, $start2, $end2, $strand, $score, $cigar) = split /\t/, $line;
	my $cne = CNE::CNE->new(assembly1 => $asm1,
				assembly2 => $asm2,
				chr1 => $chr1, 
				chr2 => $chr2,
				start1 => $start1+1,
				start2 => $start2+1,
				end1 => $end1,
				end2 => $end2,
				strand => $strand,
				similarity => $score,
				cigar => $cigar);
	$importer->store_object($cne);
    }

    $importer->finish;

    close IN;
}



