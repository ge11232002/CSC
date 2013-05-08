#!/usr/bin/perl -W

# query_seqs2sql.pl: read sequences from fasta file and store in mapping
# database.
# See usage() function for usage details.

use strict;
use Getopt::Std;
use Bio::SeqIO;
use Term::ReadKey;
use AT::DB::GenomeMapping;

# Constant definitions
my $DEFAULT_DB_HOST = "localhost";
my $MAX_BATCH_SIZE = 100000;

# Cmdline args
my %args;
getopts('eu:f:d:h:', \%args);
my $db_user = $args{'u'} or usage();
my $in_file = $args{'f'} or usage();
my $db_name = $args{'d'} or usage();
my $db_host = $args{'h'} || $DEFAULT_DB_HOST;
my $seqs_are_ests = $args{'e'};

# Password
my $db_pass = ask_password();

# Open input stream
my $seq_instream = Bio::SeqIO->new( -file => $in_file,
				    -format => 'fasta');

# Connect to database
my $db = AT::DB::GenomeMapping->connect( -dbname => $db_name,
					 -dbhost => $db_host,
					 -dbuser => $db_user,
					 -dbpass => $db_pass);
if ($db->num_query_seqs) {
    my $q = "Database $db_name already contains query sequences. Proceed (yes/no)? ";
    exit 0 if (yesno($q) eq "no");
}

# Store each sequence
my $num_in_batch = 0;
my $num_imported = 0;
my $importer = $db->create_query_seq_importer;
while (my $seq = $seq_instream->next_seq) {
    if ($num_in_batch == $MAX_BATCH_SIZE) {
	my ($num) = $importer->finish;
	$num_imported += $num;
	$importer = $db->create_query_seq_importer;
	$num_in_batch = 0;
	warn "new importer made\n";
    }
    set_acc_ver($seq, $seqs_are_ests);
    $importer->store_object($seq);
    $num_in_batch++;
}
my ($num) = $importer->finish;
$num_imported += $num;

print "$num_imported sequences imported.\n";


# SUBROUTINES

sub set_acc_ver {
    my ($seq, $is_est) = @_;
    if($is_est) {
	my ($acc1, $ver, $acc2) = 
	    ($seq->id =~ /^gi\|\d+\|\w+\|(\w+)\.(\d+)\|(\w+)/);
	unless($acc2) {
	    warn ("Could not parse id ".$seq->id);
	    return;
	}
	warn "Differing identifiers: $acc1 $acc2" unless ($acc1 eq $acc2);
	$seq->id($acc1);
	$seq->version($ver);
    }
    else {
	$seq->version(0);
    }
}


sub yesno {
    my $question = shift;
    my $answer;
    do {
	print $question;
	$answer = <STDIN>;
	chomp $answer;
	$answer = lc $answer;
    } until ($answer eq "yes" or $answer eq "no");
    return $answer;
}


sub usage {
    print "Usage: query_seqs2sql.pl ";
    print "-u <username> -f <input_file> -d <database> [-h <host>] [-e]\n\n";
    print "input_file\tSequences to import in fasta format.\n";
    print "database\tName of mapping database to import into.\n";
    print "host\t\tDatabase host (Defaults to localhost).\n";
    print "-e\t\tUse if sequences are Genbank ESTs.\n";
    exit 0;
}

sub ask_password {
    ReadMode('noecho');
    print "Password: ";
    my $pw = ReadLine(0);
    ReadMode('restore');
    print "\n";
    chomp($pw);
    return $pw;
}
