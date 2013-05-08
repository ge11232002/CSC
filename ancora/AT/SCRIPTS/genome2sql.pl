#!/usr/bin/perl -w

# genome2sql.pl: read chromosome seqs from fasta files and store in genome
# database.
# See usage() function for usage details.

use strict;
use Getopt::Std;
use Term::ReadKey;
use AT::DB::GenomeAssembly;

# Constant definitions
my $DEFAULT_DB_HOST = "localhost";

# Cmdline args
my %args;
getopt('udh:', \%args);
my $db_user = $args{'u'} or usage();
my $db_name = $args{'d'} or usage();
my $db_host = $args{'h'} || $DEFAULT_DB_HOST;
unless ($ARGV[0]) { usage(); }

# Password
my $db_pass = ask_password();

# Connect to database
my $db = AT::DB::GenomeAssembly->connect( -dbname => $db_name,
					 -dbhost => $db_host,
					 -dbuser => $db_user,
					 -dbpass => $db_pass);

# Import each chromosome file...
foreach my $filename (@ARGV) {
    $db->store_chromosome(file => $filename);
}

# SUBROUTINES

sub usage {
    print "Usage: genome2sql.pl ";
    print "-u <username> -d <database> [-h <host>] <input_file1> <input_file2> ...\n\n";
    print "database\tName of mapping database to import into.\n";
    print "input_file\tFasta-file containing the sequence of one chromosome.\n";
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
