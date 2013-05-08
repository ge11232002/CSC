#!/usr/bin/perl -w

# mappings2sql.pl: read mappings from psl/pslx file and store in database.
# See usage() function for usage details.

use strict;
use Getopt::Std;
use Term::ReadKey;
use AT::DB::GenomeMapping;
use AT::BlatpslxIO;

# Constant definitions
my $DEFAULT_DB_HOST = "localhost";

# Cmdline args
my %args;
getopts('bu:f:d:t:h:', \%args);
my $db_user = $args{'u'} or usage();
my $in_file = $args{'f'} or usage();
my $db_name = $args{'d'} or usage();
my $db_host = $args{'h'} || $DEFAULT_DB_HOST;
my $target_db = $args{'t'} or usage();
my $bin_field = $args{'b'};
my $qType = $ARGV[0];
unless($qType and ($qType eq 'mRNA' or $qType eq 'EST')) { usage(); }

# Password
my $db_pass = ask_password();

# Open input stream
my $mappings_instr = AT::BlatpslxIO->new( -file => $in_file,
					  -has_bin_field => $bin_field);

# Connect to database
my $db = AT::DB::GenomeMapping->connect( -dbname => $db_name,
					 -dbhost => $db_host,
					 -dbuser => $db_user,
					 -dbpass => $db_pass);
if ($db->num_mappings) {
    my $q = "Database $db_name already contains mappings. Proceed (yes/no)? ";
    exit 0 if (yesno($q) eq "no");
}

# Store mappings
my $importer = $db->create_mapping_importer;
my $num = 0;
while (my $mapping = $mappings_instr->next_mapping) {
    $mapping->qType($qType);
    $mapping->target_db($target_db);
    my $orig_bin = $mapping->bin;
    $importer->store_object($mapping);
    if(defined ($orig_bin) and $orig_bin != $mapping->bin) {
	warn "Bin ".$orig_bin." recalculated to ".$mapping->bin." for ".
	$mapping->tStart." - ".$mapping->tEnd;
    }
}
my @result = $importer->finish;

print join("+", @result), " records imported.\n";


# SUBROUTINES

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
    print "Usage: mappings2sql.pl ";
    print "-u <username> ";
    print "-f <input_file> -d <database> -t <target_db> [-h <host>] [-b] ";
    print "mRNA|EST\n";
    print "input_file\tMappings in psl/pslx format.\n";
    print "database\tName of mapping database to import into.\n";
    print "target_db\tName of genome database holding the target assembly.\n";
    print "host\t\tDatabase host (defaults to $DEFAULT_DB_HOST).\n";
    print "-b\t\tUse if input_file is a UCSC database dump with an extra\n";
    print "\t\tfirst field (bin field).\n";
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

