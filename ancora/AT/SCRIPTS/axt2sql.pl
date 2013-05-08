#!/usr/bin/perl -w

# axt2sql.pl: read wga in axt format and store in database.
# See usage() function for usage details.

use strict;
use warnings;
use Getopt::Std;
use Term::ReadKey;
use Net::FTP;
use IO::Handle;
use IPC::Open2;
use AT::AxtIO;
use AT::DB::GenomeAssemblyTwoBit;
use AT::DB::GenomeAssemblyNibs;
use AT::DB::GenomeAlignment;
use Bio::AlignIO;

# Constant definitions
my $DEFAULT_DB_HOST = "localhost";
my $DEFAULT_ALN_TYPE = "net";
my $FTP_USER = "anonymous";

# Global variable to keep track of child processes
my @CHILD_PIDS;

# Cmdline args
my %args;
getopts('d:e:h:u:t:q:n:a:i', \%args);
my $db_name = $args{'d'};
my $db_host = $args{'h'} || $DEFAULT_DB_HOST;
my $db_user = $args{'u'} or usage();
my $genome1 = $args{'t'} or usage();
my $genome2 = $args{'q'} or usage();
my $aln_type = $args{'a'} || $DEFAULT_ALN_TYPE;
my $genome2_asm_spec = $args{'n'} or usage();
my $create_indices = $args{'i'};
my $ftp_pass = $args{'e'};
my @in_files = @ARGV;
usage() unless(@in_files or $create_indices);
$db_name = "AT_$genome1" unless ($db_name);

# Connect to databases
my $genome_db2;
if(-d $genome2_asm_spec) {
    $genome_db2 = AT::DB::GenomeAssemblyNibs->new(id => $genome2,
						  dir => $genome2_asm_spec)
	or die "Could not find nib files at $genome2_asm_spec";
}
else {
    $genome_db2 = AT::DB::GenomeAssemblyTwoBit->new(id => $genome2,
						    file => $genome2_asm_spec);
}

my $db = AT::DB::GenomeAlignment->connect( -dbname => $db_name,
				           -dbhost => $db_host,
				           -dbuser => $db_user,
				           -dbpass => '??')
    or die "Failed to connect to database $db_name";

# Import files if any given

if(@in_files) {   
    if($in_files[0] =~ /^ftp:\/\// or $in_files[0] =~ /^http:\/\//) {
	usage() if(@in_files > 1);
	unless($ftp_pass) {
	    print "Enter e-mail address to use for anomymous FTP: ";
	    $ftp_pass = <STDIN>;
	    chomp $ftp_pass;
	}
	import_files_from_ftp($db, $genome_db2, $aln_type, $in_files[0], $FTP_USER, $ftp_pass);
    }
    else {
	import_files_from_disk($db, $genome_db2, $aln_type, \@in_files);
    }
}

# Create indices if requested

if($create_indices) {
    # Reconnect in case connection has been lost
    # (any child processes will have killed the connection upon termination)
    unless($db->ping) {
	$db->reconnect() or die "Failed to reconnect to database";
    }
    # Index database
    print STDERR "Indexing database...\n";
    $db->index_genome_alignment($genome_db2->id, $aln_type);
}

print STDERR "Done.\n";


#############
# SUBROUTINES
#############

sub import_files_from_ftp
{
    my($db, $genome_db2, $aln_type, $url, $ftp_user, $ftp_pass) = @_;

    my (undef,$host,$dir) = $url =~ /^(ftp|http):\/\/([\w\.]+)(.+)/;
    die "could not parse $url" unless($dir);

    # Pipe and fork
    pipe AXT_R, AXT_W;
    my $child_pid = fork;
    die "Couldn't fork" unless defined($child_pid);
    if($child_pid) {  # parent
	close(AXT_R);
	push @CHILD_PIDS, $child_pid;
	# Connect to ftp server
	my @ftp_params = ($host, $ftp_user, $ftp_pass, $dir);
	my $ftp = ftp_reconnect(undef, @ftp_params);
	# Get each file, pipe it through gunzip, and on to child process
	foreach my $filename (grep { /\.axt\.gz$/ } $ftp->ls) {
	    print STDERR "Fetching and parsing $filename...\n";
	    my $unzip_in = ""; #IO::Handle->new;
	    my $unzip_pid = open2('>&AXT_W', $unzip_in, 'gunzip');
	    $ftp = ftp_reconnect($ftp, @ftp_params);  # make sure we are still connected
	    $ftp->get($filename, $unzip_in) or
		kill_and_die("Could not open remote file $filename");
	    close($unzip_in);
	    waitpid($unzip_pid, 0);
	}
	close(AXT_W);
	$ftp->quit;   
	waitpid($child_pid, 0);
    }
    else {	     # child
	close(AXT_W);
	# Create importer object and use it to import all axt data
	# Note that we create the importer in the child because having it
	# in both forks would cause its DESTROY method to be called twice
	my $importer = $db->create_alignment_importer($genome_db2, $aln_type);
	import_file(\*AXT_R, $importer);
	close(AXT_R);
	finish_import($importer);
	exit;
    }

}


sub import_files_from_disk
{
    my ($db, $genome_db2, $aln_type, $filenames) = @_;
    # Create importer object and use it to import all axt data
    my $importer = $db->create_alignment_importer($genome_db2, $aln_type);
    foreach my $fn (@$filenames) {
	open IN, $fn or die "Could not open $fn";
	print STDERR "Parsing $fn...\n";
	import_file(\*IN, $importer);
    }
    finish_import($importer);
}


sub import_file
{
    my ($fh, $importer) = @_;
    my $axt_in = AT::AxtIO->new(-fh => $fh);
    my $records = 0;
    while(my $aln = $axt_in->next_aln) {
	$importer->store_object($aln);
	$records++;
    }
    return $records;
}


sub finish_import
{
    my $importer = shift;
    print STDERR "Importing data into database...\n";
    my $nr_rows = $importer->finish;
    die("Failed to import data into database") unless defined($nr_rows);
}


sub ftp_reconnect {
    my ($ftp, $host, $user, $pass, $dir) = @_;

    unless(defined($ftp) and $ftp->pwd) {
	$ftp = Net::FTP->new($host, Debug => 0)
	    or kill_and_die("could not connect to ftp server $host");
	$ftp->login($user, $pass)
	    or kill_and_die("could not login to ftp server");
	if($dir) {
	    $ftp->cwd($dir)
		or kill_and_die("Could not change to remote dir $dir");
	}
	$ftp->binary;
	print STDERR "Logged in to $host", ($dir ? ", directory $dir" : ""), "\n";
    }
	
    return $ftp;
}


sub kill_and_die {
    kill(9, @CHILD_PIDS);
    foreach my $pid (@CHILD_PIDS) {
	waitpid($pid,0);
    }
    die(@_);
}


sub usage {
    print <<EndOfUsage;

Usage: perl $0 [options] [input_files]

Options:
-d <database>        Name of alignment database to import into
                     (default: AT_<target_genome>)
-h <host>            Database host (default: $DEFAULT_DB_HOST)
-u <username>        Database username (required)
-t <target_genome>   Name of first assembly in alignments
                     (required; e.g. HS_JUL03)
-q <query_genome>    Name of second assembly in alignments
                     (required; e.g. MM_FEB03)
-n <query_assembly>  Name of .2bit-file or .nib-file directory
                     for query genome (required)
-a <aln_type>        Type of alignments (default: $DEFAULT_ALN_TYPE)
-e <e-mail>          E-mail address for anonymous FTP
                     (required if files are to be downloaded)
-i                   Index database after loading any input files
                     (default: no indexing)

Input files argument can be
* Filenames of one or more axt files to be imported to the database
 or
* A FTP or HTTP url specifying a directory from which all files ending in .axt.gz
will be downloaded by anonymous FTP and imported to the database, e.g.
http://hgdownload.cse.ucsc.edu/goldenPath/hg17/vsMm7/axtNet/

EndOfUsage
    exit 0;
}




