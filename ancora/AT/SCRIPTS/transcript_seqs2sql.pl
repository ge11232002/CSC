#!/usr/bin/perl -W

# query_seqs2sql.pl: read sequences from fasta file and store in mapping
# database.
# See usage() function for usage details.

use strict;
use Getopt::Std;
use Term::ReadKey;
use Net::FTP;
use IO::Handle;
use IPC::Open2;
use AT::Seq::Generic;
use Bio::Seq::SeqFactory;
use Bio::SeqIO;
use AT::DB::TranscriptSeq;
use AT::DB::DataSource;

# FIX NEEDED:
# The daily-nc files may be mRNA or EST. This can be checked by looking at the division field.
# Currently we record them all as mRNA, which is wrong.
# Easy solution: ignore the daily files!

# Todo:
# * Versioning support and update flag
#   Fill in the field ds_version in the @FILES_WANTED array by querying the
#   db for the most recent version (maybe version should be changed to update
#   - for update nr).
#   If not in update mode: skip any datasources for which a version exists
#   If in update mode: increment the version for all data-sources by 1

# This array specifies which files to get.
# 'handler' should be a ref to a subroutine that does the following:
#   - check whether sequence should be included, otw return undef
#   - set type and species of sequence

my @FILES_WANTED = (

		    { ds_name => 'Genbank human',
		      ds_version => 2,
		      dir => '/genbank',
		      fnre => 'gbpri\d+\.seq\.gz',
		      'format' => 'genbank',
		      organism => ['Homo sapiens'],
		      seqtype => 'mRNA',
		      handler => \&gb_handler },
		    
		    { ds_name => 'RefSeq human',
		      ds_version => 2,
		      dir => '/refseq/H_sapiens/mRNA_Prot',
		      fnre => 'human\.rna\.gbff\.gz',
		      'format' => 'genbank',
		      organism => ['Homo sapiens'],
		      seqtype => 'mRNA',
		      handler => \&gb_handler },
		    
		    { ds_name => 'dbEST human',
		      ds_version => 3,
		      dir => '/blast/db/FASTA',
		      fnre => 'est_human\.gz',
		      'format' => 'fasta',
		      organism => 'Homo sapiens',
		      seqtype => 'EST',
		      division => 'EST',
		      handler => \&fa_handler },
		    
		    { ds_name => 'Genbank mouse',
		      ds_version => 2,
		      dir => '/genbank',
		      fnre => 'gbrod\d+\.seq\.gz',
		      'format' => 'genbank',
		      organism => ['Mus musculus'],
		      seqtype => 'mRNA',
		      handler => \&gb_handler },
		    
		    { ds_name => 'RefSeq mouse',
		      ds_version => 2,
		      dir => '/refseq/M_musculus/mRNA_Prot',
		      fnre => 'mouse\.rna\.gbff\.gz',
		      'format' => 'genbank',
		      organism => ['Mus musculus'],
		      seqtype => 'mRNA',
		      handler => \&gb_handler },
		    
		    { ds_name => 'dbEST mouse',
		      ds_version => 2,
		      dir => '/blast/db/FASTA',
		      fnre => 'est_mouse\.gz',
		      'format' => 'fasta',
		      organism => 'Mus musculus',
		      seqtype => 'EST',
		      division => 'EST',
		      handler => \&fa_handler },
		    
		    { ds_name => 'GenBank HTC',
		      ds_version => 3,
		      dir => '/genbank',
		      fnre => 'gbhtc\d+\.seq\.gz',
		      'format' => 'genbank',
		      organism => ['Homo sapiens', 'Mus musculus'],
		      seqtype => 'mRNA',
		      handler => \&gb_handler },
		    
		    { ds_name => 'GenBank Patent',
		      ds_version => 3,
		      dir => '/genbank',
		      fnre => 'gbpat\d+\.seq\.gz',
		      'format' => 'genbank',
		      organism => ['Homo sapiens', 'Mus musculus'],
		      seqtype => 'mRNA',
		      handler => \&gb_handler },
		    
		    { ds_name => 'GenBank Daily',
		      ds_version => 3,
		      dir => '/genbank/daily-nc',
		      fnre => 'nc\d+\.flat\.gz',
		      'format' => 'genbank',
		      organism => ['Homo sapiens', 'Mus musculus'],
		      seqtype => 'mRNA',
		      handler => \&gb_handler },
		    
		    );


# Oher 'constant' definitions
my $DEFAULT_DB_NAME = 'TRANSCRIPT_SEQ';
my $DEFAULT_DB_HOST = "localhost";
my $MAX_BATCH_SIZE = 100000;
my $FTP_URLBASE = 'ftp.ncbi.nih.gov';
my $FTP_USER = 'anonymous';

# Cmdline args
my %args;
getopts('u:d:h:e:p:w', \%args);
my $db_user = $args{'u'} or usage();
my $db_name = $args{'d'} || $DEFAULT_DB_NAME;;
my $db_host = $args{'h'} || $DEFAULT_DB_HOST;
my $db_pass = $args{'p'} || '?';
my $ftp_pass = $args{'e'} or usage();
my $use_wanted_list = $args{'w'};

my $db = AT::DB::TranscriptSeq->connect( -dbname => $db_name,
					 -dbhost => $db_host,
					 -dbuser => $db_user,
					 -dbpass => $db_pass);

if (!$use_wanted_list and $db->nr_seqs) {
    my $q = "Database $db_name already contains sequences and not in update mode. Proceed (yes/no)? ";
    exit 0 if (yesno($q) eq "no");
}

# Connect to ftp server
my $ftp = ftp_reconnect();

# Scan each file
foreach my $file (@FILES_WANTED) {
    unless($ftp->cwd($file->{dir})) {
	warn "Could not change to remote dir ".$file->dir;
	next;
    }
    my $fnre = $file->{fnre};
    foreach my $filename (grep { /^${fnre}$/ } $ftp->ls) {
	unless($db->dbh->ping) {
	    warn "transcrpt_seqs2sql: disconnected from database; trying to reconnect...\n";
	    unless($db->reconnect()) {
		die "transcript_seqs2sql: failed to reconnect to database\n";	
	    }
	    warn "transcript_seqs2sql: reconnected to database.\n";
	}
	my $data_source = $db->get_data_source
	    ($file->{'ds_name'}, $file->{'ds_version'}, $filename);
	if($data_source) {
	    warn( defined($data_source->load_end_time) ?
		  "File $filename already loaded; skipping\n" :
		  "Part of file $filename already loaded; skipping" );
	    next;
	}
	$data_source = AT::DB::DataSource->new
	    (name => $file->{'ds_name'},
	     version => $file->{'ds_version'},
	     component => $filename);
	warn "Getting ".$file->{dir}."/$filename\n";
	my ($unzip_out, $unzip_in) = (IO::Handle->new, IO::Handle->new);
	my $unzip_pid = open2($unzip_out, $unzip_in, 'gunzip');
	my $child_pid = fork;
	die "Couldn't fork" unless defined($child_pid);
	if($child_pid) {  # parent
	    close($unzip_out);
	    $ftp->get($filename, $unzip_in) or
		warn("Could not open remote file ".$file->{dir}."/$filename");
	    close($unzip_in);
	    waitpid($unzip_pid, 0);
	    waitpid($child_pid, 0);
	}
	else {		  # child
	    close($unzip_in);
	    process_file($unzip_out, $file, $data_source, $db, $use_wanted_list);
	    close($unzip_out);
	    exit;
	}
	$ftp = ftp_reconnect($ftp);  # make sure we are still connected
    }
}

$ftp->quit;
exit 0;


# SUBROUTINES

sub process_file {
    my ($in_fh, $file_info, $data_source, $db, $use_wanted_list) = @_;

    $db->add_data_source($data_source);

    # Create sequence stream on supplied filehandle
    my $seq_factory = Bio::Seq::SeqFactory->new(-type => 'AT::Seq::Generic');
    my $seq_in = Bio::SeqIO->new( -fh => $in_fh,
				  -format => $file_info->{'format'},
				  -seqfactory => $seq_factory
				  );
    $seq_in->sequence_builder->add_unwanted_slot('features','annotation');

    # Get each seq
    my $num_in_batch = 0;
    my $num_imported = 0;
    my $importer = $db->create_seq_importer;
    while(my $seq = $seq_in->next_seq) {
	#print $seq->accession_number, "\n";
	next unless($file_info->{'handler'}->($file_info, $seq));
	next if($use_wanted_list and
		!$db->is_wanted($seq->accession_number, $seq->version));
	$seq->data_source_id($data_source->data_source_id);
	if ($num_in_batch == $MAX_BATCH_SIZE) {
	    my ($num) = $importer->finish;
	    unless (defined $num) {
		warn "Error importing ".$data_source->id_str;
		return undef;
	    }
	    $num_imported += $num;
	    $importer = $db->create_seq_importer;
	    $num_in_batch = 0;
	    warn "new importer made\n";
	}
	$importer->store_object($seq);
	$num_in_batch++;
	#print "  OK\n";
    }
    my ($num) = $importer->finish;
    unless (defined $num) {
	warn "Error importing ".$data_source->id_str;
	return undef;
    }
    $num_imported += $num;

    $db->set_data_source_load_end_time($data_source);

    warn "$num_imported sequences imported\n";
}


sub gb_handler
{
    my ($file_info, $seq) = @_;
    my $organism = $seq->species->binomial;
    return undef unless (grep { $_ eq $organism } @{$file_info->{'organism'}}
		 and
		 $seq->molecule eq "mRNA");
    $seq->type($file_info->{'seqtype'});
    return 1;
}


sub fa_handler {
    my ($file_info, $seq) = @_;
    my ($acc, $ver) = 
	($seq->id =~ /^gi\|\d+\|\w+\|(\w+)\.(\d+)\|/);
    unless($acc) {
	warn ("Could not parse id ".$seq->id);
	return undef;
    }
    $seq->accession_number($acc);
    $seq->id('');
    $seq->version($ver);
    $seq->species(Bio::Species->new(-classification =>
		[ reverse split /\s+/, $file_info->{'organism'} ]));
    $seq->type($file_info->{'seqtype'});
    $seq->division($file_info->{'division'});
    return 1;
}


sub ftp_reconnect {
    my ($ftp) = @_;

    unless(defined($ftp) and $ftp->pwd) {
	$ftp = Net::FTP->new($FTP_URLBASE, Debug => 1) or
	    die "could not connect to ftp server";
	$ftp->login($FTP_USER, $ftp_pass) or
	    die "could not login to ftp server";
	$ftp->binary;
    }
	
    return $ftp;
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
    print "[-d <database>] [-h <host>] ";
    print "-u <username> -e <email> ";
    print "[-w]\n\n";
    print "email\t\tE-mail address to use for anonymous ftp.\n";
    print "database\tName of sequence database to import into.\n";
    print "host\t\tDatabase host (Defaults to localhost).\n";
    print "username\tUsername for sequence database.\n";
    print "-w\t\tUpdate mode. Add sequences named in WANTED table\n";
    exit 0;
}
