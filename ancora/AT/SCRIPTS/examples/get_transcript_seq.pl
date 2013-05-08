#!/usr/bin/perl -w
#
# Script:  get_transcript_seq.pl
# Usage: perl get_transcript_seq.pl [options] [accession...]
#
# Give a list of accessions on the commandline or on standard input.
#
# Options:
#   -h  database host 
#   -d  database name
#   -u  database username
#   -p  database password
#

use strict;
use AT::DB::TranscriptSeq;
use Bio::SeqIO;
use Getopt::Std;

my %args;
getopts('h:u:p:d:', \%args);
my $DBNAME = $args{'d'} || 'TRANSCRIPT_SEQ';
my $DBHOST = $args{'h'} || 'nautilus.cgb.ki.se';
my $DBUSER = $args{'u'} || 'at_read';
my $DBPASS = $args{'p'} || ($args{'u'} ? '' : 'tittut');

my $out = Bio::SeqIO->new( -fh => \*STDOUT,
				 -format => "fasta");

my $seqdb = AT::DB::TranscriptSeq->connect( -dbname => $DBNAME,
					       -dbhost => $DBHOST,
					       -dbuser => $DBUSER,
					       -dbpass => $DBPASS);

if(@ARGV) {
    foreach my $id (@ARGV) {
	get_seq($id);
    }
}
else {
    warn "No accessions on commandline; reading accessions from STDIN...\n";
    while (my $line = <STDIN>) {
        chomp $line;
        get_seq($line);
    }
}

sub get_seq
{
    my $id = shift;
    my ($acc, $version) = $id =~ /^(\w+)\.(\d+)$/;
    $acc = $id unless($acc);    
    my $seq = $seqdb->get_seq($acc, $version);
    unless($seq) {
	warn "Sequence $id not found\n";
	return;
    }
    $seq->id(($seq->accession_number).'.'.($seq->version));
    $out->write_seq($seq);
}


