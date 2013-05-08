#!/usr/bin/env perl -w

use AT::MySQLdb;
use strict;
use Bio::AlignIO;

my $db = AT::MySQLdb->connect(-dbhost => "sql.cgb.ki.se",
			              -dbname => "AT",
                          -dbuser => "mihaela",
                          -dbpass => "baudolino2");

my $sth = $db->dbh->prepare("SELECT DISTINCT acc FROM QUERY_SEQUENCES");
$sth->execute;
my $i=0;
my $outstream = Bio::AlignIO->new(-fh=>\*STDOUT, -format=>"clustalw");
my $START_AFTER = ($ARGV[0] or 0);
my $after_start = !$START_AFTER;

while (my ($acc) = $sth->fetchrow_array) {
    #print STDERR ++$i." $acc\n";
    if ($acc eq $START_AFTER) {
	$after_start=1;
	next;
    }
    next unless $after_start;
    my @alignments = $db->get_alignments_for_acc({genome=>"MOUSE", compact=>1,
					      max=>1},
                                            $acc);
        if (scalar @alignments) {
		$outstream->write_aln($alignments[0]);
		print "//\n";
        }
}
$sth->finish;

exit(0);
