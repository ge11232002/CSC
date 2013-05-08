#!/usr/bin/env perl -w

use AT::MySQLdb;
use Test;
plan(tests => 4);

my $db = AT::MySQLdb->connect(-dbhost => "sql.cgb.ki.se",
			              -dbname => "AT",
                          -dbuser => "mihaela",
                          -dbpass => "baudolino2");

my $genomeseq = $db->get_genome_seq (genome=>"mouse",
                                     chr=>1,
                                     start =>4222301,
                                     end => 4282300);
if (ref($genomeseq) and $genomeseq->isa("Bio::LocatableSeq"))  {
    ok(1);
    ok($genomeseq->length, 60000);
    ok($genomeseq->start, 4222301);
    ok($genomeseq->end, 4282300);
}
else { ok(0); }

