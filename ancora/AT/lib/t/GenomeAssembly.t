#!/usr/bin/env perl -w

#use AT::DB::GenomeMapping;
use AT::DB::GenomeAssembly;
use Test;
plan(tests => 4);

my $db = AT::DB::GenomeAssembly->connect(-dbhost => "nautilus.cgb.ki.se",
					 -dbname => "MM_FEB02",
					 -dbuser => "at_read",
					 -dbpass => "tittut");

my $genomeseq = $db->get_genome_seq (
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

