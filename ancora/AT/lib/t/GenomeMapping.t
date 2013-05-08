#!/usr/bin/env perl -w

use AT::DB::GenomeMapping;
use Test;
plan(tests => 3);

my $db = AT::DB::GenomeMapping->connect(-dbhost => "nautilus.cgb.ki.se",
					-dbname => "AT_MM_FEB02",
					-dbuser => "at_read",
					-dbpass => "tittut");
my $acc = "0610007F07";
my ($mapping) = $db->get_mappings_for_acc ($acc);
ok (ref($mapping) and $mapping->isa("AT::Mapping"));
my $seq = $db->get_query_seq($acc);

if (ref($seq) and $seq->isa("Bio::SeqI"))  {
    ok(1);
    ok($seq->length, 1294);
}
else { ok(0); }

my @mappings =  $db->get_mappings_in_region(-chr=>"chr3",
					    -target_db => "MM_FEB02",
					     -start=>108921750,
					     -end => 108943988);
foreach my $mp (@mappings) {
    print join("\t", $mp->qName, $mp->matches, $mp->misMatches, 
	       $mp->tStart."-".$mp->tEnd."[".$mp->strand."]")."\n";
}
