#!/usr/bin/env perl -w

use Test;
use AT::BlatpslxIO;
plan(tests => 2);
my $stream =  AT::BlatpslxIO->new(-file=>"t/test.pslx");
my $mapping = $stream->next_mapping;
ok(ref $mapping, "AT::Mapping");
ok(!scalar($mapping->all_HSPs), 1);
while ($mapping=$stream->next_mapping)  {
    print $mapping->qStart." ".$mapping->qEnd."\n";
}
    
    
    
    
