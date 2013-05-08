#!/usr/bin/perl -w
use strict;

use AT::BlatpslxIO;
use AT::DB::GenomeAssemblyNibs;

unless(@ARGV == 1) {
    print "This script select the best mappings for each query, based on the\n";
    print "score AT::Mapping::score_A1. If for one query there are several mappings\n";
    print "with the best score, the script will output all of them.\n";
    print "Usage:\n";
    print "$1 <nibDir>\n";
    print "Mappings (psl-format) are read from STDIN and written to STDOUT.\n";
    exit;
}

my $assembly_dir = $ARGV[0];

my $gendb = AT::DB::GenomeAssemblyNibs->new(dir => $assembly_dir);

my $map_in = AT::BlatpslxIO->new(-fh => \*STDIN);
my $map_out = AT::BlatpslxIO->new(-fh => \*STDOUT);

my (@best, $best_score);
while(my $m = $map_in->next_mapping) {
    $m->set_hsp_flanks($gendb);
    my $score = $m->score_A1;
    if(@best and $best[0]->qName eq $m->qName) {
	if($score > $best_score) {
	    @best = ($m);
	    $best_score = $score;
	}
	elsif($score == $best_score) {
	    push @best, $m;
	}
    }
    else {
	foreach my $b (@best) { $map_out->write_mapping($b); }
	@best = ($m);
	$best_score = $score;
    }
}
foreach my $b (@best) { $map_out->write_mapping($b); }

# this method is currently not used
sub pcid {
    my ($m) = @_;
    return ($m->matches + $m->repMatches) /
	($m->matches + $m->repMatches + $m->misMatches + $m->qNumInsert);
}
