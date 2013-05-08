#!/usr/bin/env perl -w

use AT::DB::GenomeMapping;
use AT::DB::GenomeAssembly;
use AT::MapAlignmentFactory;
use Bio::AlignIO;

my $mapdb = AT::DB::GenomeMapping->connect(-dbhost => "nautilus.cgb.ki.se",
					   -dbname => "AT_MM_FEB02",
					   -dbuser => "at_read",
					   -dbpass => "tittut");
my $genomedb = AT::DB::GenomeAssembly->connect(-dbhost => "nautilus.cgb.ki.se",
					       -dbname => "MM_FEB02",
					       -dbuser => "at_read",
					       -dbpass => "tittut");
my $af = AT::MapAlignmentFactory->new(query_db => $mapdb,
				   target_db=>$genomedb);

my @alignments = $af->get_alignments(acc=> "AF058956", compact => 1);
my $outstream = Bio::AlignIO->new(-fh=>\*STDOUT, -format=>"clustalw");
foreach my $aln (@alignments)  {
    $outstream->write_aln($aln);
}


foreach my $aln (@alignments) {
    foreach my $exon ($aln->target_exon_list)  {
	printf("Start: %10d  End: %10d\n", 
	       $exon->entire_seq->location_from_column($exon->start), 
	       $exon->entire_seq->location_from_column($exon->start));
    }
    print "------------------\n";
}

      
exit(0);
