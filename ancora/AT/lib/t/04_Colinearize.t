#include <string.h>


use AT::Alignment::Composite;
use AT::Alignment::Subalignment;
use AT::Tools::Run::Blastz;
use AT::AlignmentSet;
use AT::Tools::Colinearizer;
use Bio::SeqIO;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Test;
plan(tests => 4);

my $file1 = "t/examples/pai1_human";
my $seq1IO     = Bio::SeqIO -> new(-file => $file1 , '-format' => 'fasta');
my $seq1Obj = $seq1IO -> next_seq();
my $file2 = "t/examples/pai1_mus";
my $seq2IO     = Bio::SeqIO -> new(-file => $file2 , '-format' => 'fasta');
my $seq2Obj = $seq2IO -> next_seq();
my $runblast = AT::Tools::Run::Blastz ->new (-seq1 => $seq1Obj, -seq2 => $seq2Obj);
$runblast->run();
$runblast->parse();

my $alignmentset= $runblast->get_alignmentset();
my $alignments  = $alignmentset->get_alignments();
my $plusstrand  = $alignments->[0];
my $minusstrand = $alignments->[1];

#Colinearization:
my @colinearal;
if($plusstrand){
  my $P_colinearizer =AT::Tools::Colinearizer->new(-CA => $plusstrand);
  $P_colinearizer -> colinearize();
  $colinearal[0]=$P_colinearizer->get_colinear();
}
if($minusstrand){
  my $M_colinearizer = AT::Tools::Colinearizer->new(-CA => $minusstrand);
  $M_colinearizer -> colinearize();
  $colinearal[1]=$M_colinearizer->get_colinear();
}

my $colinear_alignmentset = AT::AlignmentSet -> new(-seq1 => seqObj1, -seq2 => seqObj2);
for $ca(@colinearal){
  $colinear_alignmentset->add_alignment(-alignment => $ca);
}

my $best_al=$colinear_alignmentset->get_best_alignment();

ok ($best_al->get_score(), 305624);

$best_al->set_top(-top => 0.3);
$best_al->set_tws(-tws => 50);
$best_al->set_threshold();
$best_al->set_cws(-cws => 50);

ok ($best_al->get_threshold(), 0.63);

$conserved_regions = $best_al->get_conserved_regions();
$posseq1=$conserved_regions->[7]->[0];
$posseq2=$conserved_regions->[7]->[2];

ok ($posseq1, 7299);
ok ($posseq2, 46189);

