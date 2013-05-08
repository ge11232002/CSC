#include <string.h>

use AT::Tools::Run::Masker;
use Bio::SeqIO;
use Bio::Seq;
use Test;

plan(tests => 3);

my $file= "t/examples/shortseq";
my $seq1IO  = Bio::SeqIO -> new(-file => $file , '-format' => 'fasta');
my $seq1Obj = $seq1IO -> next_seq();
my $masker  = AT::Tools::Run::Masker -> new(-seq => $seq1Obj);
$masker ->run();
my $maskedseq = $masker ->get_masked_seq();

$maskedlength= $maskedseq->length();
$seqlenght   = $seq1Obj-> length();

ok($seqlenght, $maskedlength);


my $subseq1  = $maskedseq -> subseq(292,300);
my $subseq2   = $maskedseq -> subseq(1200,1203);

ok($subseq1, "NNNCAACAG");
ok($subseq2, 'NNNN');
