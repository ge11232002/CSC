

use Bio::SeqIO; #OK

use AT::Tools::PrimerFinder;
use Test;
plan(tests=>4);

my $file1 = "t/examples/pai1_human";
my $seq1IO     = Bio::SeqIO -> new(-file => $file1 , '-format' => 'fasta');
my $seq1Obj = $seq1IO -> next_seq();

#print "Length of input sequence: ".$seq1Obj->length."\n";
ok ($seq1Obj->length, 48682);

#Setting prameters for the primerfinder:
my $pcrl="400-1200";
my $xp=300;
my $xs=50;
my $pcr_nret=5;

#creating a primerfinder object
my $primerfinder = AT::Tools::PrimerFinder -> new (-seq => $seq1Obj,  -pcrl => $pcrl, -sl => $pcrl, -xp => $xp, -xs => $xs, -pcr_nret => $pcr_nret);

#Speccifying target regions:
my @target=([5000,5100],
	    [10000, 10100],
	    [15000, 15100],
	   );
my @nr_products = (4,3,3);

for $t(@target){
  my $target=$t;
  # print "Target: ".$target->[0]."-".$target->[1]."\n";
  #Searching for PCR primers for this region:
  my @pcr = $primerfinder->pcr_primers(-target => $target);
  ok (scalar(@pcr), shift @nr_products);
  my @seqprimers;
  my @seqcounter;
  
  #Looking for sequencing primers for each PCR primer pair:
  for $i(0 .. scalar(@pcr)-1){
    @seqpp=$primerfinder->seq_primers(-target => $target, -pp => $pcr[$i]);
    $seqcounter[$i]=scalar(@seqpp);	
    push @seqprimers, \@seqpp;
  }
  
}


