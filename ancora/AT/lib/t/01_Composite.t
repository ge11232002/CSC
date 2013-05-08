use AT::Alignment::Composite;
use AT::Alignment::Subalignment;
use AT::Tools::Run::Blastz;
#Dont know if all these are needed...
use Bio::SeqIO;
use Bio::Seq;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Bio::AlignIO;
use AT::Tools::Run::Primer3;
use AT::Tools::PrimerFinder;

use Test;
plan(tests => 5);

my $blastz_path;
eval {$blastz_path = `which blastz 2> /dev/null`;};
if (!$blastz_path or ($blastz_path =~ / /)) { # if space, then error message :)
    print "ok \# Skipped: (no blastz executable found)\n"x5;
    exit(0);
}


my $file1 = "t/examples/pai1_human";

my $seq1IO     = Bio::SeqIO -> new(-file => $file1 , '-format' => 'fasta');
my $seq1Obj = $seq1IO -> next_seq();

my $file2 = "t/examples/pai1_mus";
my $seq2IO     = Bio::SeqIO -> new(-file => $file2 , '-format' => 'fasta');
my $seq2Obj = $seq2IO -> next_seq();

my $runblast = AT::Tools::Run::Blastz ->new (-seq1 => $seq1Obj, -seq2 => $seq2Obj);


$runblast->run();

$runblast->parse();


my $alnset = $runblast->get_alignmentset(); 

$alnset->set_tws(-tws => 50);
$alnset->set_cws(-cws => 50);
$alnset->set_top(-top => 0.3);
$alnset->set_threshold();

ok ($alnset->get_best_alignment()->get_score(), 546192);


my @target=(($alnset->get_best_alignment()->get_subalignments()->[9]->get_conserved_regions()->[0][0]),($alnset->get_best_alignment()->get_subalignments()->[9]->get_conserved_regions->[0][1]));

ok($target[0], 17889);
ok($target[1], 17948);

my $target=\@target;
my $pcrl="400-2000";
my $xp=300;
my $xs=50;
my $pcr_nret=10;

my $primer3_path;
eval {$primer3_path = `which primer3_core 2> /dev/null`;};
if (!$primer3_path or ($primer3_path =~ / /)) { # if space, then error message :)
    print "ok \# Skipped: (no primer3_core executable found)\n"x2;
    exit(0);
}


my $primerfinder = AT::Tools::PrimerFinder -> new (-seq => $seq1Obj,  
						   -pcrl => $pcrl, 
						   -sl => $pcrl, 
						   -xp => $xp, 
						   -xs => $xs, 
						   -pcr_nret => $pcr_nret);

ok($primerfinder->get_xp(),300);


my @pcr = $primerfinder->pcr_primers(-target => $target);
my @seqprimers;
my @seqcounter;

for $i(0 .. scalar(@pcr)-1){
  @seqpp=$primerfinder->seq_primers(-target => $target, -pp => $pcr[$i]);
  $seqcounter[$i]=scalar(@seqpp);	
  push @seqprimers, \@seqpp;
}

ok(scalar(@pcr), 3);

