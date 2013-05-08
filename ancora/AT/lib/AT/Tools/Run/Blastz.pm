############################################################################
# Run::Blastz.pm -  module for running Blastz and 
#                   parsing the output
#
# Copyright
############################################################################
 
package AT::Tools::Run::Blastz;

=head1 NAME

Blastz

=head1 SYNOPSIS

    use Run::Blastz;
    my $runblastz = AT::Tools::Run::Blastz -> new(-seq1 => Bio::Seq object,
						  -seq2 => Bio::Seq object,
						  All BlastZ options:
						  -m    => 80M (bytes of space  for trace-back information),
						  -c    => 0 or 1 print census,
						  -v    => 0 or 1 verbose progress report to stderr,
						  -r    => 0 (old stle report) or 1 (new style report),
						  -B    => 0 (single strand) >0(both strands) ,
						  -C    => 0 (no chaining) 1 (just output chain) 2 (chain and extend),
						  -E    => 30 (gap extrension penalty),
						  -G    => 0 or higher (diagonal chaining penalty),
						  -H    => 0 or higher (prechain: use hsp, chain, retry strategy),
						  -K    => 3000 (threshold for MSPs),
						  -L    => K (Threshold for gapped alignments),
						  -M    => 50 (Mask threshold for seq1, if a bp is hit this many times),
						  -O    => 400 (gap-open penalty),
						  -P    => 1 (entropy not used) 1 (entropy used) >1 entropy with feedback
						  -Q    => (load the scoring matrix from a file, specify the filename),
						  -R    => 0 or higher (antidiagonal chaining penalty),
						  -T    => 0 (use W for wordsize) 1 (use 12of19),
						  -W    => 8 or other (word size),
						  -Y    => 0+300E (X-drop parameter for gapped extension));

    $runblastz -> run();
    $runblastz -> parse();

=head1 DESCRIPTION

B<Run::Blastz()> is a perl module for aligning two DNA sequences using 
the alignment tool Blastz.The class has two B<Bio::Seq> objects attributes.
The method B<Run::Blastz::Run()> executes BlastZ, and the method 
B<Run::Blastz::parse()> parses the output into simpleAlign objects. 
BlastZ tries to align the two provided sequences using the direct and 
reversed strand of sequence two.Up to two alignments are therefore 
produced in one execution of Blastz. Alignments can consist of several 
subalignments, which can consist of several ungapped blocks.

=head1 METHODS DESCRIPTION

=cut


############################################################################


use strict;

use AT::AlignmentSet;
use AT::Alignment::Composite;
use AT::Alignment::Subalignment;
use Bio::SeqIO;
use Bio::Seq;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Getopt::Long;
use File::Temp qw/ tempfile tempdir /;


=head2 B<new()>

I<Title:> new()

I<Usage:>

              $blastzrun = AT::Tools::Run::Blast->new(-seq1 => Bio::Seq object,
                                                      -seq2 => Bio::seq object,
                                                      -m    => $m,
						      -c    => $c,
						      -v    => $v,
						      -r    => $r,
						      -B    => $B,
						      -C    => $C,
						      -E    => $E,
						      -G    => $G,
						      -H    => $H,
						      -K    => $K,
						      -L    => $L,
						      -M    => $M,
						      -O    => $O,
						      -P    => $P,
						      -Q    => $Q,
						      -R    => $R,
						      -T    => $T,
						      -W    => $W,
						      -Y    => $Y);

I<input: > 

$seq1  Bio::Seq object to be aligned with -seq2

$seq2  the second Bio::seq object in the alignment

All possble BlastZ options, the letters stand for the same options as in the command-line call to BlastZ.


I<output:> Returns a new B<Run::Blastz> object.

=cut

sub new
  {
  my ($class, %args) = @_;
  my ($seq1,$seq2, $m,$c, $v, $r, $B, $C, $E, $G, $H, $K, $L,$M, $O, $P, $Q, $R, $T, $W, $Y);
  my $runcommand="";
  if (defined $args{-seq1}){
    $seq1=$args{-seq1};
    delete $args{-seq1};
  }
  if (defined $args{-seq2}){
    $seq2=$args{-seq2};
    delete $args{-seq2};
  }
  if (defined $args{ -m }){
    $runcommand=$runcommand." m=".$args{-m};
      delete $args{-m};
  }
  if (defined $args{-c}){
    $runcommand=$runcommand." c=".$args{-c};
    delete $args{-c};
  }
  if (defined $args{-v}){ 
    $runcommand=$runcommand." v=".$args{-v};
    delete $args{-v};
  }

  if (defined $args{-r}){
    $runcommand=$runcommand." r=".$args{-r};
    delete $args{-r};
  }
  if (defined $args{-B}){
    $runcommand=$runcommand." B=".$args{-B};
    delete $args{-B};
  }

  if (defined $args{-C}){
    $runcommand=$runcommand." C=".$args{-C};
    delete $args{-C};
  }
  if (defined $args{-E}){
    $runcommand=$runcommand." E=".$args{-E};
    delete $args{-E};
  }
  if (defined $args{-G}){
  $runcommand=$runcommand." G=".$args{-G};
    delete $args{-G};
  }
  if (defined $args{-H}){
    $runcommand=$runcommand." H=".$args{-H};
    delete $args{-H};
  }
  if (defined $args{-K}){
    $runcommand=$runcommand." K=".$args{-K};
    delete $args{-K};
  }
  if (defined $args{-L}){
    $runcommand=$runcommand." L=".$args{-L};
    delete $args{-L};
  }
  if (defined $args{-M}){
    $runcommand=$runcommand." M=".$args{-M};
    delete $args{-M};
  }
  if (defined $args{-O}){
    $runcommand=$runcommand." O=".$args{-O};
    delete $args{-O};
  }
if (defined $args{-P}){
  $runcommand=$runcommand." P=".$args{-P};
    delete $args{-P};
  }
  if (defined $args{-Q}){
    $runcommand=$runcommand." Q=".$args{-Q};
    delete $args{-Q};
  }
  if (defined $args{-R}){
    $runcommand=$runcommand." R=".$args{-R};
    delete $args{-R};
  }
  if (defined $args{-T}){
    $runcommand=$runcommand." T=".$args{-T};
    delete $args{-T};
  }if (defined $args{-W}){
    $runcommand=$runcommand." W=".$args{-W};
    delete $args{-W};
  }
  if (defined $args{-Y}){
    $runcommand=$runcommand." Y=".$args{-Y};
    delete $args{-Y};
  }



  return bless {
      _blastz_executable  => ($args{-executable} or "blastz"),
		_seq1     => $seq1,
		_seq2     => $seq2,
		_runcommand => $runcommand,
	       }, ref $class || $class;
}

=head2 B<get_seq1()>

I<Title:> get_seq1()

I<usage:> $seqobj1 = $blastzrun -> get_seq1();

I<output:> Returns the Bio::Seq object correspondning to sequence 1.

=cut

sub get_seq1          {$_[0] -> {_seq1} }


=head2 B<get_seq2()>

I<Title:> get_seq2()

I<Usage:> $seqobj2 = $blastzrun -> get_seq1();

I<output:> Returns the Bio::Seq object correspondning to sequence 2.

=cut

sub get_seq2          {$_[0] -> {_seq2} }


=head2 B<get_alignmentset()>

I<Title:> get_alignmentset()

I<Usage:> $alignmentset = $blastzrun -> get_alignmentset();

I<output:> Returns the alignmentset object created in B<Run::Blastz::run()>.

=cut

sub get_alignmentset  {$_[0] -> {_alignmentset}}


=head2 B<get_alignmentfile()>

I<Title:> get_alignmentfile

I<Usage> $alignfile = $blastzrun -> get_alignmentfile();

I<Output:> Returns the filename of the Blastz output.

=cut

sub get_alignmentfile {$_[0] -> {_alignmentfile}}

=head2 B<get_tempfile1()>

I<Title:> get_tempfile1()

I<Usage:> Used internally by B<Run::Blastz::parse()>.

I<Output:> Returns the filename of a temporary file holding sequence 1 in fasta format.

=cut

sub get_tempfile1 {$_[0] -> {_tempfile1}}

=head2 B<get_tempfile2()>

I<Title:> get_tempfile2()

I<Usage> Used internally by B<Run::Blastz::parse()>.

I<Output:> Returns the filename of a temporary file holding sequence 2 in fasta format.

=cut

sub get_tempfile2 {$_[0] -> {_tempfile2}}


=head2 B<run()>

I<Title:> run()

I<Usage:> $blastzrun -> run();

I<Output:> Nothing returned

I<Description:> Aligns two Bio::Seq objects using the program BlastZ. A call to this function sets the value for the _alignmentfile attribute of the B<Run::Blastz> object. 

=cut

sub run{
  my ($self,%args)=@_;
  if($args{-tmpdir}) {
    $self->{_tmpdir} = $args{-tmpdir};
    mkdir $self->{_tmpdir};
  }
  else {
    $self->{_tmpdir} = tempdir("BlastzXXXXXXXX", TMPDIR => 1, CLEANUP => 1);
  }
  $self->print_seq_to_file();
  my $seqfile1=$self->get_tempfile1();
  my $seqfile2=$self->get_tempfile2();
  my $alignmentfile = File::Temp::tempnam($self->{_tmpdir}, "align");
  my $runcommand= $self->{_runcommand};
  $runcommand = " ".$seqfile1." ".$seqfile2." ".$runcommand." > $alignmentfile";
  system($self->{'_blastz_executable'}.$runcommand);
  $self->{_alignmentfile}=$alignmentfile;
}

=head2 B<print_seq_to_file()>

I<Title:> print_seq_to_file()

I<Usage:> used internally by B<Run::Blastz::run()>

I<Output:> Nothing returned

I<Description:> Prints the two B<Bio::Seq> objects to be alignemd to two temporary files that can be read by Blastz in B<Run::Blastz::run>, and thereafter deleted in the destructor.

=cut

sub print_seq_to_file{
  my ($self, %args)=@_;
  my $tempfile1= File::Temp::tempnam($self->{_tmpdir}, "seq");
  my $tempfile2= File::Temp::tempnam($self->{_tmpdir}, "seq");
  my $output_stream1 = Bio::SeqIO->new(-file => ">$tempfile1", -format=>'Fasta');
  $output_stream1->write_seq($self->get_seq1());

  my $output_stream2 = Bio::SeqIO->new(-file => ">$tempfile2", -format=>'fasta');
  $output_stream2->write_seq($self->get_seq2());
  $self->{_tempfile1}=$tempfile1;
  $self->{_tempfile2}=$tempfile2;
}

=head2 B<parse()>

I<Title:> parse()

I<Usage:> $blastzrun -> parse();

I<Output:> Nothing returned

I<Description:> Parsing the alignmentfile created in B<Run::Blasrz::run()>.A call to B<Run::Blast::parse()> sets the Alignementset attribute of the B<Run::Blastz> object to point to an Alignmentset object.

=cut

sub parse{
  my ($self, %args)=@_;
  my $alignmentset = AT::AlignmentSet -> new (-seq1 => $self->get_seq1(), -seq2 => $self->get_seq2());
  open(ALIGNMENT, $self->get_alignmentfile());
  while(<ALIGNMENT>) {
    chomp($_);
    #for each strand of sequence two, if both are used in alignments...
    if ($_=~/^s \{/){
      #a new alignment object  is created:
      my $alignment = AT::Alignment::Composite->new(-seq1 => $self->get_seq1(), -seq2 => $self->get_seq2());
      my $firstseq=<ALIGNMENT>;
      my $secondseq=<ALIGNMENT>;
   	
      #explores the strand of seq2
      my $tempfile2=$self->get_tempfile2();

      if (index($secondseq, "$tempfile2-")> -1){ 
	$alignment->set_strand(-strand => "-");	
      }
      else{
	$alignment->set_strand(-strand => "+");
      }
     
      #Outer loop for looping through the alignment file row by row until
      #the strand of seq2 is changed (which is indicated by a #: on one row).
      until($_ =~ /^\#\:/){
	$_=<ALIGNMENT>;

	#This section regards a subalignment.
	my $subalignment;
	if($_=~/^a \{/){
	  my $line = <ALIGNMENT>;
	  my @ungapped_blocks;
	  chomp($line);
	  my $score;
	  while($line ne "}" ){
	    my @pos   = split(" ",$line);
	    $score = $pos[1] if $pos[0] eq "s";
	    # Start and end indexes of all ungapped blocks,
	    # for both sequences are stored in array.
	    push @ungapped_blocks, [@pos] if $pos[0] eq "l";
	    $line = <ALIGNMENT>;
	    chomp($line);
	  }
	  my @pos_seq1=($ungapped_blocks[0][1],$ungapped_blocks[scalar(@ungapped_blocks)-1][3]);
	  my @pos_seq2=($ungapped_blocks[0][2],$ungapped_blocks[scalar(@ungapped_blocks)-1][4]);


	  #Re-writes the aligned subsequences to get equal length.
	  my ($alignableseq1, $alignableseq2) = $self->make_alignable(-ub =>\@ungapped_blocks, -ps1 => \@pos_seq1, -ps2 => \@pos_seq2, -strand => $alignment->get_strand());


	  #produces a simple align object
	  my $simple_align_obj = $self->make_alignment(-ls1 => $alignableseq1, -ls2 => $alignableseq2, -ps1 => \@pos_seq1, -ps2 => \@pos_seq2, -strand => $alignment->get_strand());
	  my @lengths=($self->get_seq1()->length(),$self->get_seq2()->length());

	  #creates a new subalignment object.
	  $subalignment = AT::Alignment::Subalignment -> new (-ps1 => \@pos_seq1, -ps2 => \@pos_seq2, -score => $score, -ao => $simple_align_obj, -strand => $alignment->get_strand(), -lengths =>\@lengths, -ub => \@ungapped_blocks);
	  $alignment -> add_subalignment(-sub => $subalignment);
	}
      }
      $alignmentset->add_alignment(-alignment => $alignment);
    }
  }
  $self->{_alignmentset}=$alignmentset;
 # $self->delete_alignmentfile();
}

=head2 B<make_alignable()>


I<Title:> make_alignable()

I<Usage:> Called internally by B<Run::Blastz::parse()>

I<Input:> 

$ungapped_blocks - List in which each row contains the startpositions and endpositions in the corresponding B<Bio::Seq> objects of one ungapped block in a subalignment.

$pos_seq1 - Start and end positions in the first B<Bio::Seq> object of the entire subalignment.

$pos_seq2 - Start and end positions in the second B<Bio::Seq> object of the entire subalignment (If the second sequence was reverse complemented by blastz, the positions correspond to the reverse complemented sequence.)

I<Output:> 

$alseq1 - String containing the sequence corresponding to the aligned part of the first B<Bio::Seq> object, but in which gaps have been inserted at sertain positions in the sequence in between ungapped blocks.

$alseq1 - String containing the sequence corresponding to the aligned part of the second B<Bio::Seq> object, but in which gaps have been inserted at certain positions in the sequence in between ungapped blocks. If the second sequence was reverse complemented by blastz, $alseq1 will correspond to the reverse complementation.

I<Description:> $alseq1 and $alseq2 corresponds to substrings of the B<Bio::Seq> objects in the Blastz run. The sequences have been adjusted to get equal lengths, in order to be aligned using a Bio::Simplealign object. If the second sequence was reverse complemented in the alignment, the $alseq2 corresponds to the apropriate part of the I<reverse complemented> sequence.

=cut

sub make_alignable{
  my ($self, %args)=@_;
  my ($ungapped_blocks, $pos_seq1, $pos_seq2, $strand);
  if (defined $args{-ub}){
    $ungapped_blocks=$args{-ub};
    delete $args{-ub};
  }
  if (defined $args{-ps1}){
    $pos_seq1=$args{-ps1};
    delete $args{-ps1};
  }
  if (defined $args{-ps2}){
    $pos_seq2=$args{-ps2};
    delete $args{-ps2};
  }
  
  if (defined $args{-strand}){
    $strand=$args{-strand};
    delete $args{-strand};
  }
  #the original sequences:
  my $seq1=$self->get_seq1();
  my $seq2=$self->get_seq2();

  #If the second sequence was reverse complemented in the alignment by blastz, 
  #the positions above are calculaled on the reverse compemented sequence. The 
  #alignamble sequence should correspond to the reverse complementation, therefore
  #the reverse complementation of the original sequence will be used as a template.
  if($strand eq "-"){
    $seq2=$seq2->revcom;
  }
 

  #Creating two sequences of equal length that can be aligned using SimpleAlign
  #First ungapped block:
  my $alseq1=$seq1->subseq($ungapped_blocks->[0][1],$ungapped_blocks->[0][3]);
  my $alseq2=$seq2->subseq($ungapped_blocks->[0][2],$ungapped_blocks->[0][4]);
  for my $i(1  .. (scalar (@$ungapped_blocks)-1)){
    #adding  regions in between ungapped blocks:
    if ($ungapped_blocks->[$i][1]>($ungapped_blocks->[($i-1)][3]+1)){
      $alseq1=$alseq1.$seq1->subseq($ungapped_blocks->[($i-1)][3]+1,$ungapped_blocks->[$i][1]-1);
    }
    if ($ungapped_blocks->[$i][2]>($ungapped_blocks->[($i-1)][4]+1)){
      $alseq2=$alseq2.$seq2->subseq($ungapped_blocks->[($i-1)][4]+1,$ungapped_blocks->[$i][2]-1);
    }
    my $difference=(($ungapped_blocks->[($i)][2]-$ungapped_blocks->[($i-1)][4])-($ungapped_blocks->[$i][1]-$ungapped_blocks->[($i-1)][3]));
    if($difference>0){
      $alseq1=$alseq1.'-' x $difference;
      $alseq2=$alseq2;
    }
    elsif($difference<0){
      $alseq2=$alseq2.'-' x (-$difference);
      $alseq1=$alseq1;
    }
    #adds aligned regions
  if($ungapped_blocks->[$i][1]<=$ungapped_blocks->[$i][3]){
      $alseq1=$alseq1.$seq1->subseq($ungapped_blocks->[$i][1],$ungapped_blocks->[$i][3]);
    }
    #In case of strange (wrong?) alignments where the start position of an ungapped block is larger than the end position, the region between the two adjacent ungapped blocks is added.(This problem has ocurred for at least one alignment, perhaps t is a bug in BlastZ???...)
    else{
      die "Errors in the BlastZ alignment";
    }
    
    if($ungapped_blocks->[$i][2]<=$ungapped_blocks->[$i][4]){
      $alseq2=$alseq2.$seq2->subseq($ungapped_blocks->[$i][2],$ungapped_blocks->[$i][4]);
    }
    #In case of strange alignment where start position is larger than the end position
    #for an ungapped block.
    else{
      die "Errors in the BlastZ alignment";
    }
  }
  return($alseq1, $alseq2);
}


=head2 B<delete_alignmentfile()>

Removes the temporary file holding the result from Blastz.

=cut


sub delete_alignmentfile {
my ($self,%args) = @_;
my $file= $self->get_alignmentfile();
#system ("rm $file");
}


=head2 B<make_alignment()>

I<Title:> make_alignment

I<Usage:> Used internally by B<Run::Blastz::parse()>


I<Input:>

$locatableseq1 - String containing the first sequence to be included in the B<Bio::Simplealign> object.

$locatableseq1 - String containing the second sequence to be included in the B<Bio::Simplealign> object.

$pos_seq1 - Start and end positions in the first B<Bio::Seq> object of the entire subalignment.

$pos_seq2 - Start and end positions in the second B<Bio::Seq> object of the entire subalignment. The start and endpositions apply to the reverse complemented sequence if the strand is negative, which means that the positions in the original sequence are converted.

I<Output:>

$AlignObj - A B<Bio::SimpleAlign> object containg the two input sequences.

=cut



sub make_alignment{
  my ($self, %args) = @_;
  my ($locatableseq1, $locatableseq2, $pos_seq1, $pos_seq2);
  if (defined $args{-ls1}) {
    $locatableseq1=$args{-ls1};
    delete $args{-ls1};
  }
  if(defined $args{-ls2}){
    $locatableseq2=$args{-ls2};
    delete $args{-ls2};
  }
  if (defined $args{-ps1}) {
    $pos_seq1=$args{-ps1};
    delete $args{-ps1};
  }
  if(defined $args{-ps2}){
    $pos_seq2=$args{-ps2};
    delete $args{-ps2};
  }
 

  #LocatableSeq objects are used to make an alignment:
  $locatableseq1 = Bio::LocatableSeq -> new(-seq => $locatableseq1, -id => $self->get_seq1()->id(), -start =>$pos_seq1->[0], -end => $pos_seq1->[1]);
  $locatableseq2 = Bio::LocatableSeq -> new(-seq =>$locatableseq2, -id => $self->get_seq2()->id(), -start=>$pos_seq2->[0], -end =>$pos_seq2->[1]);


  
  #An alignment object is created:
  my $AlignObj= Bio::SimpleAlign->new();
  $AlignObj->add_seq($locatableseq1);
  $AlignObj->add_seq($locatableseq2);
  $AlignObj->gap_char('-');
  return ($AlignObj);
}

sub DESTROY  {
    my $self = shift;
    my $tmp = $self->{_tmpdir};
    print STDERR "Cleaning $tmp...\n";
    unlink <$tmp/*>;
    rmdir $tmp;
}
    

1;
