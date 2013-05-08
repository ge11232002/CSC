############################################################################
# Alignmentset.pm - perl module for handeling a set of composite alignments 
#
# Copyright
############################################################################
 

=head1 NAME

Alignmentset

=head1 SYNOPSIS

    use Alignmentset;
    my $alignmenset = Alignmentset -> new(-seq1 => Bio::Seq object, 
                                          -seq2 => Bio::Seq object);
    $alignmentset  -> add_alignment(-alignment => $alignment); 
    $alignmentset  -> set_thresholds(-th = $threshold);
    $alignmentset  -> set_top(-top = $top);
    $alignmentset  -> set_tws(-tws = $tws);
    $alignmentset  -> set_cws(-cws = $cws);
    $bestalignment = $alignmentset->get_best_alignment();
    $alignments = $alignmentset->get_alignments();
    $seq1 = $alignmentset->get_seq1();
    $seq2 = $alignmentset->get_seq2();



=head1 DESCRIPTION

B<AlignmentSet> is a perl module for handeling a set of composite alignments. An B<AlignmentSet> object holds a list of composite alignments. New composite alignments can be added to the list. The composite alignments are automatically sorted according to their score values when new alignments are added, and the best alignment can be extracted. A dynamic threshold for similarity can be calculated for all composite alignments of a B<AlignmentSet> object, and the module can be used for finding conserved regions in the composite alignments.

=head1 METHODS DESCRIPTION

=cut

############################################################################


package AT::AlignmentSet;

use strict;
use Bio::SeqIO;
use Bio::Seq;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Getopt::Long;


=head2 B<new()>

I<Title:> new()

I<Usage:> $alignmentset = AT::AlignmentSet->new(-seq1 => Bio::Seq object,
                                                -seq2 => Bio::Seq object);

I<Input:>

$seq1 = the first Bio::Seq used in the alignment

$seq1 = the second Bio::Seq used in the alignment


I<Output:> Returns a new B<Alignmentset> object. 

I<Description:> 

=cut

sub new
{
  my ($class, %args) = @_;
  my ($seq1, $seq2);
  if (defined $args{-seq1}){
    $seq1=$args{-seq1};
    delete $args{-seq1};
  }
  if (defined $args{-seq2}){
    $seq2=$args{-seq2};
    delete $args{-seq2};
  }

  return bless {
		_seq1          => $seq1,
		_seq2          => $seq2,
	
	       }, ref $class || $class;
}


=head2 B<get_Seq1()>

I<Title:> get_seq1()

I<Usage:> seq1 = $alignmentset -> get_seq1();

I<Input:> 

I<Output:> Returns a B<Bio::Seq> object 

=cut

sub get_seq1              {$_[0] -> {_seq1} }

=head2 B<get_Seq2()>

I<Title:> get_seq2()

I<Usage:> seq2 = $alignmentset -> get_seq2();

I<Input:> 

I<Output:> Returns a B<Bio::Seq> object 

=cut

sub get_seq2              {$_[0] -> {_seq2} }

=head2 B<get_alignments()>

I<Title:> get_alignments()

I<Usage:> $alignments = $alignmentset -> get_alignments();

I<Input:> 

I<Output:> Returns a reference to a list conatining all alignments in the alignmentset.

=cut

sub get_alignments        {$_[0] -> {_alignments}}

=head2 B<get_sorted_alignments()>

I<Title:> get_sorted_alignments()

I<Usage:> $sorted_alignments = $alignmentset -> get_sorted_alignments();

I<Input:> 

I<Output:> Returns a reference to a list conatining all alignments in the alignmentset sorted according to score.

=cut

sub get_sorted_alignments {$_[0] -> {_sorted_alignments}}

=head2 B<get_best_alignment()>

I<Title:> get_best_alignment()

I<Usage:> $best_alignment = $alignmentset -> get_best_alignment();

I<Input:> No inputs

I<Output:> A reference to a B<Composite_alignment> object

I<Description:> The first alignment in the list of sorted alignments is returned.

=cut

sub get_best_alignment{
  my ($self, %args) = @_;
  return $self->get_sorted_alignments()->[0];
}


=head2 B<add_alignment()>

I<Title:> add_alignment()

I<Usage:> $alignmentset -> add_alignment(-alignment => $composite_alignment);

I<Input:> 

$composite_alignment is a reference to a B<Composite_Alignment> object.

I<Output:> No output

I<Description:> The input B<Composite_Alignment> object is added to the list of alignments in the alignmentset.

=cut

sub add_alignment{
  my ($self, %args)=@_;
  my $alignment;
  if (defined $args{-alignment}){
    $alignment=$args{-alignment};
    delete $args{-alignment};
  }
  my @temp_al;
  if($self->get_alignments()){
    my $alref=$self->get_alignments();
    @temp_al=@$alref;
  }
  push @temp_al, $alignment;
  $self->{_alignments}=\@temp_al;
  $self->sort_by_score();
}


=head2 B<set_threshold()>

I<Title:> set_threshold()

I<Usage:> $alignmentset -> set_thershold(-th => $threshold);

I<Input:> 

$threshold = Integer representing the manually chosen threshold (optional).

I<Output:> No output

I<Description:>The threshold is set automatically for every subalignment of every composite alignment in the alignmentset. If the argument -th was specified, the threshold will be set to this manually chosen value. If -th is not specified, the thresholds will be calculated according to the values of $top in the composite alignments, and set to this walue in all subalignments respectively.

=cut

sub set_threshold{
  my ($self, %args) = @_;
  my $th;
  my $alignments=$self->get_alignments();
  if (defined $args{-th}){
    $th = $args{-th};
    foreach my $alignment(@$alignments){
      $alignment->set_threshold(-th => $th);
    }
    delete$args{-th};
  }
  else{
    foreach my $alignment(@$alignments){
      $alignment->set_threshold();
    }
  }
}

=head2 B<set_top()>

I<Title:> set_top()

I<Usage:> $alignmentset -> set_top(-top => $top);

I<Input:> 

$top = Integer representing a value of top (stringency for conservation) (Optional).

I<Output:> No output

I<Dercription:> Sets the value of -top for all composite alignments of the alignmentset.

=cut

sub set_top{
  my ($self, %args) = @_;
  my ($top);
  my $alignments=$self->get_alignments();
  if (defined $args{-top}){
    $top= $args{-top};
    my $alignments=$self->get_alignments();
    foreach my $alignment(@$alignments){
      $alignment->set_top(-top => $top);
    }
    delete$args{-top};
  }
}

=head2 B<set_tws()>

I<Title:> set_tws()

I<Usage:> $alignmentset -> set_tws(-tws => $tws);

I<Input:> 

$tws = Integer representing a value of top (stringency for conservation) (Optional).

I<Output:> No output

I<Description:> Sets the value of -tws for all composite alignments (and their subalignments) of the alignmentset.

=cut

sub set_tws{
  my ($self, %args) = @_;
  my ($tws);
  my $alignments=$self->get_alignments();
  if (defined $args{-tws}){
    $tws= $args{-tws};
    my $alignments=$self->get_alignments();
    foreach my $alignment(@$alignments){
      $alignment->set_tws(-tws => $tws);
    }
    delete$args{-tws};
  }
}

=head2 B<set_cws()>

I<Title:> set_cws()

I<Usage:> $alignmentset -> set_cws(-tws => $cws);

I<Input:> 

$cws = Integer representing a value of top (stringency for conservation) (Optional).

I<Output:> No output

I<Description:> Sets the value of -cws for all composite alignments (and their subalignments) of the alignmentset.

=cut

sub set_cws{
  my ($self, %args) = @_;
  my ($cws);
  my $alignments=$self->get_alignments();
  if (defined $args{-cws}){
    $cws= $args{-cws};
    my $alignments=$self->get_alignments();
    foreach my $alignment(@$alignments){
      $alignment->set_cws(-cws => $cws);
    }
    delete$args{-cws};
  }
}





=head2 B<sort_by_score()>

I<Title:> sort_by__score()

I<Usage:> Used internally when new alignments are added.

I<Input:> No inputs

I<Output:> No output

I<Description:> Sorting the list of composite alignments in the alignmentset according to decreasing score values.A call to this functions sets a reference to the sorted_alignments list of the alignmentset.

=cut

sub sort_by_score{
  my ($self, %args) = @_;
  my @sorted_alignments;
  my $alignments=$self->get_alignments();
  my @alignments=@$alignments;
  @sorted_alignments= sort by_score @alignments;
  $self->{_sorted_alignments}=\@sorted_alignments;
}


=head2 B<by_score()>

I<Title:> by__score()

I<Usage:> used internally by B<AlignmentSet::sort_by_score>

I<Input:> No inputs

I<Output:> No output

I<Description:> Function to allow sorting according to score.

=cut

sub by_score{
$b->get_score() <=> $a->get_score();
}

1;
