############################################################################
# Composite_Alignment.pm - perl module for handeling a set of subalignments 
#
# Copyright
############################################################################
 
package AT::Alignment::Composite;


=head1 NAME

Composite_Alignment

=head1 SYNOPSIS

    use AT::Alignment::Composite;
    my $composite_alignment = AT::Alignment::Composite-> new(-seq1   => $seq1,
					   	             -seq2   => $seq2,
							     _tws    => $tws,
							     _cws    => $cws,
							     _top    => $top,
							     _strand => $strand);

    $conserved_regions    = $composite_alignment -> get_conserved_regions();
    $conserved_alignments = $composite_alignment -> get_conserved_alignments();
    $local_alignment      = $composite_alignment -> get_local_alignment(-p => $p, -e => $e, -s => $s);

=head1 DESCRIPTION

B<AT::Alignment::Composite> is a perl module for handeling a set of subalignments.A B<AT::Alignment::Composite> object holds a list of subalignments. Subalignments can be added to the list, threshold for conservation can be calculated and the  subalignments can be searched for conserved regions.

=head1 METHODS DESCRIPTION

=cut

############################################################################


use vars '@ISA';
use Carp;
use strict;
use Bio::SeqIO;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::Location::Fuzzy;
use AT::Alignment::CompositeColumn;
use Getopt::Long;
#use AT::Alignment::Subalignment;

=head2 B<new()>

I<Title:> new()

I<Usage:> $composite_alignment = AT::Alignment::Composite->new(-seq1   => $seq1,
                                                               -seq2   => $seq2,
                                                               _tws    => $tws,
		                                               _cws    => $cws,
		                                               _top    => $top,
		                                               _strand => $strand);
);



I<Input:> 

$seq1 = B<Bio::Seq> object representing the first aligned sequence.

$seq2 = B<Bio::Seq> object representing the second aligned sequence If the second sequnce is reverse complemented in the alignment, $seq2 should represent the reverse complementation instread of the original inpustequence to BlastZ.

$strand is the strand of sequence two in the alignment, can take the values "-" and "+".

I<Optional:>

$tws it the window size used in calculation of threshold.

$cws is the window size used when finding conserved regions.

$top is the stringency of the dynamic threshold.

I<Output:> Returns a new B<AT::Alignment::Composite> object.

=cut

sub new
  {
    my ($class, %args) = @_;
    my ($seq1,$seq2, $strand);
    my $tws=100;
    my $top=0.1;
    my $cws=100;
    if (defined $args{-seq1}){
      $seq1=$args{-seq1};
      delete $args{-seq1};
    }
    if (defined $args{-seq2}){
      $seq2=$args{-seq2};
    delete $args{-seq2};
    }
    if (defined $args{-tws}){
      $tws=$args{-tws};
      delete $args{-tws};
    }
    if (defined $args{-cws}){
      $cws=$args{-cws};
      delete $args{-cws};
    }
    if (defined $args{-top}){
      $top=$args{-top};
      delete $args{-top};
    }
    if (defined $args{-strand}){
      $strand=$args{-strand};
      delete $args{-strand};
    }
    return bless {
		  _seq1          => $seq1,
		  _seq2          => $seq2,
		  _score         => 0,
		  _tws           => $tws,
		  _cws           => $cws,
		  _top           => $top,
		  _strand        => $strand,
		 }, ref $class || $class;
  }

=head2 B<add_subalignment()>

I<Title:> add_subalignment()

I<Usage:> $composite_alignment -> add_subalignment(-sub => $sub);

I<Input:> 

$sub is a reference to a B<AT::Alignment::Subalignment> object.

I<Output:> No output

I<Description:> The input B<AT::Alignment::Subalignment> object is added to the list of alignments in the alignmentset.

=cut

sub add_subalignment{
  my ($self, %args )=@_;
  my $subalignment;
  if(defined $args{-sub}){
    $subalignment=$args{-sub};
    delete $args{-sub};
  }
  my $score = $self->get_score()+$subalignment->get_score();
  my @temp_sub;
  if($self->get_subalignments()){
    my $subref=$self->get_subalignments();
    @temp_sub=@$subref;
  }
  push @temp_sub, $subalignment;
  $self->{_subalignments}=\@temp_sub;
  #Sets the new score value for the composite alignment:
  $self->{_score}=$score;
  #Signify that the similarity threshold needs to be recalculated:
  $self->{_threshold} = undef;
  #Sets the new similarity threshold for the composite alignment:
  #$self->set_threshold();
}

=head2 B<set_strand()>

I<Title:> set_strand()

I<Usage:> $composite_alignment ->set_strand(-strand => $strand);

I<Input:> 

$strand is a string representing the strand of sequence two in the alignment. Can take the values "+" or "-".

I<Output:> No output

I<Description:> The strand of the second sequence in current composite alignment is set to "+" or "-".

=cut

sub set_strand{
  my ($self, %args)=@_;
  my $strand;
  if(defined $args{-strand}){
    $strand=$args{-strand};
    delete $args{-strand};
    $self->{_strand}=$strand;
  }
}

=head2 B<set_tws()>

I<Title:> set_tws()

I<Usage:> $composite_alignment ->set_tws(-tws => $tws);

I<Input:> 

$tws is the window size for calculating threshold.

I<Output:>

I<Description:> Sets the value of window size for calculating threshold. The default value is 100 bp. Automaticly sets the new value of $tws for all subalignments.

=cut

sub set_tws{
  my ($self, %args)=@_;
  my $tws;
  if(defined $args{-tws}){
    $tws=$args{-tws};
    delete $args{-tws};
    my $subalignments=$self->get_subalignments();
    foreach my $subalignment(@$subalignments){
     $subalignment->set_tws(-tws => $tws);
    }
    $self->{_tws}=$tws;
  }
  #$self->{_tws}=$tws;
}

=head2 B<set_cws()>

I<Title:> set_cws()

I<Usage:> $composite_alignment ->set_cws(-cws => $cws);

I<Input:> 

$cws is the minimum length (window size) of conserved regions.

I<Output:> No output

I<Description:> Sets the walue of the window size for identification of conserved regions. The default value is 100 bp. Automaticly sets the new value of $cws for all subalignments.


=cut

sub set_cws{
  my ($self, %args)=@_;
  my $cws;
  if(defined $args{-cws}){
    $cws=$args{-cws};
    delete $args{-cws};
    my $subalignments=$self->get_subalignments();
    foreach my $subalignment(@$subalignments){
      $subalignment->set_cws(-cws => $cws);
    }
  }
  $self->{_cws}=$cws;
}

=head2 B<set_top()>

I<Title:> set_top()

I<Usage:> $composite_alignment ->set_top(-top => $top);

I<Input:>

$top is the stringency for dynamic calculation of threshold.

I<Output:> No output

I<Description:> The threshold for similarity will be set to a value so that the best "$top"% scoring words of lenght $cws will be repported as conserved.

=cut

sub set_top{
  my ($self, %args)=@_;
  my $top;
  if(defined $args{-top}){
    $top=$args{-top};
    delete $args{-top};
  }
  $self->{_top}=$top;
}

=head2 B<set_threshold()>

I<Title:> set_threshold()

I<Usage:> $composite_alignment -> set_thershold(-th => $threshold);

          Also used internally without input parameters to calculate the threshold dynamically.

I<Input:> 

optional: integer representing a manually defined threshold.

I<Output:> No output

I<Description:> Sets the threshold used for identification of conserved region in the composite alignment.
When calling this function, the threshold is automatically set for all subalignments in the composite alignment.
If the method is called with the argument -th representing a manually chosen threshold, the threshold will be set to this value.
When calling B<set_threshold> without inputparameters, the threshold will be calculated dynamically according to the quality of
the alignment, and the values of _cws and _top. The threshold is then chosen so that a certain fraction (= top) of words of
length >= $window size in the composite alignment are more similar than the threshold. B<set_threshold()> uses the similarity
list for the composite alignment to find the apropriate threshold. B<set_threshold> is used internally when new subalignments
are added to the composite alignments.

=cut


sub set_threshold{
  my ($self, %args)=@_; 
  #Default value
  my $threshold=0.3;

  #All subalignments;
  my $subalignments = $self->get_subalignments();

  #If manually chosen threshold:
  if (defined $args{-th}){
    $threshold=$args{-th};
    delete $args{-th};
  }

  #Else, compute threshold dynamically:
  else{
    my $ws = $self->{_tws};
    my $top = $self->{_top};
    my @similarity;
    #looping through all subalignments
    foreach my $subalignment(@$subalignments){
      my $subsim = $subalignment->get_similarity();
      push @similarity,@$subsim;
    }
    #sorting the array of "percent identity":
    my @sorted_similarity=sort{$a <=> $b} @similarity;
    my $no_of_words=scalar(@sorted_similarity);
    my $threshold_position=sprintf("%.0f", ($no_of_words-$top*$no_of_words));
    $threshold=$sorted_similarity[$threshold_position];
  }
  
  for my $sub(@$subalignments){
    $sub->set_threshold(-th => $threshold);
  }
  $self->{_threshold}=$threshold;
}





=head2 B<get_seq1()>

I<Title:> get_seq1()

I<Usage:> $seq1 = $composite_alignment -> get_seq1();

I<Input:> Nothing

I<Output:> Returns a B<Bio::Seq> object. 

I<Description:> 

=cut

sub get_seq1          {$_[0] -> {_seq1} }

=head2 B<get_seq2()>

I<Title:> get_seq2()

I<Usage:> $seq2 = $composite_clignment -> get_seq2();

I<Input:> Nothing

I<Output:> Returns a B<Bio::Seq> object. 

I<Description:> 

=cut

sub get_seq2          {$_[0] -> {_seq2} }

=head2 B<get_strand()>

I<Title:> get_strand()

I<Usage:> $strand = $composite_alignment -> get_strand();

I<Input:> Nothing

I<Output:> Returns a string with value "+" or "-" depending on the strand of sequence two in the alignment. 

I<Description:> 

=cut

sub get_strand        {$_[0] -> {_strand}}

=head2 B<get_subalignments()>

I<Title:> get_subalignments()

I<Usage:> $subaligments = $composite_alignments -> get_subalignments();

I<Input:> Nothing

I<Output:> Returns a reference to a list of B<Subalignment> objects.

I<Description:> 

=cut

sub get_subalignments {$_[0] -> {_subalignments}}


=head2 B<get_threshold()>

I<Title:> get_threshold()

I<Usage:> $threshold = $composite_alignment -> get_threshold();

I<Input:> Nothing

I<Output:> Returns the threshold of similarity for the alignment

I<Description:> 

=cut

sub get_threshold     {
  my ($self, %args )=@_;
  unless($_[0] -> {_threshold}){
    $self->set_threshold();
  }
  $_[0] -> {_threshold}
}



=head2 B<get_score()>

I<Title:> get_score()

I<Usage:> $score = $composite_alignment->get_score();

I<Input:> Nothing

I<Output:> Returns the score value for the alignment. 

I<Description:> 

=cut

sub get_score         {$_[0] -> {_score}}


=head2 B<get_conserved_regions()>

I<Title:> get_conserved_regions()

I<Usage:> $conserved_regions = $composite_alignment -> get_conserved_regions();

I<Input:> Nothing

I<Output:> Returns a reference to a  two dimentional list containing all the positions in the two original sequences of the conserved regions in the subalignment, as well as the similarities between the two sequences in the conserved regions. Every row in the output array is a reference to an array corresponding to one conserved region, having the format:

(start pos seq1, end pos seq1, start pos seq2, end pos seq2, percent identity).

Note, the positions in sequence two here refere to positions in the B<original sequence> that was used as input to BlastZ, indepented of the strand of sequence 2 in the blastz alignment.

=cut

sub get_conserved_regions            {
  my ($self, %args )=@_;
  unless($_[0] -> {_threshold}){
    $self->set_threshold();
  }
  unless($_[0] -> {_conserved_regions}){
    $self->conserved_regions();
  }
  $_[0] -> {_conserved_regions}
}


=head2 B<get_conserved_alignments()>

I<Title:> get_conserved_alignments()

I<Usage:> $conserved_alignments = $composite_alignment -> get_conserved_alignments();

I<Input:> Nothing

I<Output:> Returns a reference to a list of all b<Bio::SimpleAlign> object corresponding to the conserved regions.

I<Description:> 

=cut

sub get_conserved_alignments            { 
  my ($self, %args )=@_;
  unless($_[0] -> {_threshold}){
    $self->set_threshold();
  }
  unless($_[0] -> {_conserved_regions}){
    $self->conserved_regons();
  }
  $_[0] -> {_conserved_alignments}
}





=head2 B<get_local_alignment()>

I<Title:> get_local_alignment()

I<Usage:> $local_alignment = $composite_alignment ->get_local_alignment(-p => $p, -e => $e, -s => $s);

I<Input:> 

$p = Integer representing the position of interest.

$e = Integer representing extension of alignment on both side of position.

$s = Integer representing reference species.

I<Output:> A list of Alignment objects from all subalignments covering the specified region.

I<Description:> Loops through all subalignments to see wheather the position of interest is covered by a subalignment. (The position should refere to the original input sequences.) If it is, the alignment For E bases on both sides of P is extracted. The position can be covered by several subalignments, in wich case several subalignments will be returned.

=cut


sub get_local_alignment{
  my ($self, %args) = @_;
  my ($p, $e, $s);
  if(defined $args{-p}){
    $p = $args {-p};
  }
  if(defined $args{-e}){
    $e = $args {-e};
  }
  if(defined $args{-s}){
    $s = $args{-s};
  }
  my @subAlignObj;

  my $subalignments = $self->get_subalignments();
  foreach my $subalignment(@$subalignments){
    my $positions_seq1 = $subalignment->get_pos_seq1();
    my $positions_seq2 = $subalignment->get_pos_seq2();
    my $subAlignObj=0;

    if ($s eq "h"){
      if ( ($positions_seq1->[0] < $p) && ($positions_seq1->[1] > $p)){
	if($positions_seq1->[0] < ($p-$e) && $positions_seq1->[1]>($p+$e)){
	  my $alignmentObj = $subalignment->get_alignment_obj();
	  my $name=$self->get_seq1()->id;
	  my $start = $alignmentObj -> column_from_residue_number( $name, $p-$e );
	  my $end = $alignmentObj -> column_from_residue_number( $name, $p+$e ); 
	  $subAlignObj = $alignmentObj->slice( $start, $end);
	}
	elsif($positions_seq1->[0] < ($p-$e) && $positions_seq1->[1] < ($p+$e)){
	  my $alignmentObj = $subalignment->get_alignment_obj();
	  my $name=$self->get_seq1()->id;
	  my $start = $alignmentObj -> column_from_residue_number( $name, $p-$e );
	  my $end = $alignmentObj -> column_from_residue_number( $name, $positions_seq1->[1] );
	  $subAlignObj = $alignmentObj->slice( $start, $end);
	}
	elsif($positions_seq1->[0] > ($p-$e) && $positions_seq1->[1] > ($p+$e)){
	  my $alignmentObj = $subalignment->get_alignment_obj();
	  my $name=$self->get_seq1()->id;
	  my $start = $alignmentObj -> column_from_residue_number( $name,$positions_seq1->[0]);
	  my $end = $alignmentObj -> column_from_residue_number( $name, $p+$e );
	  $subAlignObj = $alignmentObj->slice( $start, $end);
	}	
	elsif($positions_seq1->[0] > ($p-$e) && $positions_seq1->[1] < ($p+$e)){
	  my $alignmentObj = $subalignment->get_alignment_obj();
	  my $name=$self->get_seq1()->id;
	  my $start = $alignmentObj -> column_from_residue_number( $name,$positions_seq1->[0]);
	  my $end = $alignmentObj -> column_from_residue_number( $name, $positions_seq1->[0]);
	  $subAlignObj = $alignmentObj->slice( $start, $end);
	}
      }
    }
    elsif($s eq "m"){ 
      my ($seq2_start, $seq2_end);
      #Recalculation of position if sequence 2 is reverse coplemented:
      if($self->get_strand eq "-"){
	my $length=$self->get_seq2()->length;
	$p = $length - $p +1;
      }
      if ( ($positions_seq2->[0] < $p) && ($positions_seq2->[1] > $p)){
	if($positions_seq2->[0] < ($p-$e) && $positions_seq2->[1]>($p+$e)){
	  my $alignmentObj = $subalignment->get_alignment_obj();
	  my $name=$self->get_seq2()->id;
	  my $start = $alignmentObj -> column_from_residue_number( $name, $p-$e );
	  my $end = $alignmentObj -> column_from_residue_number( $name, $p+$e );
	#   print "Start: ".$start."\tend: ".$end."\n";
	  $subAlignObj = $alignmentObj->slice( $start, $end);
	}
	elsif($positions_seq2->[0] < ($p-$e) && $positions_seq2->[1] < ($p+$e)){
	  my $alignmentObj = $subalignment->get_alignment_obj();
	  my $name=$self->get_seq2()->id;
	  my $start = $alignmentObj -> column_from_residue_number( $name, $p-$e );
	  my $end = $alignmentObj -> column_from_residue_number( $name, $positions_seq2->[1] );
	  $subAlignObj = $alignmentObj->slice( $start, $end);
	}
	elsif($positions_seq2->[0] > ($p-$e) && $positions_seq2->[1] > ($p+$e)){
	  my $alignmentObj = $subalignment->get_alignment_obj();
	  my $name=$self->get_seq2()->id;
	  my $start = $alignmentObj -> column_from_residue_number( $name, $positions_seq2->[0]);
	  my $end = $alignmentObj -> column_from_residue_number( $name, $p+$e );
	  $subAlignObj = $alignmentObj->slice( $start, $end);
	}
	elsif($positions_seq2->[0] > ($p-$e) && $positions_seq2->[1] < ($p+$e)){
	  my $alignmentObj = $subalignment->get_alignment_obj();
	  my $name=$self->get_seq1()->id;
	  my $start = $alignmentObj -> column_from_residue_number( $name,$positions_seq1->[0]);
	  my $end = $alignmentObj -> column_from_residue_number( $name, $positions_seq1->[0]);
	  $subAlignObj = $alignmentObj->slice( $start, $end);
	}
      }
    }
    if($subAlignObj !=0){
      push @subAlignObj, $subAlignObj;
    }
  }
  return \@subAlignObj;
}



=head2 column_from_residue_number

 Title     : column_from_residue_number
 Usage     : my $column = $composite_alignment->column_from_residue_number(1, 120);
 Function  : Gives the position in the alignment of a given residue
             in a given sequence. (Meant to be similar to the
             Bio::SimpleAlign method.)
 Returns   : An AT::Alignment::Composite::Column object
 Args      : 1. position of the sequence in the alignment (1 or 2)
                A small and probably irrelevant note:
                Bio::SimpleAlign uses sequence id instead of
                position. We cannot do that since we work with
                relative coords. If the aligned sequences are
                different subsequences of one larger seq, it is
                likely that they will have the same id.
                Bio::SimpleAlign would the distinguish between the
                seqs using their absolute coords, which we cannot do.
             2. residue number (value btw 1 and length of seq)

=cut


sub column_from_residue_number
{
    my ($self, $seq_pos, $resnumber) = @_;

    croak("First argument sequence position missing or invalid")
        unless ($seq_pos and $seq_pos >= 1 and $seq_pos <= 2);
    croak("Second argument residue number missing") unless $resnumber;

    my @simplealns =
        map {$_->get_alignment_obj} @{$self->get_subalignments};
    my $seq;
    # get seq and check that it is long enough to contain $resnr
    if($seq_pos == 1){
      $seq = $self->get_seq1();
    }
    elsif($seq_pos == 2){
      $seq = $self->get_seq2();
    }
    if ($resnumber > $seq->length)
    { croak ("Second argument residue number [$resnumber] ".
        "larger than sequence length [".$seq->length."]"); }
    
    # find column of the resnumber in the seq
    my $column = 0;
    my $i = 0;
    for(;$i < @simplealns; $i++) {
        my $alnseq = $simplealns[$i]->get_seq_by_pos($seq_pos);
        if ($resnumber <= $alnseq->end) {
            if ($resnumber >= $alnseq->start) {
                $column = $alnseq->column_from_residue_number($resnumber);
            }
            last;
        }
    }
    return AT::Alignment::CompositeColumn->new( subaln => $i+1,
						column => $column);
}


=head2 location_from_column

 Title     : column_from_residue_number
 Usage     : my $location = $composite_alignment->location_from_column(1, $column);
 Function  : Gives the coordinates for a given sequence at a given
             column in the alignment. If the column is in a gap,
             a range is returned.
             (Meant to be similar to the Bio::LocatableSeq method.)
 Returns   : A Bio::LocationI-compliant object
 Args      : 1. position of the sequence in the alignment (1 or 2)
                See column_from_residue_number for why we use position
                instead of sequence id.
             2. an AT::Alignment::Composite::Column object

=cut


sub location_from_column
{
    my ($self, $seqpos, $column) = @_;
    my $loc;

    unless ($column and $column->isa('AT::Alignment::CompositeColumn'))
    { croak 'Second arg missing or wrong object type'; }
    
    my $subalns = $self->get_subalignments();

    if($column->is_within_subalignment) {
        unless (0 < $column->{subaln} and $column->{subaln} <= @$subalns)
        { croak 'subaln out of range'; }
        $loc = $subalns->[$column->{subaln}-1]
            ->get_alignment_obj->get_seq_by_pos($seqpos)
            ->location_from_column($column->{column});
    }
    else {
        unless  (0 < $column->{subaln} and $column->{subaln} <= @$subalns+1) 
        { croak 'subaln out of range'; }
        my ($begin, $end);
        if($column->{subaln} > 1) {
            $begin = $subalns->[$column->{subaln}-2]
                ->get_alignment_obj->get_seq_by_pos($seqpos)->end;
        }
        if($column->{subaln} <= @$subalns) {
            $end = $subalns->[$column->{subaln}-1]
                ->get_alignment_obj->get_seq_by_pos($seqpos)->start; 
        }
	if(defined($begin) and defined($end)) { 
	    $loc = Bio::Location::Simple->new
		(-start => $begin,
		 -end => $end,
		 -location_type => ($begin == $end) ? 'IN_BETWEEN' : 'EXACT');
	}
	else {
	    $loc = Bio::Location::Fuzzy->new(-start => ($begin or '<'.$end),
					-end => ($end or '>'.$begin),
					-location_type => 'BETWEEN');
	}	
    }
    return $loc;
}


=head2 B<conserved_regions()>

I<Title:> conserved_regions()

I<Usage:> Used internally

I<Input:> No input

I<Output:> No output

I<Description:> Every subalignment of the composite alignment is searched for regions that are longer than $window_size and more similar than $threshold. The corresponding positions in the original sequence of these conserved regions are then calculated. Alignment objects are created for every conserved region.

=cut

sub conserved_regions{
  my ($self, %args) =@_;
  my $ws=$self->{_cws};
  my @conserved_regions;
  my @conserved_alignments;
  my $subalignments=$self->get_subalignments();
  foreach my $subalignment(@$subalignments){
    my $temp_conserved_regions= $subalignment->get_conserved_regions();
    my $temp_alignments = $subalignment -> get_conserved_alignments();
    for my $i(0 .. scalar (@$temp_conserved_regions)-1){
      push @conserved_regions, $temp_conserved_regions->[$i];
      push @conserved_alignments, $temp_alignments->[$i];
    }
    $self->{_conserved_regions}=\@conserved_regions;
    $self->{_conserved_alignments}=\@conserved_alignments;
  }
}


sub DESTROY { }



1;
