############################################################################
# Subalignment.pm - perl module for handeling a subalignment 
#
# Copyright
############################################################################
 
package AT::Alignment::Subalignment;


=head1 NAME

B<Subalignment>

=head1 SYNOPSIS

    use Subalignment;
    my $subalignment = AT::Alignment::Subalignment->new(-alignment_obj => $alignment_obj,
                                                        -score         => $score,
                                                        -strand        => $strand);

    my $conserved_regions   = $subalignment -> get_conserved_regions();
    my $conserved_alignment = $subalignment -> get_conserved_alignments();

=head1 DESCRIPTION

B<Subalignment> is a module for handeling a subalignment. A Subalignment in this context contains one B<Bio::SimpleAlign> object, and can be part of a longer composite alignment. Using this module it is possible to extract highly conserved regions in the alignment and B<Bio::SimplaAlign> objects corresponding to these conserved regions.

=head1 METHODS DESCRIPTION

=cut

############################################################################

use strict;
use Bio::SeqIO;
use Bio::Seq;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Getopt::Long;
use Carp;

=head2 B<new()>

I<Title:> new()

I<Usage:> $subalignment = AT::Alignment::Subalignment->new(-ao      => $alignment_obj,
                                                           -lengths => $lengths
                                                           -score   => $score,
                                                           -strand  => $strand);

I<Input:> 

$alignment_obj = B<Bio::SimpleAlign> object holding the actual alignment. The start and end positions in the original sequences of the sequences in the alignment object must be specified.

$lengths is a list containing the full lengths of the sequence1 and 2 respectively.

$score integer representing the score value of the subalignment.

$strand  string representing the strand of sequence 2 in the subalignment, can take the values "+" or "-".

I<Output:> Returns a new B<Subalignment> object.

=cut

sub new
{
  my ($class, %args) = @_;
  my ($alignment_obj, $score, $strand, $pos_seq1, $pos_seq2, $lengths, $tws, $cws);
  $tws=100;
  $cws=100;
  if (defined $args{-ps1}){
    #$pos_seq1=$args{-ps1};
    delete $args{-ps1};
  }
  if (defined $args{-ps2}){
    #$pos_seq2=$args{-ps2};
    delete $args{-ps2};
  }
  if (defined $args{-ao}){
    $alignment_obj=$args{-ao};
    my @seqs=$alignment_obj->each_seq();
 
    if ($seqs[0]->start && $seqs[0]->end()&& $seqs[1]->start && $seqs[1]->end()){
      $pos_seq1->[0]=$seqs[0]->start;
      $pos_seq1->[1]=$seqs[0]->end();
      $pos_seq2->[0]=$seqs[1]->start;
      $pos_seq2->[1]=$seqs[1]->end();
    }
    else{
      croak"Subalignment: The alignment objects must consist of LocatableSeqs with defined start and end positions! \n";
    }
    delete $args{-ao};
  }    
  if (defined $args{-score}){
    $score=$args{-score};
    delete $args{-score};
  }
   if (defined $args{-strand}){
    $strand=$args{-strand};
    delete $args{-strand};
  }
  if (defined $args{-lengths}){
    $lengths=$args{-lengths};
    delete $args{-lengths};
  }

  if (defined $args{-tws}){
    $tws=$args{-tws};
    delete $args{-tws};
  
  } 
if (defined $args{-cws}){
    $cws=$args{-cws};
    delete $args{-cws};
  
  }

  return bless {
		    _pos_seq1        => $pos_seq1,
		    _pos_seq2        => $pos_seq2,
		    _alignment_obj   => $alignment_obj,
		    _score           => $score,
		    _strand          => $strand,
		    _lengths         => $lengths,
		    _tws             => $tws,
		    _cws             => $cws,
		   }, ref $class || $class;
    }

=head2 B<get_pos_seq1()>

I<Title:> get_pos_seq1()

I<Usage:> $pos_seq1 = $subalignment -> get_pos1();

I<Input:> Nothing

I<Output:> Returns a list containing the start position and the end position in sequence 1 of the subalignment. 

=cut

sub get_pos_seq1          {$_[0] -> {_pos_seq1}}

=head2 B<get_pos_seq2()>

I<Title:> get_pos_seq2()

I<Usage:> $pos_seq1 = $subalignment -> get_pos2();

I<Input:> Nothing

I<Output:> Returns a list containing the start position and the end position in sequence 2 of the subalignment. These positions always refere to the B<LocatableSeq> object that was used to build the alignment object of the subalignment. This means that if the second sequence was reverse complemented by blast, the positions are calculated on the reverse complemented sequence. (Start and end positions should be recalculated using the length of sequence 2 for the corresponding positions in the original sequence.)

=cut

sub get_pos_seq2          {$_[0] -> {_pos_seq2}}


=head2 B<get_alignment_obj()>

I<Title:> get_alignment_obj()

I<Usage:> $alignment_obj = $subalignment -> get_alignment_obj();

I<Input:> Nothing

I<Output:> Returns the B<Bio::SimpleAlign> object corresponding to the SimpleAlign object.

=cut

sub get_alignment_obj     {$_[0] -> {_alignment_obj}}

=head2 B<get_score()>

I<Title:> get_score()

I<Usage:> $score = $subalignment -> get_score();

I<Input:> Nothing

I<Output:> Returns an integer representing the score value of the subalignment.

=cut

sub get_score             {$_[0] -> {_score}}

=head2 B<get_similarity()>

I<Title:> get_similarity()

I<Usage:> $similarity = $subalignment -> get_similarity();

I<Input:> Nothing

I<Output:> Returns a list containing the similarities of all possible "windows" longer than $window_size in the subalignment.

I<Description:> Calculates the similarities between the two sequences in the alignment for every possible "window" of size $window_size in the actual alignment. The similarities for all "windows" are pushed into a list which is used to calculate the dynamic threshold for similarity of the composite alignment in wich the subalignment object is included.

=cut

sub get_similarity{
  my ($self, %args)=@_;
  my $ws=$self->{_tws};
  unless($self->{_similarity}){
    
    # Separates the matchline of the alignment into an array,
    # where each element contains either match(*) or missmatch/gap(blank)
    my $matchline = $self->get_alignment_obj()->match_line();
    my @matchline=split //, $matchline;
    my @similarity;
    
    #If matchline is shorter than the definded window_size,
    #theshold gets the default value.
    if(scalar(@matchline)>$ws){
      #slides through the matchline...
      for my $i(0 .. (scalar(@matchline)-$ws)){
	$similarity[$i]=0;
	for my $j($i .. ($i+$ws-1)){
	  if($matchline[$j] eq '*'){
	    $similarity[$i]=$similarity[$i]+1;
	}
	}
	$similarity[$i]=$similarity[$i]/$ws;
      }
    }
    $self->{_similarity}=\@similarity;
  }
  return $self -> {_similarity};
}



=head2 B<get_threshold()>

I<Title:> get_threshold()

I<Usage:> $threshold = $subalignment -> get_threshold();

I<Input:> Nothing

I<Output:> Returns a real number representing the threshols for  similarity of the subalignment.

=cut

sub get_threshold         {$_[0] -> {_threshold}}



=head2 B<get_strand()>

I<Title:> get_strand()

I<Usage:> $conservedregs = $subalignment -> get_strand();

I<Input:> Nothing

I<Output:> Returns a string representing the trand of sequence 2, can be either "+" or "-".

=cut

sub get_strand            {$_[0] -> {_strand}}

=head2 B<get_lenghts()>

I<Title:> get_lengths()

I<Usage:> $lengths = $subalignment -> get_lengths();

I<Input:> Nothing

I<Output:> Returns a list containg the length of sequence 1 and the length of sequence 2.

=cut

sub get_lengths           {$_[0] -> {_lengths}}


=head2 B<get_conserved_regions()>

I<Title:> get_conserved_regions()

I<Usage:> $conserved = $subalignment -> get_conserved_regions();

I<Input:> Nothing

I<Output:> 

Returns a reference to a  two dimentional list containing all the positions in the two original sequences of the conserved regions in the subalignment, as well as the similarities between the two sequences in the conserved regions. Every row in the output array is a reference to an array corresponding to one conserved region, having the format:

(start pos seq1, end pos seq1, start pos seq2, end pos seq2, percent identity).

Note, the positions in sequence two here refere to positions in the B<original sequence> that was used as input to BlastZ, indepented of the strand of sequence 2 in the blastz alignment.

=cut


sub get_conserved_regions     {
  my ($self, %args)=@_;
  unless($self->{_conserved_regions}){
    $self->conserved_regions();
  }
  return $self -> {_conserved_regions};
}

=head2 B<get_conserved_alignments()>

I<Title:> get_conserved_alignments()

I<Usage:> $conserved_alignments = $subalignment -> get_conserved_alignments();

I<Input:> Nothing

I<Output:> a reference to a list of b<Bio::SimpleAlign> objects.

I<Description:> Returns a reference to an array containing alignment objects corresponding to the conserved regions.

=cut

sub get_conserved_alignments    { 
  my ($self, %args)=@_;
  unless($self->{_conserved_alignments}){
    $self->conserved_regions();
  }
  return $self -> {_conserved_alignments};
}



=head2 B<set_threshold()>

I<Title:> set_threshold()

I<Usage:> $subalignment -> set_threshold (-th = $threshold);

I<Input:> Real number representing the threshold for similarity when searching for highly similar regions in the subalignment.

I<Output:> Nothing is returned

I<Description:> Sets the threshold attribute of the B<Subalignment> object.

=cut

sub set_threshold{
  my ($self, %args) = @_;
  my $threshold;
  if (defined $args{-th}){
    $threshold=$args{-th};
    delete $args{-th};
  }
  $self->{_threshold}=$threshold;
}

=head2 B<set_cws()>

I<Title:> set_cws()

I<Usage:> $subalignment -> set_cws( -ws = $window_size);

I<Input:> Real number representing the lenght of the window used when searching the alignment for conserved regions.

I<Output:> Nothing is returned

=cut

sub set_cws{
  my ($self, %args) = @_;
  my $cws;
  if (defined $args{-cws}){
    $cws=$args{-cws};
    delete $args{-cws};
  }
  $self->{_cws}=$cws;
}

=head2 B<set_tws()>

I<Title:> set_tws()

I<Usage:> $subalignment -> set_tws( -ws = $window_size);

I<Input:> Real number representing the lenght of the window used when calculating similarities.

I<Output:> Nothing is returned

=cut

sub set_tws{
  my ($self, %args) = @_;
  my $tws;
  if (defined $args{-tws}){
    $tws=$args{-tws};
    delete $args{-tws};
  }
  $self->{_tws}=$tws;
}

=head2 B<conserved_regions()>

I<Title:> conserved_regions()

I<Usage:> $subalignment -> conserved_regions(-ws => $window_size);

I<Input:> No input

I<Output:> No output

I<Description:> The subalignment is searched for regions that are longer than $window_size and more similar than $threshold. The positions in the original input sequences for every conserved region in the subalignment is calculated, and pushed into a two-dimentional array. Creates alignment objects for every conserved region in the subalignment. This method sets the values of the attributes _conserved_regions and _conserved_alignments of the B<Subalignment> object.

=cut

sub conserved_regions{
  my ($self, %args)=@_;
  my $ws = $self ->{_cws};

  #For storing the positions of conserved regions:
  my (@cons_regs);
  
  #For storing the conserved alignments:
   my @cons_alignments;

  #For keeping track of highly conserved regions:
  my ($percent_id, $similar, $start, $last_equal, $end);
  my $previous_end=0;
  my $previous_start=0;

  # Separates the matchline into an array where each element contains either match(*) 
  # or missmatch/gap(blank)
  my $alignment=$self->get_alignment_obj();
  my @seqs=$alignment->each_seq();
  my $matchline = $alignment->match_line(); 
  my @matchline=split //, $matchline;

  #looping through the whole match line
  for my $i(0 .. (scalar(@matchline)-($ws+1))){
    if($matchline[$i] eq '*'){

      #When an initial match character is found, a window of length $ws is first examined.
      my $first_window=substr($matchline,$i,$ws);
      my $similar=($first_window =~tr/\*//);
      my $start=$i;
      my $j=$i+$ws;
      while (substr($matchline, $j-1, 1) ne '*'){
	$j--;
      }
      $last_equal=$j;
      $percent_id=$similar/$ws;
      $i=$i+$ws+1;

      #Extension of the window until % identity falls below the threshold.
      while($i<(scalar(@matchline)-1) && $percent_id > $self->get_threshold()){
	if($matchline[$i] eq '*'){
	  $similar++;
	  $last_equal=$i;
	}
	$percent_id=$similar/($i-$start);
	$i++;
      }

      #In the end, only the percent identity for the region
      #between the first and last matches is of interest:
      if($last_equal-$start >= $ws){
	$percent_id=$similar/($last_equal-$start);
      }

      #Is the region conserved enough, does it end after the previous conserved region?
      #Then the region is reported as conserved
      $similar++;
      if ($percent_id > $self->get_threshold() && $last_equal>$previous_end){

	#Does this region overlap whith the previous one? In that case,
	#if the combined percent identity is above threshold, 
	#the two regions should be merged
	if($previous_end>$start){
	  my $new_region     = substr ($matchline, $previous_start, ($last_equal-$previous_start));
	  my $new_sim        = ($new_region =~tr/\*//);
	  my $new_percent_id = $new_sim / ($last_equal-$previous_start);
	  if ($new_percent_id > $self->get_threshold()){
	    $start=$previous_start;
	    $percent_id=$new_percent_id;
	    #removes the last conserved region so that it can be replaced by the new merged one
	    pop @cons_regs;
	    pop @cons_alignments;
	  }
	}

	#In sequence objects, the positions are numbered from "1" instread of "0":
	$start += 1;
	$end   = $last_equal +1;

	#Extracts the high scoring alignment and Calculates the 
	#reference positions (positions in the original sequence)
	#of the high scoring regions:
	my $conserved_alignment=$alignment->slice($start, $end);
	my @conserved_seqs=$conserved_alignment->each_seq();
	my @conserved_regions = ($conserved_seqs[0]->start,$conserved_seqs[0]->end, $conserved_seqs[1]->start,$conserved_seqs[1]->end, $percent_id) ;
	
	#Recalculation of positions if sequence 2 is reverse coplemented:
	if($self->get_strand() eq "-"){
	  my $oldstart = $conserved_regions[2];
	  my $length=$self->get_lengths()->[1];
	  $conserved_regions[2] = $length-$conserved_regions[3];
	  $conserved_regions[3] = $length-$oldstart;
	}
	push @cons_regs,\@conserved_regions;
	push @cons_alignments, $conserved_alignment;	
	$previous_end=$last_equal;
	$previous_start=$start-1;
      }
    }
  }
$self->{_conserved_regions}= \@cons_regs;
$self->{_conserved_alignments}=\@cons_alignments;
}






1;
