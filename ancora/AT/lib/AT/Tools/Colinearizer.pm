############################################################################
# Colinearizer.pm - perl module for creating a colinear alignment of a 
# complosite alignment
#
# Copyright
############################################################################
 
package AT::Tools::Colinearizer;


=head1 NAME

Colinearizer

=head1 SYNOPSIS

    use AT::Tools::Colinearizer;
    my $colinearizer = AT::Tools::Colinearizer-> new(-CA => Composite_Alignment object);
    $composite_alignment -> get_noncolinear();
    $composite_alignment -> get_colinear();
    $composite_alignment -> colinearize();

=head1 DESCRIPTION

B<Colinearizer()> is a perl module for turning a Composite_Alignment object into a new Composite_Alignment object containing colinear subalignments only.

=head1 METHODS DESCRIPTION

=cut

############################################################################

use strict;
use vars '@ISA';
use Bio::SeqIO;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Getopt::Long;
#use AT::Root;
use AT::Alignment::Composite;
#use AT::Alignment::Subalignment;

@ISA = qw(AT::Root);

=head2 B<new()>

I<Title:> new()

I<Usage:> $colinearizer = AT::Tools::Colinearizer->new(-CA => Composite_Alignment object);

I<Input:> A B<AT::Alignment::Composite> object.

I<Output:> Returns a new B<AT::Tools::Colinearizer> object. 

I<Description:> 

=cut

sub new
{
    my ($class, %args) = @_;
    my ($CA);
    if (defined $args{-CA}){
      $CA=$args{-CA};
      delete $args{-CA};
    }
    my $self = bless {
		  _noncolinear => $CA,
		 max_overlap => 1000
		 }, ref $class || $class;
    $self->max_overlap($args{-max_overlap}) if (defined $args{-max_overlap});
    return $self;
}


=head2 B<get_noncolinear()>

I<Title:> get_noncolinear()

I<Usage:> $CA = $Colinearizer -> get_colinear();

I<Input:> Nothing

I<Output:> Returns a B<Composite_Alignment> object. 

I<Description:> 

=cut

sub get_noncolinear {$_[0] -> {_noncolinear} }
  

=head2 B<get_colinear()>

I<Title:> get_colinear()

I<Usage:> $colinear = $Colinearizer -> get_noncolinear();

I<Input:> Nothing

I<Output:> Returns a B<Composite_Alignment> object. 

I<Description:> 

=cut

sub get_colinear {
  my ($self,%args) = @_;
  unless($_[0] -> {_colinear}){
    $self->colinearize();
  }
$_[0] -> {_colinear}
}



=head2 B<colinearize()>

I<Title:> colinearize()

I<Usage:> $composite_alignment -> colinearize();

I<Input:> No inputs

I<Output:> No output

I<Description:> Finds the colinear path of subalignments giving the highest score of the composite alignment.

=cut


sub colinearize {
  my($self, %args) = @_;
  my $noncolin=$self->get_noncolinear();
 
  my $subalignments = $noncolin->get_subalignments();
 
  #start_end_score will keep the positions and score of each subalignment 
  #in this format:
  #start_seq1 end_seq1  start_seq2 end_seq2 score
  my @matches;
  foreach my $subalignment(@$subalignments){
    push @matches, { start1 => $subalignment->get_pos_seq1()->[0],
		     end1 => $subalignment->get_pos_seq1()->[1],
		     start2 => $subalignment->get_pos_seq2()->[0],
		     end2 => $subalignment->get_pos_seq2()->[1],
		     score => $subalignment->get_score(),
		     subalignment => $subalignment };
  }

  # sort the matches by one of the position fields
  # (which one is unimportant)
  @matches = sort {$a->{'start1'} <=> $b->{'start1'}} @matches;
  my $colin = AT::Alignment::Composite -> new(-seq1 => $noncolin->get_seq1(),
					      -seq2 => $noncolin->get_seq2(),
					      -strand => $noncolin->get_strand);

  # find optimal combination of matches
  my ($results ) = $self->_find_path ($colin, @matches);
  
  #sets the colinear argument:
  $self -> {_colinear} = $colin;

}
  

sub _find_path {
    my ($self, $colin, @matches) = @_;
    my $n = scalar @matches;
    my (@A, @B);
    # Some explanations:
    # $n is the number of matches that we consider
    # $A[$i] is the max score sum possible when using match $i and any
    # matches 1,2,..,$i-1 that are compatible with $i and eachother
    # $A[$i] is solved according to the recursion:
    #    $A[$i] = score for match $i + max($A[$k])
    # where the max is taken over all $k < $i that are compatible with $i
    # Note that if $k and $i are compatible, then $i is also compatible
    # with all matches used in the sum $A[$k] (otherwise this DP strategy
    # would not work).
    # $B[$i] is the $k used to solve $A[$i] (for backtrace).

    # Compute @A and @B
    for my $i (0..$n-1) {
      $B[$i] = -1;
      my $max_subsolution_score = 0;
      for my $j (0..$i-1) {
	if($matches[$j]->{end1} < $matches[$i]->{start1} + $self->{max_overlap}and
	   $matches[$j]->{end2} < $matches[$i]->{start2} + $self->{max_overlap}and
	   $A[$j] > $max_subsolution_score) {
	  $max_subsolution_score = $A[$j];
	  $B[$i] = $j;
	}
      }
      $A[$i] = $matches[$i]->{score} + $max_subsolution_score;
    }
    # Find highest score sum (max($A[$j]) over all $j) and backtrace from the
    # last match $j in it to generate results
    my (@results, @colinear_subalignments);
    my $j = 0;
    for my $i (1..$n-1) {
      $j = $i if($A[$i] > $A[$j]);
    }
    while ($j != -1) {
      unshift  @results, [ 
			  $matches[$j]->{'start1'},
			  $matches[$j]->{'end1'},
			  $matches[$j]->{'start2'},
			  $matches[$j]->{'end2'},
			  $matches[$j]->{'score'},
			  $matches[$j]->{'subalignment'},
			     $A[$j], $j, $B[$j]
			 ];
      push @colinear_subalignments, $matches[$j]->{'subalignment'};
      $j = $B[$j];
    }    
    foreach my $subaln (reverse @colinear_subalignments) {
	$colin->add_subalignment(-sub => $subaln);
    }

    return (\@results );
}
  

sub print_results  {
    
    foreach my $rowref (@_) {
      print join("\t", @$rowref)."\n";
    }
}



1;
