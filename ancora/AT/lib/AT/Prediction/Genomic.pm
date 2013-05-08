# AT::Prediction::Genomic module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::Prediction::Genomic - abstract base class for genomic predictions

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::Prediction::Genomic;

use strict;
use vars '@ISA';
use AT::Root;
use AT::Tools::SeqHandler;
use Carp;

@ISA = qw(AT::Root);


sub new { croak "Method new should be defined in subclass"; }


# Subclasses must implement these methods

#sub seq { $_->throw_not_implemented(); }
#sub mapping_list { $_->throw_not_implemented(); }

# There's still some work to do on the methods below:
# - Must make sure that genomic_seq returns a locatable seq
#   i.e. if seq doesn't, we should cast seq into a locatableseq object
#   with abs_start and abs_end.
# - abs_start, abs_end, abs_strand can be made to work quicker by
#   having them calculate their values from relative start, end, strand
#   and the values of the attached seq
# - likewise, chr should use the attached seq instead of 'genomic_seq'
# - for this, we need each object to have rel_start, etc... these can
#   default to 1, the length of the parent seq and strand 1.
#   but should correspond to start, end, strand of a seqfeat.


=head2 chr, abs_start, abs_end, abs_strand

 Title     : chr, abs_start, abs_end, abs_strand
 Usage     : print  "Prediction on ", $pred->chr,
		    " strand ", $pred->abs_strand,
		    " from ", $pred->abs_start,
		    " to ", $pred->abs_end,
		    "\n";
 Function  : Get absolute location of prediction
 Returns   : Scalars.
	     chr returns a string on the form 'chr1'
	     abs_strand returns -1 or 1.
 Args      : -

=cut



sub assembly
{
    my ($self) = @_;
    unless ($self->{'assembly'}) {
	($self->{'assembly'}) = $self->genomic_seq->id =~ /(.+)_chr/;
    }
    return $self->{'assembly'};
}


sub chr
{
    my ($self) = @_;
    unless ($self->{'chr'}) {
	($self->{'chr'}) = $self->genomic_seq->id =~ /(chr.+)/;
    }
    return $self->{'chr'};
}

sub abs_start  { $_[0]->{'abs_start'} || $_[0]->seq->start; }
sub abs_end    { $_[0]->{'abs_end'} || $_[0]->seq->end; }
sub abs_strand { $_[0]->{'abs_strand'} || $_[0]->seq->strand; }


=head2 loc_str

 Title     : loc_str
 Usage     : print  "Prediction at ", $pred->loc_str, "\n";
 Function  : Return the the chromosomal location of the prediction
	     as a string.
 Returns   : A string on the form 'chr14:1112000-1114000:+'
 Args      : -

=cut


sub loc_str
{
    my ($self) = @_;

    my $str = $self->chr.":".$self->abs_start."-".$self->abs_end.":".
	($self->abs_strand == 1 ? '+' : '-');
    return $str;
}


sub id_loc_str
{
    my ($self) = @_;
    my $str = $self->genomic_seq->id.":".$self->abs_start."-".$self->abs_end.":".
	($self->abs_strand == 1 ? '+' : '-');
    return $str;
}


=head2 length

 Title     : length
 Usage     : print " Prediction spans ", $pred->length, " bp\n";
 Function  : Length of prediction in bp.
 Returns   : Scalar
 Args      : -

=cut

sub length { $_[0]->seq->length; } 
  # length is probably defined in a similar way in Bio::SeqFeature; it should
  # not matter that we define it here as well.

=head2 genomic_seq

 Title     : genomic_seq
 Usage     : my $seq = $pred->genomic_seq();
 Function  : Gets the genomic sequence spanned by the prediction
	     At the moment, I'm not sure why this method is not
	     simply called seq.
 Returns   : Bio::LocatableSeq
 Args      : -

=cut

sub genomic_seq { $_[0]->seq; }

=head2 oriented_genomic_seq

 Title     : oriented_genomic_seq
 Usage     : my $seq = $pred->oriented_genomic_seq();
 Function  : Like genomic_seq, except that the
	     sequence is revcom'd if its
	     strand is negative.
 Returns   : Bio::LocatableSeq
 Args      : -

=cut

sub oriented_genomic_seq {
    my ($self) = @_;
    my $seq = $self->genomic_seq;
    if($seq->strand == -1) {
	my $revseq = $seq->revcom;
	$revseq->id($seq->id);
	#$revseq->id($seq->id.'(revcom)');
	$revseq->start($seq->start);
	$revseq->end($seq->end);
	$revseq->strand($seq->strand);
	$seq = $revseq;
    }
    return $seq;
}

=head2 rel2abs

 Title     : rel2abs
 Usage     : my $abs_coord =
		$pred->rel2abs($rel_coord);
	     my ($abs_start, $abs_end) =
		$pred->rel2abs($rel_start, $rel_end);
 Function  : Conversion from relative to absolute coords.
	     Relative coords start at 1 and increase in the
	     direction of the prediction (e.g. direction of
	     transcription for a gene prediction).
	     Absolute coords represent chromosomal location
	     and always refer to the chromosomal plus strand.
	     As exemplified, the method can be called with
	     one or two coordinates. If two coords are given,
	     they are taken to represent a range. The first
	     coord must be less than or equal to the second.
	     For the returned range, the first coord
	     will always be less than or equal to the second.
 Returns   : One or two scalars
 Args      : One or two scalars

=cut


sub rel2abs
{
    my ($self, $rel_start, $rel_end) = @_;
    if($rel_start < 1 or $rel_start > $self->length or
       (defined($rel_end) and $rel_end > $self->length)) {
	croak('Invalid relative coord(s) '.$rel_start.
		(defined($rel_end) ? ('-'.$rel_end) : '').
		' (valid range = 1-'.$self->length.')');
    }
    return rel2abs_nobounds(@_);
}


sub rel2abs_nobounds
{
    my ($self, $rel_start, $rel_end) = @_;
    if(defined($rel_end) and $rel_end < $rel_start) {
	croak('Invalid relative coord(s) '.$rel_start.
		(defined($rel_end) ? ('-'.$rel_end) : '').
		' (valid range = 1-'.$self->length.')');
    }
    if($self->strand == 1) {
	my $abs_start = $self->abs_start + $rel_start - 1;
	return $abs_start unless (defined $rel_end);
	my $abs_end = $self->abs_start + $rel_end - 1;
	return ($abs_start, $abs_end);
    }
    else {
	my $abs_end = $self->abs_end - $rel_start + 1;
	return $abs_end unless (defined $rel_end);
	my $abs_start = $self->abs_end - $rel_end + 1;
	return ($abs_start, $abs_end);
    }
}


=head2 abs2rel

 Title     : abs2rel
 Usage     : my $rel_coord =
		$pred->abs2rel($abs_coord);
	     my ($rel_start, $rel_end) =
		$pred->abs2rel($abs_start, $abs_end);
 Function  : Conversion from absolute to relative coords.
	     See the the previous method for details.
 Returns   : One or two scalars
 Args      : One or two scalars

=cut


sub abs2rel
{
    my ($self, $abs_start, $abs_end) = @_;
    if($abs_start < $self->start or $abs_start > $self->end or
       (defined($abs_end) and $abs_end > $self->end)) {
	croak('Invalid absolute coord(s) '.$abs_start.
		(defined($abs_end) ? ('-'.$abs_end) : '').
		' (valid range = '.$self->start.'-'.$self->end.')');
    }
    return abs2rel_nobounds(@_);
}


sub abs2rel_nobounds
{
    my ($self, $abs_start, $abs_end) = @_;  
    if(defined($abs_end) and $abs_end < $abs_start) {
	croak('Invalid absolute coord(s) '.$abs_start.
		(defined($abs_end) ? ('-'.$abs_end) : '').
		' (valid range = '.$self->start.'-'.$self->end.')');
    }
    if($self->strand == 1) {
	my $rel_start = $abs_start - $self->abs_start + 1;
	return $rel_start unless (defined $abs_end);
	my $rel_end = $abs_end - $self->abs_start + 1;
	return ($rel_start, $rel_end);
    }
    else {
	my $rel_end = $self->abs_end - $abs_start + 1;
	return $rel_end unless (defined $abs_end);
	my $rel_start = $self->abs_end - $abs_end + 1;
	return ($rel_start, $rel_end);
    }
}


=head2 overlaps

 Title     : overlaps
 Usage     : if($pred1->overlaps($pred2)) { ... }
 Function  : Check for overlap with another prediction.
 Returns   : 0 if there is no overlap
	     1 if predictions overlap on the same strand
	     -1 if predictions overlap on opposite strands
 Args      : An AT::Prediction::Genomic-compliant object

=cut

sub overlaps
{
    my($self, $another) = @_;

    if($self->genomic_seq->id eq $another->genomic_seq->id and
       $self->abs_start <= $another->abs_end and
       $self->abs_end >= $another->abs_start) {
	return 1 if ($self->abs_strand == $another->abs_strand);
	return -1;
    }
    return 0;
}




