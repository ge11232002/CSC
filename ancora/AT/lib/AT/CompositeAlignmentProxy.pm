# AT::CompositeAlignmentProxy module
#
# Copyright Boris Lenhard
# You may distrubute this code under the same terms as perl itself
#
# POD documentation:

=head1 NAME

AT::CompositeAlignmentProxy

=head2 SYNOPSIS

=head2 DESCRIPTION

Adds some functionality to AT::Alignment:Composite.
These methods will probably be put in AT::Alignment::Composite.

=head2 TO DO

=head1 APPENDIX

=cut

# The code begins HERE

package AT::CompositeAlignmentProxy;

use strict;
use vars '@ISA';
use vars '$AUTOLOAD';
use Carp;


sub new  {
    my ($caller, %args) = @_;
    my $self = bless {
	aln => ($args{aln} or croak 'No aln arg'),
	seqs => ($args{seqs} or croak 'No seqs arg'),
    }, ref $caller || $caller;
    return $self;
}


=head2 column_from_residue_number

 Title     : column_from_residue_number
 Usage     : my $column = $aln->column_from_residue_number(1, 120);
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
	map {$_->get_alignment_obj} @{$self->aln->get_subalignments};
    
    # get seq and check that it is long enough to contain $resnr
    my $seq = $self->{seqs}->[$seq_pos-1];
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
    return AT::Alignment::Composite::Column->new( subaln => $i+1,
						  column => $column);
}


=head2 location_from_column

 Title     : location_from_columnq
 Usage     : my $location = $aln->location_from_column(1, $column);
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

    unless ($column and $column->isa('AT::Alignment::Composite::Column'))
    { croak 'Second arg missing or wrong object type'; }
    
    my $subalns = $self->aln->get_subalignments;

    if($column->is_within_subalignment) {
	unless (0 < $column->subaln and $column->subaln <= @$subalns)
	{ croak 'subaln out of range'; }
	$loc = $subalns->[$column->subaln-1]
	    ->get_alignment_obj->get_seq_by_pos($seqpos)
	    ->location_from_column($column->column);
    }
    else {
	unless  (0 < $column->subaln and $column->subaln <= @$subalns+1) 
	{ croak 'subaln out of range'; }
	my ($begin, $end);
	if($column->subaln > 1) {
	    $begin = $subalns->[$column->subaln-2]
		->get_alignment_obj->get_seq_by_pos($seqpos)->end;
	}
	if($column->subaln <= @$subalns) {
	    $end = $subalns->[$column->subaln-1]
		->get_alignment_obj->get_seq_by_pos($seqpos)->start;
	}
	if(defined($begin) and defined($end)) {    # added just this if clause
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


# This just forwards everything it doesn't know of to the actual alignment
# object
sub AUTOLOAD  {
    my ($self, @params) = @_;
    my ($arg) = $AUTOLOAD =~ /:([^:]+)$/;
    if (exists $self->{$arg})  {
	$self->{$arg} = $params[0] if (@params);
	return $self->{$arg};
    }
    else {
	return $self->aln->$arg(@params);
    }
}


sub DESTROY { }

1;


# Note: new package; should be but in separate file

package AT::Alignment::Composite::Column;


use strict;
use vars '@ISA';
use AT::Root;
use Carp;

@ISA = qw/AT::Root/;

# Conventions:
# Subalignment and column numbering start at 1.
# If column is 0, the position is between the indicated subaln and the
# previous.
# (subaln, column) = (1,0) means that the position is before the first
# subalignment.
# To allow indication of a position after the last subalignment, subalignment
# nubering can exceed the number of subalignments by 1.

sub new  {
    my ($caller, %args) = @_;
    croak 'No subaln arg' unless defined ($args{subaln});
    croak 'No column arg' unless defined ($args{column});
    my $self = bless {
	subaln => $args{subaln},
	column => $args{column},
    }, ref $caller || $caller;
    return $self;
}

# return true if column is within a subalignment
sub is_within_subalignment { $_[0]->{'column'} }
