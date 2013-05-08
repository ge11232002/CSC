# AT::CMP::CorrespCoords module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::CMP::CorrespCoords - corresponding coordinate ranges in different
genomic predictions

=head1 SYNOPSIS

 print  "the range ",
	join '-', $cc->rel_coords(pos => 1),
	" in gene 1 ",
	" corresponds to the range ",
	join '-', $cc->rel_coords(pos => 2),
	" in gene 2.\n";

 # See AT::CMP::ComparisonFactory::GenomicAln for a larger usage example.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::CMP::CorrespCoords;

use strict;
use vars '@ISA';
use AT::Root;
use Carp;

@ISA = qw(AT::Root);

=head2 new

 Title     : new
 Usage     : 
 Function  : Constructor
	     Do not call this method directly. Create object using
	     a Comparison.
 Returns   : AT::CMP::CorrespCoords
 Args      : 

=cut

sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
	_coords => ($args{'coords'} or croak 'No coords arg')
    }, ref $caller || $caller;

    return $self;
}


=head2 nr_coords

 Title     : nr_coords
 Usage     : print "I have ", $cc->nr_coords, " coordinate pairs";
 Function  : Returns number of coordinate pairs (equal to number
	     of objects in the comparison). Some of the coordinate
	     pairs may be undefined.
 Returns   : Scalar
 Args      : -

=cut

sub nr_coords
{
    my ($self) = @_;
    return scalar @{$self->{'_coords'}};
}


=head2 abs_coords

 Title     : abs_coords
 Usage     : my ($start, $end) = $cc->abs_coords(pos => 1);
 Function  : Get absolute coords for the range in a given bject
 Returns   : Array of two integers
 Args      : Supply one of:
	     obj   Reference to object to get coords for
	     pos   Position in comparison of object to
		   get coords for (1 - nr of objects)

=cut

sub abs_coords
{
    my ($self, %args) = @_;
    return $self->_get_coords(%args, abs => 1);
}


=head2 rel_coords

 Title     : rel_coords
 Usage     : my ($start, $end) = $cc->rel_coords(pos => 1);
 Function  : Get relative coords for the range in a given object
 Returns   : Array of two integers
 Args      : Supply one of:
	     obj   Reference to object to get coords for
	     pos   Position in comparison of object to
		   get coords for (1 - nr of objects)

=cut


sub rel_coords
{
    my ($self, %args) = @_;
    return $self->_get_coords(%args, abs => 0);
}


sub _get_coords
{
    my ($self, %args) = @_;
    my $pos;
    if($args{'pos'}) {
	$pos = $args{'pos'};
    }
    elsif($args{'obj'}) {
	$pos = $self->_pos_by_obj($args{'obj'});
    }
    else {
	croak 'No pos or obj argument';
    }
    croak "Invalid position $pos" if ($pos < 1 or $pos > $self->nr_coords);
    my @c = (@{$self->{'_coords'}->[$pos-1]})[1,2];
    @c = $self->{'_coords'}->[$pos-1]->[0]->rel2abs(@c) if ($args{'abs'});
    return @c;
}


sub _pos_by_obj
{
    my ($self, $obj) = @_;
    for(my $i = 0; $i < $self->nr_coords; $i++) {
	return $i+1 if($self->{'_coords'}->[$i]->[0] eq $obj);
    }
    croak "No coords for object $obj";
}


sub debug_dump
{
    my ($self) = @_;
    my $str = "";
    foreach my $se_pair (@{$self->{'_coords'}}) {
	$str .= $se_pair->[0]->id_loc_str."\t".$se_pair->[1]."-".$se_pair->[2]."\n";
    }
    return $str;
}

