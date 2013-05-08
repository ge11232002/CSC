# AT::Seq::Generic module
#
# Copyright Boris Lenhard
# You may distrubute this code under the same terms as perl itself
#
# POD documentation:

=head1 NAME

AT::Seq::Generic

=head2 SYNOPSIS


=head2 DESCRIPTION


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# The code begins HERE

package AT::Seq::Generic;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use Bio::Seq::RichSeq;


@ISA = qw(Bio::Seq::RichSeq AT::Root);


=head2 new

 Title     : new
 Usage     : my $seq = AT::Seq::Generic->new();
 Function  : Constructor
 Returns   : An AT::Seq::Generic object
 Args      :
 
=cut

sub new  {
    my ($caller, %args) = @_;

    my $self = $caller->SUPER::new(%args);

    $self->{'type'} = $args{'type'} || $args{'-type'} || 'unknown';
    $self->{'replaced'} = $args{'replaced'} || 0;
    $self->{'data_source_id'} = $args{'data_source_id'} || 0;
    
    return $self;
}




1;
