# AT::FT::PAS module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::FT::PAS

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::FT::PAS;

use strict;
use vars '@ISA';
use AT::Root;
use Bio::SeqFeature::Generic;
use Carp;


@ISA = qw/Bio::SeqFeature::Generic AT::Prediction::Genomic/;


sub new
{
    my ($caller, %args) = @_;

    my $self = $caller->SUPER::new(%args);

    $self->primary_tag($args{'primary'} || $args{'-primary'} || 'PAS');
    $self->{abs_start} = ($args{'abs_start'} || $args{'-abs_start'} ||
			  croak 'No abs_start argument');
    $self->{abs_end} = ($args{'abs_end'} || $args{'-abs_end'} ||
			croak 'No abs_end argument');

    return $self;
}


sub abs_pos
{
    my ($self) = @_;
    return $self->strand == 1 ? $self->abs_start : $self->abs_end;
}


