# AT::Prediction::Exon module
#
# Copyright Boris Lenhard
# 
#

# POD

=head1 NAME

AT::Prediction::Exon - abstract class

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::Prediction::Exon;

use strict;
use vars '@ISA';
use AT::Root;
use Bio::SeqFeature::Generic;
use Carp;


@ISA = qw(Bio::SeqFeature::Generic AT::Prediction::Genomic);


sub new
{
    my ($caller, %args) = @_;

    my $self = $caller->SUPER::new(%args);

    $self->primary_tag($args{'primary'} || $args{'-primary'} || 'exon');
    $self->{abs_start} = ($args{'abs_start'} || $args{'-abs_start'} ||
			  croak 'No abs_start argument');
    $self->{abs_end} = ($args{'abs_end'} || $args{'-abs_end'} ||
			croak 'No abs_end argument');

    return $self;
}


