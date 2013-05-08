# AT::Prediction::UndefGenePart module
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

package AT::Prediction::UndefGenePart;

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

    $self->primary_tag($args{'primary'} || $args{'-primary'} || 'undef_gene_part');
    $self->{abs_start} = ($args{'abs_start'} || $args{'-abs_start'} ||
			  croak 'No abs_start argument');
    $self->{abs_end} = ($args{'abs_end'} || $args{'-abs_end'} ||
			croak 'No abs_end argument');

    return $self;
}


