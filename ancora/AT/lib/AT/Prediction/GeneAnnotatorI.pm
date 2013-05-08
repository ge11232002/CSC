# AT::Prediction::GeneAnnotatorI module
#
# Copyright Boris Lenhard
# 
#

# POD

=head1 NAME

AT::Prediction::GeneAnnotatorI module

=head1 SYNOPSIS

=head1 DESCRIPTION

This class defines the interface for classes that annotate features (exons,
introns etc) on predicted genes and spliceforms.

=head1 APPENDIX

The rest of the documentation details the methods that subclasses
must implement.
Internal methods are usually preceded with a _

=cut

package AT::Prediction::GeneAnnotatorI;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;

@ISA = qw(AT::Root);


=head2 new

 Title     : new
 Usage     : my $annotator = AT::Prediction::ConcreteGeneAnnotator->new();
 Function  : Constructor
 Returns   : An object of the class requested.
 Args      : No arguments are currently required by the general
             interface definition.

=cut

sub new { croak "Method new not defined"; }


=head2 annotate_genes

 Title     : annotate_genes
 Usage     : $annotator->annotate_genes(\@genes);
 Function  : Annotate the genes in some way.
 Returns   : -
 Args      : Reference to an array of genes to be annotated.

=cut


sub annotate_genes { croak "Method annotate_genes not defined"; } 


1;
