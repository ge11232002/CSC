# AT::Prediction::MapPartitionerI module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::Prediction::MapPartitionerI module

=head1 DESCRIPTION

This class defines the interface for classes that partition mappings into
AT::Prediction::Gene objects and further into
AT::Prediction::SpliceForm objects.

=head1 APPENDIX

The rest of the documentation details the methods that subclasses
must implement.
Internal methods are usually preceded with a _

=cut

package AT::Prediction::MapPartitionerI;

use Carp;
use AT::Root;

@ISA = qw(AT::Root);


=head2 new

 Title     : new
 Usage     : my $partitioner = AT::Prediction::ConcreteMapPartitioner->new();
 Function  : Constructor
 Returns   : An object of the class requested.
 Args      : No arguments are currently required by the general
             interface definition.

=cut

sub new { croak "Method new not defined"; }


=head2 partition

 Title     : partition
 Usage     : my @genes = $partitioner->partition(mappings => \@mappings,
			             	         target_seq => $seq);
 Function  : 
 Returns   : An array of AT::Prediction::Gene objects
 Args      : mappings    Ref to array of AT::Mapping objects
             target_seq	 A sequence spanning all the mappings
			 (Bio::LocatableSeq-compliant)

=cut

sub partition { croak "Method partition not defined"; } 


1;
