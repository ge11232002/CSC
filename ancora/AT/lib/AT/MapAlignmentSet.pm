# AT::MapAlignmentSet module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::MapAlignmentSet - represents a set of cDNA-genomic mappings with
homologous or overlapping genomic regions.

=head1 SYNOPSIS

=over 4

=item * get the AT::MapAlignment objects that make up the set

my @alns = $alnset->mapping_alignments();

=item * get the alignment of genomic sequences

my $aln = $alnset->genomic_alignment();

=item * get all disjunct genomic sequences of the set (overlapping sequences
merged)

my $aln = $alnset->genomic_sequences();

=back

Se AT::MapAlignmentSetFactory for further details.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::MapAlignmentSet;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;


@ISA = qw(AT::Root);


sub new  {
    my ($caller, %args) = @_;
    my $self = bless { _alignments => $args{alignments} || [],
		       _genomic_seqs => $args{genomic_seqs} || [],
		       _genomic_aln => $args{genomic_aln},
		   }, ref $caller || $caller;
    return $self;
}


sub mapping_alignments {
    my ($self) = @_;
    return @{$self->_alignments};
}


sub genomic_alignment {
    my ($self) = @_;
    return $self->_genomic_aln;
}


sub genomic_sequences {
    my ($self) = @_;
    return @{$self->_genomic_seqs};
}


1;

