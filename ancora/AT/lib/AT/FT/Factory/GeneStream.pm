# AT::FT::Factory::GeneStream module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

package AT::FT::Factory::GeneStream;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;

@ISA = qw/AT::Root/;


sub new
{
    my ($caller, %args) = @_;
    my $self = bless {
	_gene_factory => ($args{gene_factory} || croak "No gene_factory arg"),
	_tName => ($args{tName} || croak "No tName arg"),
	_tStarts => ($args{tStarts} || croak "No tStarts arg"),
	#_tEnd => ($args{tEnd} || croak "No tEnd arg"),
	_tStarts_idx => 0
    }, ref($caller) || $caller;
    return $self;
}


# currently this only works for whole chromosomes

sub next_gene_set
{
    my ($self) = @_;
    my $tName = $self->{_tName};        
    my $tStarts = $self->{_tStarts};
    my $i = $self->{_tStarts_idx};
    my $gene_list_ref = [];
    while ($i < @$tStarts and @$gene_list_ref == 0) {
        my $tStart = $tStarts->[$i];
	my $cluster_end;
        ($gene_list_ref, $cluster_end) =
            $self->{_gene_factory}->_create_genes_for_mapping_cluster($tName, undef, $tStart, $tStart+99);
        $i++ while($i < @$tStarts and $tStarts->[$i] <= $cluster_end);
    }
    $self->{_tStarts_idx} = $i;
    return undef unless (@$gene_list_ref);
    return wantarray ? @$gene_list_ref : $gene_list_ref;    
}
