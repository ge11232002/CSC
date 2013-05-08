# AT::Prediction::SpliceForm module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::Prediction::SpliceForm - represents a distinct spliceform

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::Prediction::SpliceForm;

use AT::Prediction::ExonHolder;

@ISA = qw(AT::Prediction::ExonHolder);


sub new
{
    my ($caller, %args) = @_;

    my ($self) = $caller->SUPER::new(%args);

    #$self->{TU} = $args{TU};
    $self->{_mappingL} = $args{mapping_list} || [];
       # ^^^ should perhaps be moved to superclass
    $self->{_hsp_clusterL} = $args{hsp_cluster_list} || [];

    return $self;
}


#sub mappings { return @{$_[0]->_mappings}; }
