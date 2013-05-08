# AT::Prediction::Gene module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::Prediction::Gene - represents a predicted gene

=head1 SYNOPSIS



=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::Prediction::Gene;

use strict;
use vars '@ISA';
use AT::Prediction::ExonHolder;

@ISA = qw(AT::Prediction::ExonHolder);


=head2 new

 Title     : new
 Usage     : my $gene = AT::Prediction::Gene->new
		(seq => $seq,
		 spliceform_list => [ \@spliceforms ] );				
 Function  : Constructor.
	     This method is not intended to be used directly.
	     Use AT::Prediction::GenePredictor to create gene
	     objects.
 Returns   : An AT::Prediction::Gene object
 Args      : 

=cut

sub new
{
    my ($caller, %args) = @_;
    
    my ($self) = $caller->SUPER::new(%args);

    my $seq = $args{seq};
    unless ($seq) {
	die "No sequence provided. Cannot make gene object.";
    }

    $self->{seq} = $seq;
    $self->{_TSSL} = $args{TSS_list} || [];
    $self->{_PASL} = $args{PAS_list} || [];
    $self->{_spliceformL} = $args{spliceform_list} || [];
    $self->{_hsp_clusterL} = $args{hsp_cluster_list} || [];
    $self->{_unique_mappingL} = $args{mapping_list} || $args{unique_mapping_list} || [];
    $self->{_nonunique_mappingL} = $args{nonunique_mapping_list} || [];
    $self->{_excluded_mappingL} = $args{excluded_mapping_list} || [];
    $self->{_commentL} = $args{comment_list} || [];

    return $self;
}
			 

=head2 mapping_list

 Title     : mapping_list
 Usage     : my @mappings = $gene->mapping_list;
 Function  : Returns all the mappings supporting this gene prediction.
 Returns   : Array of AT::Mapping objects
 Args      : None

=cut

sub mapping_list
{
    my ($self) = @_;
    return ($self->unique_mapping_list, $self->nonunique_mapping_list);
}


sub mRNA_mapping_list
{
    my ($self) = @_;
    return grep {$_->qType eq 'mRNA'} $self->mapping_list;
}


sub EST_mapping_list
{
    my ($self) = @_;
    return grep {$_->qType eq 'EST'} $self->mapping_list;
}


sub print_borders
{
    my ($self, $stream) = @_;
    $stream = \*STDOUT unless (defined $stream);
    print $stream "Borders for ", $self->loc_str, "\n";
    foreach my $b (sort {$a->abs_start <=> $b->abs_start}
		   ($self->TSS_list, $self->PAS_list)) {
        my ($rel_start, $rel_end) = $self->abs2rel($b->abs_start, $b->abs_end);
        print $stream (ref($b) =~ /TSS/) ? 'TSS' : 'PAS';
	print $stream " ", $rel_start, "-", $rel_end, "\t";
	print $stream $b->abs_start, "-", $b->abs_end, " ", $b->abs_pos, "\n";
    }
}


=head2 want more methods?

This class inherits AT::Prediction::ExonHolder, which in turn inherits
AT::Prediction::Genomic. See those modules for more method descriptions.

=cut

sub _nr
{
    my ($self, @redundant) = @_;
    my %non_redundant;
    foreach my $item (@redundant) {
	$non_redundant{$item} = $item;
    }
    return values(%non_redundant);
}
