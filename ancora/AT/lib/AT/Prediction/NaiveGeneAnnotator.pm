# AT::Prediction::NaiveGeneAnnotator module
#
# Copyright Boris Lenhard
# 
#

# POD

=head1 NAME

AT::Prediction::NaiveGeneAnnotator module

=head1 SYNOPSIS

=head1 SYNOPSIS

This is class is meant to provide basic gene structure annotation without
questioning the supporting mappings much.

It is still being developed and thus subject to change.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::Prediction::NaiveGeneAnnotator;

use strict;
use vars '@ISA';
use Carp;
use AT::Prediction::GeneAnnotatorI;
use AT::Prediction::Gene;
use AT::Prediction::SpliceForm;
use AT::Prediction::Exon;
use AT::Prediction::Intron;

use constant DEF_MIN_INTRON_SIZE => scalar 10;

@ISA = qw(AT::Prediction::GeneAnnotatorI);

=head2 new

 Title     : new
 Usage     : my $annotator = AT::Prediction::NaiveGeneAnnotator->new
		(min_intron_size => 6);
 Function  : Constructor
 Returns   : AT::Prediction::NaiveGeneAnnotator
 Args      : min_intron_size    min intron size in bp's (optional)
 
=cut


sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
	min_intron_size => ($args{min_intron_size} || DEF_MIN_INTRON_SIZE)
    }, ref $caller || $caller;
    
    return $self;
}


=head2 annotate_genes

 Title     : annotate_genes
 Usage     : $annotator->annotate_genes(\@genes);
 Function  : Adds exons and introns to the genes in a naive way
	     (just trusts the supporting mappings).
	     Target inserts must exceed a certain length to be
	     annotated as introns (otherwise they are considered
	     as part of an exon).
 Returns   : -
 Args      : Reference to an array of genes to be annotated.

=cut


sub annotate_genes
{
    my ($self, $genes) = @_;

    foreach my $gene (@$genes) {
	$self->_add_exons_introns($gene);
    }
}


sub _add_exons_introns
{
    my ($self, $gene) = @_;

    my @hsp_cl = $gene->hsp_cluster_list;
    my (@exons, @introns);

    my $i = 0;
    while($i < @hsp_cl) {
	my $j = $i+1;
	while($j < @hsp_cl and
	      $hsp_cl[$j-1]->{end} + $self->min_intron_size
	      >= $hsp_cl[$j]->{start})
	{
	    $j++;
	} 
	push @exons, AT::Prediction::Exon->new
	    ( -start => $hsp_cl[$i]->{start} - $gene->start + 1,
	      -end => $hsp_cl[$j-1]->{end} - $gene->start + 1,
	      -strand => $gene->strand,
	      -abs_start => $hsp_cl[$i]->{start},
	      -abs_end => $hsp_cl[$j-1]->{end}
	      );
	if ($j < @hsp_cl) {
	    push @introns, AT::Prediction::Intron->new
		( -start => $hsp_cl[$j-1]->{end}+1 - $gene->start + 1,
		  -end => $hsp_cl[$j]->{start}-1 - $gene->start + 1,
		  -strand => $gene->strand,
		  -abs_start => $hsp_cl[$j-1]->{end}+1,
		  -abs_end => $hsp_cl[$j]->{start}-1
		  );
	}
	$i = $j;	
    }

    $gene->set_exons_introns(\@exons, \@introns);

    # let i be beginning cluster, let j be ending cluster
    # start at i=0; let j = i+1; try to increment j as much as possible
    # make exon from i to j-1
    # make intron from j-1 to j (if j exists)
    # let i=j
}

1;
