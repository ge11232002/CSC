# AT::FT::Factory::MappingScorer module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::FT::Factory::MappingScorer module

=head1 SYNOPSIS

=head1 APPENDIX

NOTE: THIS MODULE HAVE BEEN REPLACED BY MappingScorer.pm

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::FT::Factory::MappingScorer;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;


@ISA = qw/AT::Root/;


=head2 new

 Title     : new
 Usage     : my $partitioner = AT::Prediction::NaiveMapPartitioner->new();
 Function  : Constructor
 Returns   : AT::Prediction::NaiveMapPartitioner
 Args      : exclute_on_revcom_ss  If true, excludes mappings with
	                           more revcomd canonical than
				   canonical splice sites.

=cut


sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
	EST_score => ($args{EST_score} || 10),
	high_mRNA_score => ($args{high_mRNA_score} || 25),
	low_mRNA_score => ($args{low_mRNA_score} || 22),
	good_mRNA_libs => ($args{good_mRNA_libs} || {}), # do $machine->trusted_mRNA_libs({22 => 1}) to trust mRNAs from lib #22
	refSeq_score => ($args{refSeq_score} || 0),
	revRefSeq_score => ($args{revRefSeq_score} || 10),
	canonical_jnc_score => ($args{canonical_jnc_score} || 10),
	_mapping_scores => undef
    }, ref $caller || $caller;
   
    return $self;
}


sub run
{
    my ($self, %args) = @_;
    my $mappings = $args{mappings} || croak "No mappings arg";
    foreach my $mapping (@$mappings) {
	$self->_score_mapping($mapping);
    }
}


sub _score_mapping {
    my ($self, $m) = @_;

    # calc score1
    my $score1;
    if(my @primary_m = $m->primary_mRNA_mapping_list) {
	my ($mRNA_score, $refSeq_score) = (0,0);
        foreach my $primary_m (@primary_m) {
            if(substr($primary_m->qName,0,3) eq 'NM_') {
        	    $refSeq_score =
		    ($primary_m->mRNAInfo->refSeqStatus eq 'Reviewed') ?
		    $self->revRefSeq_score : $self->refSeq_score;
	    }
	    else {
	        $mRNA_score = $self->_calc_mRNA_mapping_score($m, $primary_m);
	    }
        }
        $score1 = $mRNA_score + $refSeq_score;
    }
    elsif($m->primary_EST_mapping_list) {
        $score1 = $self->EST_score;
    }
    else {
        warn "No recognized qTypes for mapping";
        $score1 = 0;
    }

    # calc score2
    my $score2 = $m->is_stretch_adjacent ? 0 : $score1;

    # set mapping scores
    $m->score1($score1);
    $m->score2($score2);
}


sub _calc_mRNA_mapping_score {
    my ($self, $gfm, $m) = @_;
    
    my $low_score = $self->low_mRNA_score;
    my $max_offset = 16;
    my $min_utr = 20;

    $self->_debug_msg("\nmRNA_scorer: ",$m->qName,"\n");

    # check that the mapping is not revcom'd
    my $strand = $gfm->strand;
    return $low_score if($strand != $m->strand_numeric);
    $self->_debug_msg("mRNA_scorer: ori $strand ok\n");

    # check that entire mrna is mapped (except for tail and up to 15nt on either side)
    my ($a_tail) = $m->query_seq->seq =~ /(a*)$/i;
    my ($map_start, $map_end) = ($m->qStart, $m->qEnd);
    return $low_score if($map_start > $max_offset or
			 $map_end < $m->qSize -$max_offset -length($a_tail));
    $self->_debug_msg("mRNA_scorer: limits $map_start-$map_end ok w length=",$m->qSize," tail=",length($a_tail),"\n");
  
    # check that mrna is from a trusted lib or is non-HTC with a good ORF
    my ($cds_start, $cds_end) = $m->mRNAInfo->cds =~ /^(\d+)\.\.(\d+)$/;     # get ORF limits
    $self->_debug_msg("mRNA_scorer: ", ($cds_end ? "cds=$cds_start-$cds_end\n" : "no cds\n"));
    return $low_score
	unless($self->{good_mRNA_libs}->{$m->mRNAInfo->library} or
	       ($m->query_seq->division ne 'HTC'
		and defined($cds_end)
		and $cds_start >= $map_start + $min_utr
		and $cds_end <= $map_end - $min_utr
		and $cds_end - $cds_start >= 299)
	       );
    $self->_debug_msg("mRNA_scorer: trusted mRNA\n");

    # check for long target gaps
    for my $i (2..$m->nr_HSPs) {
        return $low_score if($m->HSP($i)->qStart - $m->HSP($i-1)->qEnd > 4); # fail if gap of at least 4 bases
    }
    $self->_debug_msg("mRNA_scorer: no long target gaps\n");

    # check for unrecognized splice junctions
    for my $i (1..$gfm->nr_HSPs-1) {
    my $sep = $gfm->HSP($i)->right_gap;
    return $low_score if($sep->is_broken or
			$sep->jnc_type eq 'other' or
			$sep->jnc_strand ne $strand);
    }
    $self->_debug_msg("mRNA_scorer: jnx OK!\n");

    # return high score if passed all tests
    return $self->high_mRNA_score;
}


1;
