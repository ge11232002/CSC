# AT::Prediction::NaiveMapPartitioner module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::Prediction::MappingEndMover module

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut



package AT::FT::Factory::MappingEndValidator;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use Bio::Location::Fuzzy;

@ISA = qw/AT::Root/;


=head2 new

 Title     : new
 Usage     : 
 Function  : Constructor
 Returns   : 
 Args      : 

=cut


sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
	min_UTR_length => ($args{'min_UTR_length'} || 20),
	verbose => ($args{'verbose'} || 0)
    }, ref $caller || $caller;
    
    return $self;
}


sub run
{
    my ($self, %args) = @_;
    my $mappings = $args{mappings} || croak "No mappings arg";
    #my $target_seq = $args{target_seq} || croak "No target_seq arg";    
    $self->_debug_msg("\n-- VALIDATOR BEGIN --\n");
    foreach my $mapping (@$mappings) {
	$self->_validate_ends($mapping);
    }
    $self->_debug_msg("-- VALIDATOR END --\n\n");
}


# Ends are trusted if:
# query mapped up to 10 bp from tail

sub _validate_ends
{
    my ($self, $gfmapping) = @_;

    if(my ($m) = $gfmapping->primary_mRNA_mapping_list) {
      if($m->strand_numeric == 1) {
	$gfmapping->trust_start(1) if($m->qStart <= $m->qInitTs + 10);
	$gfmapping->trust_end(1) if($m->qEnd >= $m->qSize - 10 - $m->qTermAs);
      }
      else {
	$gfmapping->trust_end(1) if($m->qStart <= $m->qInitTs + 10);
	$gfmapping->trust_start(1) if($m->qEnd >= $m->qSize - 10 - $m->qTermAs);
      }
    }
    else {
      foreach my $m ($gfmapping->primary_EST_mapping_list) {
	if($m->strand_numeric == 1) {
	  $gfmapping->trust_start(1) if($m->tStart <= $gfmapping->start+5 and $m->qStart <= $m->qInitTs + 10);
	  $gfmapping->trust_end(1) if($m->tEnd >= $gfmapping->end-5 and $m->qEnd >= $m->qSize - 10 - $m->qTermAs);
	}
	else {
	  $gfmapping->trust_end(1) if($m->tEnd >= $gfmapping->end-5 and $m->qStart <= $m->qInitTs + 10);
	  $gfmapping->trust_start(1) if($m->tStart <= $gfmapping->start+5 and $m->qEnd >= $m->qSize - 10 - $m->qTermAs);
	}
      }
    }
}


######################
# OLD CODE STARTS HERE
######################

# Ends are trusted if:
# a. the transcript is a reviewed RefSeq
# b. the transcript is an mRNA, and
#     any CDS annotation leaves 20 bp of UTR
# c. the transcript is an EST annotated as sequenced from 5' or 3' end, then
#     the end it was sequenced from is trusted

sub _OLD_validate_ends
{
    my ($self, $m) = @_;
    my $refseq = $m->primary_refSeq_mapping;
    if($refseq and $refseq->mRNAInfo->refSeqStatus eq 'Reviewed') {
	$m->trust_start(1);
	$m->trust_end(1);
    }
    elsif($m->strand) {
	if($m->primary_mRNA_mapping_list) {
	    my ($start_ok, $end_ok) = $self->_check_mrnas($m);
	    if($m->strand == 1) {
		$m->trust_start(1) if($start_ok);
		$m->trust_end(1) if($end_ok);
	    }
	    else {
		$m->trust_start(1) if($end_ok);
		$m->trust_end(1) if($start_ok);
	    }
	}
	elsif($m->strand == 1) {
	    $m->trust_start(1) if($m->has_5p_EST);
	    $m->trust_end(1) if($m->has_3p_EST); # or $tail_ori);
	}
	elsif($m->strand == -1) {
	    $m->trust_start(1) if($m->has_3p_EST); # or $tail_ori);
	    $m->trust_end(1) if($m->has_5p_EST);
	}
	$m->trust_start(0) if($m->query_bp_past_left_border);
	$m->trust_end(0) if($m->query_bp_past_right_border);
    }
#    print STDERR join("\t",
#	($refseq ? $refseq->mRNAInfo->refSeqStatus : '-',
#	 scalar($m->primary_mRNA_mapping_list) ? 'mRNA' : '-',
#	 $m->strand,
#	 $m->has_3p_EST ? '3pEST' : '-',
#	 $m->has_5p_EST ? '5pEST' : '-',
#	 $m->query_bp_past_left_border,
#	 $m->query_bp_past_right_border,
#	 $m->trust_start ? 'S' : '-',
#	 $m->trust_end ? 'E' : '-'),
#	 $m->qName_string), "\n";
}



sub _check_mrnas
{
    my ($self, $m) = @_;
    my $min_utr_len = $self->min_UTR_length;
    my ($start_ok, $end_ok) = (1,1);
    foreach my $mrna ($m->primary_mRNA_mapping_list) {
	my $cds = $mrna->mRNAInfo->cds;
	if($cds ne 'n/a') {
	    my ($cds_start, $cds_end) = ($cds =~ /^(.*)\.\.(.*)$/);
	    unless ($cds_end) {
	        warn "!! Could not parse $cds\n";
	        return (0,0);
	    }
	    if($cds_start =~ /</ or $cds_start < $mrna->qStart + $min_utr_len) {
		if($mrna->strand_numeric == $m->strand) {
		    $start_ok = 0;
		}
		else {
		    $end_ok = 0;
		}
	    }
	    if($cds_end =~ />/ or $cds_end > $mrna->qEnd - $min_utr_len) {
		if($mrna->strand_numeric == $m->strand) {
		    $end_ok = 0;
		}
		else {
		    $start_ok = 0;
		}
	    }
	    $self->_debug_msg($mrna->qName," $cds_start-$cds_end ",$mrna->qStart,"-",$mrna->qEnd);
	}
    }
    $self->_debug_msg("$start_ok $end_ok\n");
    return ($start_ok, $end_ok);
}

1;
