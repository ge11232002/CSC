# AT::Prediction::AntisensePair module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::Prediction::AntisensePair

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::Prediction::AntisensePair::Mapping2Mapping;

use vars '@ISA';
use strict;
use AT::Root;
use AT::Tools::RangeHandler;
use Carp;

@ISA = qw/AT::Prediction::AntisensePair/;


sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
	member1 => ($args{'member1'} || croak 'No member1 argument'),
	member2 => ($args{'member2'} || croak 'No member2 argument'),
	#acc1 => ($args{'acc1'} || croak "No acc1 argument"),
	#acc2 => ($args{'acc2'} || croak "No acc2 argument"),
    }, ref $caller || $caller;

#    # check that there is overlap
#    if($self->member1->strand == $self->member2->strand or !($self->ol_region)) {
#	warn "Attempt to create AntisensePair from non-overlapping genes ".
#	    $self->member1->loc_str. " and ". $self->member2->loc_str;	
#	return undef;
#    }  

    return $self;
}


sub category
{
    my ($self) = @_;

    unless ($self->{category}) {
	my $cat = 1;
	$cat += 2 unless($self->max_exon_overlap >= 20);
	#print "nr introns: ", $self->member1->nr_introns, " ",
	#$self->member2->nr_introns, "\n";
	#_debug_print_gs($self->member1); print "\n";
	#_debug_print_gs($self->member2); print "\n";
	$cat += 1 if($self->member1->is_spliced and
		     $self->member2->is_spliced);
	if($cat == 4) {
	    my @m;
	    if($self->member1->start > $self->member2->start and
	       $self->member1->end < $self->member2->end) {
		@m = ($self->member1, $self->member2);
	    }
	    elsif($self->member2->start > $self->member1->start and
		  $self->member2->end < $self->member1->end) {
		@m = ($self->member2, $self->member1);
	      }
	    if(@m) {
		my @exon_list = $m[1]->HSP_list;
		foreach my $exon (@exon_list[1..@exon_list-2]) {
		    if($m[0]->start < $exon->start and
		       $m[0]->end > $exon->end) {
			    $cat = 5;
			    last;
			}
		}
	    }
	    else {
		$cat = 5;
	    }
	}
	$self->{category} = $cat;
    }
    return $self->{category};
}


sub overlapper_ids
{
    my ($self, $member_nr) = @_;
    die "overlapper_ids() undefined";
}


sub overlapper_gfmappings
{
    my ($self, $member_nr) = @_;
    die "overlapper_gfmappings() undefined";
}


sub representative_ids
{
    my ($self) = @_;
    die "representative_ids() undefined";
}


sub representative_gfmappings
{
    my ($self) = @_;
    die "representative_gfmappings() undefined";
}


sub query_seq_overlap_string
{
    my ($self, $acc, $member_nr) = @_;
    my $gfmap = $member_nr==1 ? $self->member1 : $self->member2;
    my $m;
    foreach my $primap ($gfmap->primary_mRNA_mapping_list) {
	if($primap->qName eq $acc) {
	    $m = $primap;  
	    last;
	}
    }
    unless($m) {
        warn "No primary mapping of $acc found in member $member_nr\n";
        return "";
    }
    my @trseq_ol_list;
    foreach my $genome_ol ($self->exon_overlap_regions) {
	my ($tStart, $tEnd) = @$genome_ol;
	my ($qStart, $qEnd) = $m->qRange_from_tRange($tStart, $tEnd);
	push @trseq_ol_list, "$qStart..$qEnd";
    }
    return join(",",@trseq_ol_list);
}


sub _calc_exon_overlap_regions
{
    my ($self) = @_;
    my @e1 = $self->member1->HSP_list;
    my @e2 = $self->member2->HSP_list;
    my @regions;

    my($i, $j) = (0, 0);

    while($i < @e1 and $j < @e2) {
	my $ol_start = ($e1[$i]->start > $e2[$j]->start) ?
	    $e1[$i]->start : $e2[$j]->start;
	my $ol_end = ($e1[$i]->end < $e2[$j]->end) ?
	    $e1[$i]->end : $e2[$j]->end;
	push @regions, [$ol_start, $ol_end] if($ol_start <= $ol_end);
	if ($e1[$i]->end < $e2[$j]->end) {
	    $i++;
	}
	else {
	    $j++;
	}
    }

    return \@regions;
}


1;
