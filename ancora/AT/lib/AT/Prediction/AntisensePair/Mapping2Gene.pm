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

package AT::Prediction::AntisensePair::Mapping2Gene;

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
		if($m[1]->isa("AT::FT::GFMapping")) {
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
		    my @exon_list = $m[1]->exon_list;
		    foreach my $exon (@exon_list[1..@exon_list-2]) {
		        if($m[0]->start < $exon->abs_start and
		           $m[0]->end > $exon->abs_end) {
				$cat = 5;
			    last;
			}
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
    my ($mrna_exon_ol_ids, $mrna_nonexon_ol_ids,
	$est_exon_ol_ids, $est_nonexon_ol_ids) = $self->_get_overlapper_ids($member_nr);
    return ([@$mrna_exon_ol_ids, @$mrna_nonexon_ol_ids],
	    [@$est_exon_ol_ids, @$est_nonexon_ol_ids]
	    );
}


sub exon_overlapper_ids
{
    my ($self, $member_nr) = @_;
    my ($mrna_exon_ol_ids, undef, $est_exon_ol_ids) = $self->_get_overlapper_ids($member_nr);
    return ($mrna_exon_ol_ids, $est_exon_ol_ids);
}


sub nonexon_overlapper_ids
{
    my ($self, $member_nr) = @_;
    my (undef, $mrna_nonexon_ol_ids, undef, $est_nonexon_ol_ids) = $self->_get_overlapper_ids($member_nr);
    return ($mrna_nonexon_ol_ids, $est_nonexon_ol_ids);
}


sub _get_overlapper_ids
{
    my ($self, $member_nr) = @_;
    my $ol_ids = $self->{"overlapper_ids$member_nr"};
    return @$ol_ids if($ol_ids);

    my (@mrna_exol_ids, @mrna_nonexol_ids, @est_exol_ids, @est_nonexol_ids);
    my ($mrna_mappings, $est_mappings) = $self->overlapper_gfmappings($member_nr);
    foreach my $m (@$mrna_mappings, @$est_mappings) {
	if($self->_gfmapping_has_exon_overlap($m)) {
	    push @mrna_exol_ids, $m->mRNA_acc_list;
	    if($m->acc_list == 1) {
		push @est_exol_ids, $m->EST_acc_list;
	    }
	    else {
		unless($m->primary_mapping_list) {
		    print STDERR "WARNING: attempt to call _get_overlapper_ids without loading primary mappings\n";
		}
		foreach my $primap ($m->primary_EST_mapping_list) {
		    push @est_exol_ids, $primap->qName
			if($self->_range_has_exon_overlap($primap->tStart, $primap->tEnd));
		}
	    }
	}
	else 
	{
	    push @mrna_nonexol_ids, $m->mRNA_acc_list;
	    push @est_nonexol_ids, $m->EST_acc_list;
	}
    }
    
    return (\@mrna_exol_ids, \@mrna_nonexol_ids, \@est_exol_ids, \@est_nonexol_ids);
}


sub overlapper_gfmappings
{
    my ($self, $member_nr) = @_;
    my $ol_gfmap = $self->{"overlapper_gfmappings$member_nr"};
    return @$ol_gfmap if($ol_gfmap);

    my ($memA, $memB) = $member_nr==1 ? ($self->member1,$self->member2) : ($self->member2,$self->member1);
    my @overlapper_mrna_mappings;
    my @overlapper_est_mappings;

    if($memA->isa("AT::FT::GFMapping")) {
	if($memA->mRNA_acc_list) {
	    push @overlapper_mrna_mappings, $memA;
	}
	else {
	    push @overlapper_est_mappings, $memA;
	}
    }   
    else {
	my @all_mappings = $memA->gfmapping_list;
	foreach my $m (@all_mappings) {
	    my $supp_start = $m->start > $memA->start ? $m->start : $memA->start;
	    my $ol_start = ($supp_start > $memB->start) ? $supp_start : $memB->start;
	    my $supp_end = $m->end < $memA->end ? $m->end : $memA->end;
	    my $ol_end = ($supp_end < $memB->end) ? $supp_end : $memB->end;
	    my $ol = $ol_end - $ol_start + 1;
	    if($ol > 0) {
		if($m->mRNA_acc_list) {
		    push @overlapper_mrna_mappings, $m;
		}
		else {
		    push @overlapper_est_mappings, $m;
		}
	    }
	}
    }
    $ol_gfmap = $self->{"overlapper_gfmappings$member_nr"} = [\@overlapper_mrna_mappings, \@overlapper_est_mappings];
    return @$ol_gfmap;
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
    unless($gfmap->isa("AT::FT::GFMapping")) {
	warn "Member $member_nr is not a GFMapping";
	return "";
    }
    my $m;
    foreach my $primap ($gfmap->primary_mRNA_mapping_list) {
	if($primap->qName eq $acc) {
	    $m = $primap;  
	    last;
	}
    }
    unless($m) {
        warn "No primary mapping of $acc found in member $member_nr";
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
    my $member1 = $self->member1;
    my $member2 = $self->member2;
    my (@e1,@e2);
    if($member1->isa("AT::FT::GFMapping")) {
	@e1 = $member1->HSP_list;
	@e2 = $member2->exon_list;
        @e2 = reverse @e2 if($member2->strand == -1);
    }
    else {
	@e1 = $member2->HSP_list;
	@e2 = $member1->exon_list;
        @e2 = reverse @e2 if($member1->strand == -1);
    }
    my @regions;

    my($i, $j) = (0, 0);

    while($i < @e1 and $j < @e2) {
	my $ol_start = ($e1[$i]->start > $e2[$j]->abs_start) ?
	    $e1[$i]->start : $e2[$j]->abs_start;
	my $ol_end = ($e1[$i]->end < $e2[$j]->abs_end) ?
	    $e1[$i]->end : $e2[$j]->abs_end;
	push @regions, [$ol_start, $ol_end] if($ol_start <= $ol_end);
	if ($e1[$i]->end < $e2[$j]->abs_end) {
	    $i++;
	}
	else {
	    $j++;
	}
    }

    return \@regions;
}



1;
