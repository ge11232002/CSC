# AT::Prediction::AntisensePair module
#
# Copyright Par Engstrom and Boris Lenhard
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

package AT::Prediction::AntisensePair;

use vars '@ISA';
use strict;
use AT::Root;
use AT::Tools::RangeHandler;
use AT::Prediction::AntisensePair::Mapping2Mapping;
use AT::Prediction::AntisensePair::Mapping2Gene;
use Carp;

@ISA = qw/AT::Root/;


sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
	pair_id => ($args{'pair_id'} || undef),
	member1 => ($args{'member1'} || croak 'No member1 argument'),
	member2 => ($args{'member2'} || croak 'No member2 argument'),
	source_acc1 => $args{'source_acc1'},
	source_acc2 => $args{'source_acc2'},
	source_mem1 => ($args{'source_mem1'} || 0),
	source_mem2 => ($args{'source_mem2'} || 0),
	source_cat => $args{'source_cat'},
	source_overlap_support => ($args{'source_overlap_support'} || []),
	supported_eol_regions => ($args{'supported_eol_regions'} || []),
	total_repeat_exon_ol => ($args{'total_repeat_exon_ol'} || undef),
	notes => ($args{'notes'} || [])
    }, ref $caller || $caller;

    # check that there is overlap
    if($self->member1->strand == $self->member2->strand or !($self->ol_region)) {
	warn "Attempt to create AntisensePair from non-overlapping genes ".
	    $self->member1->loc_str. " and ". $self->member2->loc_str;	
	return undef;
    }  

    return $self;
}


# Here are the category definitions as we understand it:
# Categories 1-2 contain gene pairs with overlap between the exons.
# Category 1: at least one of the genes are intronless
# Category 2: both genes have introns
# Categories 3-5 contain gene pairs without overlap between the exons.
# Category 3: one gene is intronless and entirely contained within an
# intron of the other gene
# Category 4: both genes have introns; one of the genes is entirely
# contained within an intron of the other
# Category 5: both genes have introns; the genes are intertwined (here
# meaning that, for each gene, it is possible to find an exon that is
# within an intron of the other gene)

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
		my @exon_list = $m[1]->exon_list;
		foreach my $exon (@exon_list[1..@exon_list-2]) {
		    if($m[0]->start < $exon->abs_start and
		       $m[0]->end > $exon->abs_end) {
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


#sub category
#{
#    my ($self) = @_;
#
#    unless ($self->{category}) {
#	my $cat = 1;
#	$cat += 2 unless($self->max_exon_overlap);
#	#print "nr introns: ", $self->member1->nr_introns, " ",
#	#$self->member2->nr_introns, "\n";
#	#_debug_print_gs($self->member1); print "\n";
#	#_debug_print_gs($self->member2); print "\n";
#	$cat += 1 if($self->member1->nr_introns and
#		     $self->member2->nr_introns);
#	if($cat == 4) {
#	    $cat++;
#	    my @m;
#	    if($self->member1->start > $self->member2->start and
#	       $self->member1->end < $self->member2->end) {
#		@m = ($self->member1, $self->member2);
#	    }
#	    elsif($self->member2->start > $self->member1->start and
#		  $self->member2->end < $self->member1->end) {
#		@m = ($self->member2, $self->member1);
#	      }
#	    if(@m) {
#		foreach my $intron ($m[1]->intron_list) {   # we may need to rewrite this...
#		    if($m[0]->start >= $intron->abs_start and
#		       $m[0]->end <= $intron->abs_end) {
#			    $cat--;
#			    last;
#			}
#		}
#	    }
#	}
#	$self->{category} = $cat;
#    }
#    return $self->{category};
#}


sub SU_category
{
    my ($self) = @_;
    my $SU = $self->{SU_category};
    return $SU if($SU);
    $SU="";
    foreach my $g ($self->member1, $self->member2) {
	if($g->is_spliced == 0) { $SU .= 'U'; }
	else { $SU .= 'S'; }
    }
    $self->{SU_category} = $SU;
    return $SU;
# ^^ maybe we should require canonical spl jnc?
}


sub relative_orientation
{
    my ($self) = @_;

    # Let a be the gene on the plus strand; let b be the other gene
    my ($a, $b) = ($self->member1, $self->member2);
    ($a, $b) = ($b, $a) if($a->start > $b->start);

    # If a and b start at the same pos: return "F" (full overlap)
    return 'F' if($a->start == $b->start);

    # If b ends before or at the same pos as a: return "F" (full overlap)
    return 'F' if($b->end <= $a->end);

    # If a is on the plus strand: return "C" (convergent)
    return 'C' if($a->strand == 1);

    # Otherwise: return "D" (divergent)
    return 'D';
}


sub set_repeat_exon_ol
{
    my ($self, $genome_db) = @_;
    my @ol_reg = $self->exon_overlap_regions;
    my $rep_bases = 0;
    foreach my $reg (@ol_reg) {
      my $seq = $genome_db->get_genome_seq(chr => $self->member1->chr,
					   start => $reg->[0],
					   end => $reg->[1]);
      $rep_bases += scalar(my @ary = ($seq->seq =~ /[agct]/g));
    }
    return $self->{total_repeat_exon_ol} = $rep_bases;
}


sub relies_on_stretch_adjacent
{
    my ($self) = @_;
    
# one way:
# go through peps over overlap region for each gene ... messy
# another way:
# recreate genes without stretch adjacent with get_gene_wo_stretch_adjacent
# note: needs dynamic genes for this
# see if they overlap...
#  another way would be to move this to genomic1.pl (have it make 4 genes from each dynamic)
#  make the new genes in the output_pairs method
#  probably easier!

# now: test $gene->use_stretch_adjacent(0) for some good region in m2gprocess

}


sub add_source_ortho_support
{
    my ($self, $other) = @_;
    if($self->member1->loc_str eq $other->member1->loc_str) {
	$self->{source_mem1} |= $other->{source_mem1};
	if($self->member2->loc_str eq $other->member2->loc_str) {
	    $self->{source_mem2} |= $other->{source_mem2};
	}
    }
    elsif($self->member1->loc_str eq $other->member2->loc_str) {
	$self->{source_mem1} |= $other->{source_mem2};
	if($self->member2->loc_str eq $other->member1->loc_str) {
	    $self->{source_mem2} |= $other->{source_mem1};
	}
    }
    return 1;
}


sub add_supported_eol_regions
{
    my ($self, $r) = @_;
    $self->{supported_eol_regions} = AT::Tools::RangeHandler->compute_union
	($self->supported_eol_regions, $r);
}


sub add_note
{
    my ($self, $note) = @_;
    push @{$self->{notes}}, $note;
}

#sub overlapper_ids
#{
#    my ($self, $member_nr) = @_;
#    my ($geneA, $geneB);
#    if($member_nr == 1) {
#	$geneA = $self->member1;
#	$geneB = $self->member2;
#    }
#    else {
#	$geneA = $self->member2;
#	$geneB = $self->member1;
#    }
#    my (@mRNA_ids, @EST_ids);
#    my @mappings = $geneA->gfmapping_list;
#    foreach my $m (@mappings) {
#	my $supp_start = $m->start > $geneA->start ? $m->start : $geneA->start;
#	my $ol_start = ($supp_start > $geneB->start) ? $supp_start : $geneB->start;
#	my $supp_end = $m->end < $geneA->end ? $m->end : $geneA->end;
#	my $ol_end = ($supp_end < $geneB->end) ? $supp_end : $geneB->end;
#	my $ol = $ol_end - $ol_start + 1;
#	if($ol > 0) {
#	    my @my_mRNA_ids = $m->mRNA_acc_list;
#	    my @my_EST_ids = $m->EST_acc_list;   
#	    push @mRNA_ids, { id => join(',',@my_mRNA_ids), ol => $ol } if(@my_mRNA_ids);
#	    push @EST_ids, { id => join(',',@my_EST_ids), ol => $ol } if(@my_EST_ids);
#	}
#    }
#    @mRNA_ids = sort {$b->{ol} <=> $a->{ol}} @mRNA_ids;
#    @EST_ids = sort {$b->{ol} <=> $a->{ol}} @EST_ids;
#    return ( [map {$_->{id}} @mRNA_ids], [map {$_->{id}} @EST_ids] );
#}


sub overlapper_ids
{
    my ($self, $member_nr) = @_;
    my $ol_ids = $self->{"overlapper_ids$member_nr"};
    return @$ol_ids if($ol_ids);
    my ($mrna_mappings, $est_mappings) = $self->overlapper_gfmappings($member_nr);
    my @mRNA_ids = map { $_->mRNA_acc_list } @$mrna_mappings;
    my @EST_ids = map { join('/',$_->EST_acc_list) || () } (@$mrna_mappings, @$est_mappings);
    $ol_ids = $self->{"overlapper_ids$member_nr"} = [\@mRNA_ids, \@EST_ids];
    return @$ol_ids;
}


sub overlapper_gfmappings
{
    my ($self, $member_nr) = @_;
    my $ol_gfmap = $self->{"overlapper_gfmappings$member_nr"};
    return @$ol_gfmap if($ol_gfmap);
    my ($geneA, $geneB);
    if($member_nr == 1) {
	$geneA = $self->member1;
	$geneB = $self->member2;
    }
    else {
	$geneA = $self->member2;
	$geneB = $self->member1;
    }
    my @overlapper_mrna_mappings;
    my @overlapper_est_mappings;
    my @all_mappings = $geneA->gfmapping_list;
    foreach my $m (@all_mappings) {
	my $supp_start = $m->start > $geneA->start ? $m->start : $geneA->start;
	my $ol_start = ($supp_start > $geneB->start) ? $supp_start : $geneB->start;
	my $supp_end = $m->end < $geneA->end ? $m->end : $geneA->end;
	my $ol_end = ($supp_end < $geneB->end) ? $supp_end : $geneB->end;
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
    $ol_gfmap = $self->{"overlapper_gfmappings$member_nr"} = [\@overlapper_mrna_mappings, \@overlapper_est_mappings];
    #print STDERR "OL: $member_nr ", scalar(@overlapper_mrna_mappings), " ", scalar(@overlapper_est_mappings), "\n";
    return @$ol_gfmap;
}


sub representative_ids
{
    my ($self) = @_;
    my ($rep_map1, $rep_map2) = $self->representative_gfmappings();
    return unless($rep_map2);
    my $rep_ids1 = ($rep_map1->refSeq_acc_list)[0] ||
	           ($rep_map1->mRNA_acc_list)[0] ||
	           join('/',$rep_map1->EST_acc_list);
    my $rep_ids2 = ($rep_map2->refSeq_acc_list)[0] ||
	           ($rep_map2->mRNA_acc_list)[0] ||
	           join('/',$rep_map2->EST_acc_list);
    return ($rep_ids1, $rep_ids2);
}


sub representative_gfmappings
{
    my ($self) = @_;
    my $rep_gfmap = $self->{representative_gfmappings} ||= $self->_find_representative_gfmappings();
    return @$rep_gfmap;
}


sub _find_representative_gfmappings
{
    my ($self) = @_;

    my $gene1 = $self->member1;
    my $gene2 = $self->member2;
    my ($mrnas1, $ests1) = $self->overlapper_gfmappings(1);
    my ($mrnas2, $ests2) = $self->overlapper_gfmappings(2);
    my ($rep1A, $rep2A, $olA, $rep1B, $rep2B, $olB);

    ($rep1A, $rep2A, $olA) = $self->_get_rep_maps($mrnas1, $mrnas2);
    return [$rep1A, $rep2A] if($olA >= 20);

    ($rep1A, $rep2A, $olA) = $self->_get_rep_maps($mrnas1, $ests2);
    ($rep1B, $rep2B, $olB) = $self->_get_rep_maps($ests1, $mrnas2);
    if($olA >= $olB) {
	return [$rep1A, $rep2A] if ($olA >= 20);
    }
    elsif($olB >= 20) {
	return [$rep1B, $rep2B];
    }

    ($rep1A, $rep2A, $olA) = $self->_get_rep_maps($ests1, $ests2);
    return [$rep1A, $rep2A] if($olA >= 20);

    warn "No representative gmappings for pair, ", $gene1->loc_str, "/", $gene2->loc_str, "\n";
    return [];
}


sub _get_rep_maps
{
    my ($self, $mappings1, $mappings2, $preferred_id1, $preferred_id2) = @_;
    my $gene1 = $self->member1;
    my $gene2 = $self->member2;
    my $max_ol = -1;
    my $n_preferred = 0;
    my $category = $self->category;
    my ($rep1, $rep2);
    foreach my $m1 (@$mappings1) {
	foreach my $m2 (@$mappings2) {
	    my $start1 = $m1->start > $gene1->start ? $m1->start : $gene1->start;
	    my $start2 = $m2->start > $gene2->start ? $m2->start : $gene2->start;
	    my $end1 = $m1->end < $gene1->end ? $m1->end : $gene1->end;
	    my $end2 = $m2->end < $gene2->end ? $m2->end : $gene2->end;
	    my $ol_start = $start1 > $start2 ? $start1 : $start2;
	    my $ol_end = $end1 < $end2 ? $end1 : $end2;
	    my $ol = $ol_end - $ol_start;
	    if($ol > $max_ol and ($category > 2 or $self->_gfmappings_have_exon_overlap($m1,$m2))) {
		$max_ol = $ol;
	        ($rep1, $rep2) = ($m1, $m2);
	    }
	    elsif($ol == $max_ol and $n_preferred != 2) {
		#...
	    }
	}
    }
    #if($rep1) {
     #   print STDERR "REP: ", $rep1->loc_str, " ", $rep2->loc_str, " ", $max_ol+1, "\n";
    #}
    #else {
	#print STDERR "REP: -\n";
    #}
    return ($rep1, $rep2, $max_ol+1);
}


sub get_mRNA_mRNA_pairs
{
    my ($self) = @_;
    
    my ($mrnas1) = $self->overlapper_gfmappings(1);
    my ($mrnas2) = $self->overlapper_gfmappings(2);
    my @pairs;
    foreach my $m1 (@$mrnas1) {
	foreach my $m2 (@$mrnas2) {
	    if($m1->start <= $m2->end and $m1->end >= $m2->start) {
		push @pairs, AT::Prediction::AntisensePair::Mapping2Mapping->new(member1 => $m1, member2 => $m2);
	    }
	}
    }
    return @pairs;
}


sub get_mRNA_gene_pairs
{
    my ($self) = @_;
    my ($mrnas1) = $self->overlapper_gfmappings(1);
    my ($mrnas2) = $self->overlapper_gfmappings(2);
    my @pairs;
    foreach my $m (@$mrnas1) {
	    push @pairs, AT::Prediction::AntisensePair::Mapping2Gene->new(member1 => $m, member2 => $self->member2);
    }
    foreach my $m (@$mrnas2) {
	    push @pairs, AT::Prediction::AntisensePair::Mapping2Gene->new(member2 => $self->member1, member1 => $m);
    }
    return @pairs;
}


sub get_mRNA_mRNA_pairs_fallback_to_mRNA_gene_pairs
{
    my ($self) = @_;
    
    my ($mrnas1) = $self->overlapper_gfmappings(1);
    my ($mrnas2) = $self->overlapper_gfmappings(2);
    my @pairs;
    my %got_partner;
    foreach my $m1 (@$mrnas1) {
	foreach my $m2 (@$mrnas2) {
	    if($m1->start <= $m2->end and $m1->end >= $m2->start) {
		push @pairs, AT::Prediction::AntisensePair::Mapping2Mapping->new(member1 => $m1, member2 => $m2);
		$got_partner{$m1} = 1;
		$got_partner{$m2} = 1;
	    }
	}
    }
    foreach my $m (@$mrnas1) {
	unless($got_partner{$m}) {
	    push @pairs, AT::Prediction::AntisensePair::Mapping2Gene->new(member1 => $m, member2 => $self->member2);
	}
    }
    foreach my $m (@$mrnas2) {
	unless($got_partner{$m}) {
	    push @pairs, AT::Prediction::AntisensePair::Mapping2Gene->new(member2 => $self->member1, member1 => $m);
	}
    }
    return @pairs;
}


sub exon_overlap_support
{
    my ($self) = @_;

    my $result = $self->{exon_overlap_support};
    return @$result if($result);

    unless($self->max_exon_overlap >= 20) {
	$result = $self->{exon_overlap_support} = ["","",0,0,0,0];
	return @$result;
    }

    my $gene1 = $self->member1;
    my $gene2 = $self->member2;
    my $gene_splicing = $self->_pair_splicing_status($self->member1, $self->member2);

    my ($mrnas1, $ests1) = $self->overlapper_gfmappings(1);
    my ($mrnas2, $ests2) = $self->overlapper_gfmappings(2);
    $mrnas1 = [ grep { $self->_gfmapping_has_exon_overlap($_) } @$mrnas1 ];
    $mrnas2 = [ grep { $self->_gfmapping_has_exon_overlap($_) } @$mrnas2 ];
    $ests1 = [ grep { $self->_gfmapping_has_exon_overlap($_) } @$ests1 ];
    $ests2 = [ grep { $self->_gfmapping_has_exon_overlap($_) } @$ests2 ];

    my $splicing = $self->_get_support_from_mapping_sets($mrnas1, $mrnas2, $gene_splicing);
    my $seqtypes = $splicing ? 8 : 0;
    unless($splicing & $gene_splicing) {
	my $rv;
	$rv = $self->_get_support_from_mapping_sets($mrnas1, $ests2, $gene_splicing);
	$seqtypes |= 4 if($rv);
	$splicing |= $rv;
    	$rv = $self->_get_support_from_mapping_sets($ests1, $mrnas2, $gene_splicing);
	$seqtypes |= 2 if($rv);
	$splicing |= $rv;
        unless($splicing & $gene_splicing) {
	    $rv = $self->_get_support_from_mapping_sets($ests1, $ests2, $gene_splicing);
	    $seqtypes |= 1 if($rv);
	    $splicing |= $rv;   
	}
    }  

    if($splicing & 8)           { $splicing = "SS"; }
    elsif(($splicing & 6) == 6) { $splicing = "SU/US"; }
    elsif($splicing & 4)        { $splicing = "SU"; }
    elsif($splicing & 2)        { $splicing = "US"; }
    elsif($splicing & 1)        { $splicing = "UU"; }
    else                        { $splicing = ""; }

    if($seqtypes & 8)           { $seqtypes = "MM"; }
    elsif(($seqtypes & 6) == 6) { $seqtypes = "ME/EM"; }
    elsif($seqtypes & 4)        { $seqtypes = "ME"; }
    elsif($seqtypes & 2)        { $seqtypes = "EM"; }
    elsif($seqtypes & 1)        { $seqtypes = "EE"; }
    else                        { $seqtypes = ""; }

    $result = $self->{exon_overlap_support} = [$seqtypes, $splicing, scalar @$mrnas1, scalar @$mrnas2, scalar @$ests1, scalar @$ests2];
    return @$result;
}


sub _get_support_from_mapping_sets
{
    my ($self, $set1, $set2, $gene_splicing) = @_;

    my $splicing = 0;
   
    foreach my $m1 (@$set1) {
	foreach my $m2 (@$set2) {
	    if($m1->start <= $m2->end and $m1->end >= $m2->start) {
		if($self->_gfmappings_have_exon_overlap($m1, $m2)) {
		    $splicing |= $self->_pair_splicing_status($m1, $m2);
		    return $splicing if($splicing & $gene_splicing);
		}
	    }
	}
    }

    return $splicing;
}


sub _pair_splicing_status
{
    my ($self, $a, $b) = @_;
    
    my $a_stat = $a->is_spliced;
    my $b_stat = $b->is_spliced;
    if($a_stat) {
	return $b_stat ? 8 : 4;
    }
    else {
	return $b_stat ? 2 : 1;
    }
}



sub ol_region
{
    my ($self) = @_;

    my $ol_start = ($self->member1->start > $self->member2->start) ?
	$self->member1->start : $self->member2->start;
    my $ol_end = ($self->member1->end < $self->member2->end) ?
	$self->member1->end : $self->member2->end;
    return ($ol_start <= $ol_end) ? ($ol_start, $ol_end) : undef;
}


sub source_ortho_support_level
{
    my ($self) = @_;
    my $level = 0; 
    $level++ if($self->source_mem1);
    $level++ if($self->source_mem2);
    if($level == 2) {
	$level++
	    if(($self->source_mem1 | $self->source_mem2) == 3);
    }
    return $level;
}


sub _debug_print_gs
{
    my $tu = shift;

    print "GS for ", $tu->loc_str, "\n";
    for(my $i = 1; $i <= $tu->nr_exons or $i <= $tu->nr_introns; $i++) {
	if(my $exon = $tu->exon($i)) {
	    print "E ",$exon->start, "-",$exon->end, "\t";
	    print $exon->abs_start, "-",$exon->abs_end, "\t";
	    print "(", $exon->abs_end-$exon->abs_start+1, ")\n";
	}
	if(my $intron = $tu->intron($i)) {
	    print "I ",$intron->start, "-",$intron->end, "\t";
	    print $intron->abs_start, "-",$intron->abs_end, "\t";
	    print "(", $intron->abs_end-$intron->abs_start+1, ")\n";
	}
    }
}


sub overlap_pos
{
    my ($self) = @_;
    my $pos;
    unless ($pos = $self->{overlap_pos}) {
        my $a = $self->member1;
        my $b = $self->member2;
        my $start = $a->start > $b->start ? $a->start : $b->start; 
        my $end = $a->end < $b->end ? $a->end : $b->end;
	$pos = $self->{overlap_pos} = [$start,$end];
    }
    return @$pos;
}


sub overlap_length
{
    my ($start, $end) = shift->overlap_pos;
    return $end-$start+1;
}


# Return total # bases of exon overlap
sub exon_overlap_sum
{
    my ($self) = @_;

    unless ($self->{exon_overlap_sum}) {
	my $sum = 0;
	foreach my $overlap ($self->exon_overlaps) {
	    $sum += $overlap;
	}
	$self->{exon_overlap_sum} = $sum;
    }
    return $self->{exon_overlap_sum};
}


# Return largest exon overlap
sub max_exon_overlap
{
    my ($self) = @_;
    my $max_ol = 0;
    foreach my $ol ($self->exon_overlaps) {
	$max_ol = $ol if ($ol > $max_ol);
    }
    return $max_ol;
}


sub exon_overlaps
{
    my $self = shift;
    my $ol = $self->{exon_overlaps} ||= [ map { $_->[1] - $_->[0] + 1 } $self->exon_overlap_regions ];
    return @$ol;
}


sub exon_overlap_starts
{
    my $self = shift;
    return map { $_->[0] } $self->exon_overlap_regions;
}


sub exon_overlap_regions
{
    my $self = shift;
    my $regions = $self->{exon_overlap_regions} ||= $self->_calc_exon_overlap_regions();
    return @$regions;
}


sub _calc_exon_overlap_regions
{
    my ($self) = @_;

    my @e1 = $self->member1->exon_list;
    my @e2 = $self->member2->exon_list;
    my @regions;

    # sort in ascending absolute coord order (i.e. reverse array for - strand)
    if($self->member1->strand == -1) {
	@e1 = reverse @e1;
    }
    else {
	@e2 = reverse @e2;
    }

    my($i, $j) = (0, 0);

    while($i < @e1 and $j < @e2) {
	my $ol_start = ($e1[$i]->abs_start > $e2[$j]->abs_start) ?
	    $e1[$i]->abs_start : $e2[$j]->abs_start;
	my $ol_end = ($e1[$i]->abs_end < $e2[$j]->abs_end) ?
	    $e1[$i]->abs_end : $e2[$j]->abs_end;
	push @regions, [$ol_start, $ol_end] if($ol_start <= $ol_end);
	if ($e1[$i]->abs_end < $e2[$j]->abs_end) {
	    $i++;
	}
	else {
	    $j++;
	}
    }

    return \@regions;
}


sub _gfmappings_have_exon_overlap
{
    my ($self, $a, $b) = @_;
    return unless($a->start <= $b->end and $a->end >= $b->start);
    my @e1 = $a->HSP_list;
    my @e2 = $b->HSP_list;

    my($i, $j) = (0, 0);

    #print STDERR "comparing ", $a->acc_string, " vs ", $b->acc_string, "\n";

    while($i < @e1 and $j < @e2) {
	my $ol_start = ($e1[$i]->start > $e2[$j]->start) ?
	    $e1[$i]->start : $e2[$j]->start;
	my $ol_end = ($e1[$i]->end < $e2[$j]->end) ?
	    $e1[$i]->end : $e2[$j]->end;
	return 1 if($ol_end - $ol_start >= 19);
	if ($e1[$i]->end < $e2[$j]->end) {
	    $i++;
	}
	else {
	    $j++;
	}
    }
    #print STDERR "no overlap\n";
    return 0;
}


sub _gfmapping_has_exon_overlap
{
    my ($self, $m) = @_;
    my @e2 = $self->exon_overlap_regions;
    return unless(@e2 and $m->start <= $e2[-1][1] and $m->end >= $e2[0][0]);
    my @e1 = $m->HSP_list;

    my($i, $j) = (0, 0);

    #print STDERR "comparing ", $a->acc_string, " vs ", $b->acc_string, "\n";

    while($i < @e1 and $j < @e2) {
	my $ol_start = ($e1[$i]->start > $e2[$j][0]) ?
	    $e1[$i]->start : $e2[$j][0];
	my $ol_end = ($e1[$i]->end < $e2[$j][1]) ?
	    $e1[$i]->end : $e2[$j][1];
	return 1 if($ol_end - $ol_start >= 19);
	if ($e1[$i]->end < $e2[$j][1]) {
	    $i++;
	}
	else {
	    $j++;
	}
    }
    #print STDERR "no overlap\n";
    return 0;
}



sub _range_has_exon_overlap
{
    my ($self, $start, $end) = @_;
    my @e2 = $self->exon_overlap_regions;
    return 0 unless(@e2 and $start <= $e2[-1][1] and $end >= $e2[0][0]);

    foreach my $exol (@e2) {
	my $ol_start = ($start > $exol->[0]) ? $start : $exol->[0];
	my $ol_end = ($end < $exol->[1]) ? $end : $exol->[1];
	return 1 if($ol_end - $ol_start >= 19);
	return 0 if($end < $exol->[1]);
    }
    return 0;
}


## Return exon overlaps in decreasing order
## use exon features for each member
## require exactly matching junctions for merging two overlaps
#sub exon_overlaps
#{
#    my ($self) = @_;
#
#    my @e1 = $self->member1->exon_list;
#    my @e2 = $self->member2->exon_list;
#    my @overlaps;
#
#    # sort in ascending absolute coord order (i.e. reverse array for - strand)
#    if($self->member1->strand == -1) {
#	@e1 = reverse @e1;
#    }
#    else {
#	@e2 = reverse @e2;
#    }
#
#    my($i, $j) = (0, 0);
#
#    while($i < @e1 and $j < @e2) {
#	my $ol_start = ($e1[$i]->abs_start > $e2[$j]->abs_start) ?
#	    $e1[$i]->abs_start : $e2[$j]->abs_start;
#	my $ol_end = ($e1[$i]->abs_end < $e2[$j]->abs_end) ?
#	    $e1[$i]->abs_end : $e2[$j]->abs_end;
#	my $ol_len = $ol_end - $ol_start + 1;
#	if ($ol_len > 0) {
#	    if(@overlaps and
#	       $e1[$i-1]->abs_end == $e2[$j-1]->abs_end and
#	       $e1[$i]->abs_start == $e2[$j]->abs_start) {
#		$overlaps[-1] += $ol_len;
#	    }
#	    else {
#		push @overlaps, $ol_len;
#	    }
#	}
#	if ($e1[$i]->abs_end < $e2[$j]->abs_end) {
#	    $i++;
#	}
#	else {
#	    $j++;
#	}
#    }
#
#    return @overlaps;
#}


