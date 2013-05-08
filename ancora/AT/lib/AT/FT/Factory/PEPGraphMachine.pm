# AT::FT::Factory::PEPGraphMachine module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::FT::Factory::PEPGraphMachine module

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::FT::Factory::PEPGraphMachine;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use AT::FT::PEPGraph;
use AT::FT::PEP;
use AT::FT::PEPString;
use AT::FT::PEPSupport;


#use Class::Struct '_pep_string' => [
#    left_pep => '_pep',
#    right_pep => '_pep',
#    score => '$'
#];
#
#use Class::Struct '_pep_support' => [
#    mapping => 'AT::FT::GFMapping',
#    HSP_nr => '$'
#];
#
#use Class::Struct '_pep' => [
#    start => '$',
#    end => '$',
#    score => '$',
#    prev => '_pep',
#    next => '_pep',
#    staple_score => '$',      # currently not used?
#    left_strings => '%',
#    right_strings => '%',
#    open_support => '@',
#    close_support => '@',
#    visited => '$'
#];

use Class::Struct '_pos' => [
    pos => '$',
    open => '$',
    score1 => '$', # primary score
    score2 => '$',  # secondary score
    mapping => 'AT::FT::GFMapping',
    HSP_nr => '$'
];


use constant MIN_INTRON_LEN_FOR_SS_CHECK => scalar 12;


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
	canonical_jnc_score => ($args{canonical_jnc_score} || 10),
	minor_jnc_score => ($args{minor_jnc_score} || 0),
	other_jnc_score => ($args{other_jnc_score} || 0),
	break_score => ($args{break_score} || 0),
	short_intron_score => ($args{short_intron_score} || 0),
	short_intron_threshold_minor_jnc => ($args{short_intron_threshold_minor_jnc} || 100000),
	short_intron_threshold_canonical_jnc => ($args{short_intron_threshold_canonical_jnc} || 100000),
	add_map_score_to_string => ($args{add_map_score_to_string} || 1),
    }, ref $caller || $caller;
   
    return $self;
}


sub run
{
    my ($self, %args) = @_;
    my $mappings = $args{mappings} || croak "No mappings arg";
    return () unless(@$mappings);
    
    #$self->{_mapping_scores} = {};
    my $pep_graph = $self->_make_pep_graph($mappings);  
    #$self->{_mapping_scores} = undef;

    return $pep_graph;
}


sub _make_pep_graph
{
    my ($self, $group) = @_;

    # Create '_pos' objects that will help us create the graph
    # Each object is a start or end position of a HSP object
    # (the open attribute is true for starts, false for ends)
    # We sort the positions in increasing order
    my @pos;
    foreach my $mapping (@$group) {
	for my $i (1..$mapping->nr_HSPs) {
	    my $HSP = $mapping->HSP($i);
	    #my ($score1, $score2) = $self->_mapping_scores($mapping);
	    my $score1 = $mapping->score1;
	    my $score2 = $mapping->score2;
	    push @pos, _pos->new(pos => $HSP->start,
				 mapping => $mapping,
				 HSP_nr => $i,
				 open => 1,
				 score1 => $score1,
				 score2 => $score2);
	    push @pos, _pos->new(pos => $HSP->end,
				 mapping => $mapping,
				 HSP_nr => $i,
				 open => 0,
				 score1 => -$score1,
				 score2 => -$score2);
	}
    }
    @pos = sort {$a->pos <=> $b->pos or $b->open <=> $a->open} @pos;

    # Create linked list of peps by looking at each pair of adjacent positions
    # Identical positions are treated together.
    # i..j-1 is a range of identical positions in @pos to be used for a pep start
    # j..k-1 is a range of identical positions in @pos to be used for a pep end
    my $pep_head = AT::FT::PEP->new(start => 0, end => 0, score => 0);
	# ^^ dummy head for linked list
    my $prev_pep = $pep_head;
    my %HSP_opening_pep;
    my %supporting_mappings;
    my ($score1_change, $score2_change);
    my $score1 = $pos[0]->score1;
    my $score2 = $pos[0]->score2;
    my ($i, $j) = (0,1);
    for(; $j < @pos and $pos[$i]->pos == $pos[$j]->pos and $pos[$j]->open; $j++) {
	    $score1 += $pos[$j]->score1;
	    $score2 += $pos[$j]->score2;
    }
    while($j < @pos) {

	# update list of supporting mappings
	if($pos[$i]->open) {
	    for(my $l = $i; $l < $j; $l++) {
		my $m = $pos[$l]->mapping;
		$supporting_mappings{$m} = $m;
	    }
	}
	else {
	    for(my $l = $i; $l < $j; $l++) {
		delete $supporting_mappings{$pos[$l]->mapping};
	    }
	}

	# find next range; calc new score change
	$score1_change = $pos[$j]->score1;
	$score2_change = $pos[$j]->score2;
	my $k = $j+1;
	for(; $k < @pos and $pos[$j]->pos == $pos[$k]->pos and
	    $pos[$j]->open == $pos[$k]->open; $k++)
	{
	    $score1_change += $pos[$k]->score1;
	    $score2_change += $pos[$k]->score2;
	}

# attach list of other supporting mappings
# we shouldnt create a support object for each mapping & pep -- too much data
# for now just attach mapping list
# maintain currently supporting mappings as hash %sm
#  key = ref, value = ref
#  update hash when score is changed

#print STDERR $pos[$i]->pos, '-', $pos[$j]->pos, "\t", $pos[$i]->open, ' ',
#$pos[$j]->open, "\t", $score, "\t", $score_change, "\n"
#    if($pos[$j]->pos <= 191749952);

	# Create pep if
	#  - this region has a score (i.e. it is in at least one HSP)
	#  - the positions are not a closing immediately followed by an opening
	if($score1 and ($pos[$i]->open or !$pos[$j]->open or
		       $pos[$i]->pos+1 != $pos[$j]->pos))
	{
	    # create pep for this region
	    my $pep = AT::FT::PEP->new(start => $pos[$i]->pos + ($pos[$i]->open ? 0 : 1),
				end => $pos[$j]->pos - ($pos[$j]->open ? 1 : 0),
				score1 => $score1, score2 => $score2);

	    # here: need to tell peps which HSPs support them

	    # tell peps which HSPs they open & close
	    my (@open_support, @close_support);
	    @open_support =
		map { AT::FT::PEPSupport->new(mapping => $_->mapping,
					      HSP_nr => $_->HSP_nr,
				              score => $_->score1) }
		@pos[$i..$j-1] if($pos[$i]->open);
	    @close_support =
		map { AT::FT::PEPSupport->new(mapping => $_->mapping,
				              HSP_nr => $_->HSP_nr,
				              score => -$_->score1) }
		@pos[$j..$k-1] if(!$pos[$j]->open);
	    $pep->set_open_support(\@open_support);
	    $pep->set_close_support(\@close_support);
	    $pep->set_mapping_list([values %supporting_mappings]);

	    # set open and close scores (we could do this later)
	    $self->_set_pep_open_close_scores($pep);

	    # Associate this pep with each HSP opened by it
	    if($pos[$i]->open) {
		foreach my $pos (@pos[$i..$j-1]) {
		    $HSP_opening_pep{$pos->mapping->HSP($pos->HSP_nr)} = $pep;
		}
	    }

	    # add pep to doubly linked list
	    $pep->prev($prev_pep);
	    $prev_pep->next($pep);
	    $prev_pep = $pep;
	}

	# move range fwd and update score
	($i, $j) = ($j, $k);
	$score1 += $score1_change;
	$score2 += $score2_change;
    }

    # should make support an object consisting of mapping ref, hsp nr

    # Now add 'strings' (intronic connections) between peps
    #  for each pep: collect all supporting left HSPs
    #  they will give 0-more right HSPs
    #  each right HSP will give one right pep
    #  for each left-right connection: increase the score
    my $add_map_score_to_string = $self->add_map_score_to_string;
    for(my $left_pep = $pep_head; $left_pep; $left_pep = $left_pep->next) {
	foreach my $support ($left_pep->close_support_list) {
	    my $mapping = $support->mapping;
	    my $left_HSP = $mapping->HSP($support->HSP_nr);
	    my $right_HSP = $mapping->HSP($support->HSP_nr + 1);
	    next unless ($right_HSP);
	    my $score = $add_map_score_to_string ? $mapping->score1 : 0;
	    my $right_pep = $HSP_opening_pep{$right_HSP};
	    unless($right_pep) {
		if($score) { warn "No right PEP where score=$score"; }
		next;
	    }
	    my $string = $left_pep->right_strings->{$right_pep};
	    unless($string) {
		$score += $self->_string_score($mapping, $left_HSP, $right_HSP);
		$left_pep->right_strings->{$right_pep} =
		$right_pep->left_strings->{$left_pep} =
		    AT::FT::PEPString->new(left_pep => $left_pep,
				     right_pep => $right_pep,
				     score => $score); 
	    }
	    else {
		$string->score($string->score + $score);
	    }
	}
    }

    return AT::FT::PEPGraph->new(head => $pep_head);
}


sub _set_pep_open_close_scores
{
    my ($self, $pep) = @_;
    my ($open_score, $close_score) = (0,0);
    my (@open_qNames, @close_qNames);  # for debug output
    foreach my $support ($pep->open_support_list) {
        if($support->mapping->start == $pep->start and
           $support->mapping->trust_start) {
	    $open_score += $support->score;
	   push @open_qNames, $support->mapping->qName_list;
	}
    }
    foreach my $support ($pep->close_support_list) {
        if($support->mapping->end == $pep->end and
    	$support->mapping->trust_end) {
	    $close_score += $support->score;
	    push @close_qNames, $support->mapping->qName_list;
	}
    }
    $pep->open_score($open_score);
    $pep->close_score($close_score);
    $pep->{_open_qNames} = \@open_qNames; # for dbg
    $pep->{_close_qNames} = \@close_qNames; # for dbg
}


sub _string_score {
    my ($self, $m, $left_hsp, $right_hsp) = @_;
    my $score = 0;
    my $gap = $left_hsp->right_gap;
    my $jnc_type = $gap->jnc_type;
    if($jnc_type eq 'canonical') {
        $score += $self->canonical_jnc_score;
	$score += $self->short_intron_score if($right_hsp->start - $left_hsp->end < $self->short_intron_threshold_canonical_jnc);
    }
    elsif($jnc_type eq 'minor') {
	$score += $self->minor_jnc_score;
	$score += $self->short_intron_score if($right_hsp->start - $left_hsp->end < $self->short_intron_threshold_minor_jnc);
    }
    elsif($gap->is_broken) {
	$score += $self->break_score;
    }
    else {
	$score += $self->other_jnc_score;
    }
    return $score;
}


#sub _mapping_scores {
#    my ($self, $m) = @_;
#    my $scores = $self->{_mapping_scores}{$m};
#    unless ($scores) {
#	my $score1;
#        if(my @primary_m = $m->primary_mRNA_mapping_list) {
#	    my ($mRNA_score, $refSeq_score) = (0,0);
#	    foreach my $primary_m (@primary_m) {
#	        if(substr($primary_m->qName,0,3) eq 'NM_') {
#		    $refSeq_score =
#			($primary_m->mRNAInfo->refSeqStatus eq 'Reviewed') ?
#			$self->revRefSeq_score : $self->refSeq_score;
#		}
#		else {
#		    $mRNA_score = $self->_calc_mRNA_mapping_score($m, $primary_m);
#		}
#	    }
#	    $score1 = $mRNA_score + $refSeq_score;
#	}
#	elsif($m->primary_EST_mapping_list) {
#	    $score1 = $self->EST_score;
#	}
#        else {
#	    warn "No recognized qTypes for mapping";
#	    $score1 = 0;
#	}
#	my $score2 = $m->is_stretch_adjacent ? 0 : $score1;
#	$self->{_mapping_scores}{$m} = $scores = [$score1, $score2];
#    }
#    return @$scores;
#}
#
#
#sub _calc_mRNA_mapping_score {
#    my ($self, $gfm, $m) = @_;
#    
#    my $low_score = $self->low_mRNA_score;
#    my $max_offset = 16;
#    my $min_utr = 20;
#
#    $self->_debug_msg("\nmRNA_scorer: ",$m->qName,"\n");
#
#    # check that the mapping is not revcom'd
#    my $strand = $gfm->strand;
#    return $low_score if($strand != $m->strand_numeric);
#    $self->_debug_msg("mRNA_scorer: ori $strand ok\n");
#
#    # check that entire mrna is mapped (except for tail and up to 15nt on either side)
#    my ($a_tail) = $m->query_seq->seq =~ /(a*)$/i;
#    my ($map_start, $map_end) = ($m->qStart, $m->qEnd);
##    if($strand == 1) {
##	$map_start = $m->qStart;
##	$map_end = $m->qEnd;
##    }
##    else {
##	$map_start = $m->qSize - $m->qEnd + 1;
##	$map_end = $m->qSize - $m->qStart + 1;
##    }
#    return $low_score if($map_start > $max_offset or
#			 $map_end < $m->qSize -$max_offset -length($a_tail));
#    $self->_debug_msg("mRNA_scorer: limits $map_start-$map_end ok w length=",$m->qSize," tail=",length($a_tail),"\n");
#  
#    # check that mrna is from a trusted lib or is non-HTC with a good ORF
#    my ($cds_start, $cds_end) = $m->mRNAInfo->cds =~ /^(\d+)\.\.(\d+)$/;     # get ORF limits
#    $self->_debug_msg("mRNA_scorer: ", ($cds_end ? "cds=$cds_start-$cds_end\n" : "no cds\n"));
#    return $low_score
#	unless($self->{good_mRNA_libs}->{$m->mRNAInfo->library} or
#	       ($m->query_seq->division ne 'HTC'
#		and defined($cds_end)
#		and $cds_start >= $map_start + $min_utr
#		and $cds_end <= $map_end - $min_utr
#		and $cds_end - $cds_start >= 299)
#	       );
#    $self->_debug_msg("mRNA_scorer: trusted mRNA\n");
#
#    # check for long target gaps
#    for my $i (2..$m->nr_HSPs) {
#        return $low_score if($m->HSP($i)->qStart - $m->HSP($i-1)->qEnd > 4); # fail if gap of at least 4 bases
#    }
#    $self->_debug_msg("mRNA_scorer: no long target gaps\n");
#
#    # check for unrecognized splice junctions
#    for my $i (1..$gfm->nr_HSPs-1) {
#    my $sep = $gfm->HSP($i)->right_gap;
#    return $low_score if($sep->is_broken or
#			$sep->jnc_type eq 'other' or
#			$sep->jnc_strand ne $strand);
#    }
#    $self->_debug_msg("mRNA_scorer: jnx OK!\n");
#
#    # return high score if passed all tests
#    return $self->high_mRNA_score;
#}



1;
