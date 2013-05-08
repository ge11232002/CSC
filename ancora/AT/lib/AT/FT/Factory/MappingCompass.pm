# AT::FT::Factory::MappingCompass module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::FT::Factory::MappingCompass module

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::FT::Factory::MappingCompass;

use strict;
use vars '@ISA';
use Carp;
use AT::Prediction::MapPartitionerI;
use AT::Prediction::Gene;
use AT::Prediction::SpliceForm;
#use AT::Tools::PolyATailFinder;
use AT::Tools::PolyNStretchFinder;

use constant ANNORI_CONFIDENCE_LEVEL => scalar 99;
use constant ANNORI_MIN_PPT_CORRECT => scalar 990;
use constant PROBLEM_TAG_METHOD => scalar 'compass';

@ISA = qw/AT::Prediction::MapPartitionerI/;


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
	min_intron_length => ($args{min_intron_length} || 20),
	evaluation => ($args{evaluation} || 0),
	#exclude_misoriented_mRNAs => ($args{exclude_misoriented_mRNAs} || 0),
	_trusted_lib_lists => {}
    }, ref $caller || $caller;
 
    #$self->verbose(2);
 
    return $self;
}


sub run
{
    my ($self, %args) = @_;
    my $mappings = $args{mappings} || croak "No mappings arg";
    my $target_seq = $args{target_seq} || croak "No target_seq arg";
    my $strand = $args{strand} || 0;

    # add support for use of ref. set

    print STDERR "Compass: orienting ".scalar(@$mappings)." mappings\n" if($self->debug);

    return $self->_orient($mappings, $target_seq);
}


sub output_evaluation
# this should really be a method in an evaluation object
{
    my ($self,$fh) = @_;
    my $eval_ref = $self->evaluation;
    unless($eval_ref) {
	warn "No evaluation";
	return;
    }
    my %e = %$eval_ref;  # make working copy

    my @qTypes = ('RefSeq','mRNA','EST');
    my @counts = ('splOri- testOri-', 'splOri- testOri+', 'splOri+ testOri-', 'correct', 'incorrect');

    print $fh join("\t", ' 'x22, @qTypes), "\n";
    foreach my $count (@counts) {
	print $fh sprintf("%22s",$count);
	foreach my $qType (@qTypes) {
	    my $key = "${qType} $count";
	    my $number = $e{$key} || 0;
	    print $fh "\t",$number;
	    delete $e{$key};
	}
	print "\n";
    }
    
    while(my ($key, $value) = each %e) {
	print $fh $key, ":\t", $value, "\n";
    }
}


sub _orient
{
    my ($self, $mappings, $seq) = @_;
    my $evaluation = $self->evaluation();
    my (@plus, @minus, @ambiv);
    my ($plus_splices, $minus_splices) = $self->_find_introns($mappings, $seq);
    foreach my $mapping (@$mappings) {
	#$self->_set_EST_read_dir($mapping);
	$self->_set_stretch_adjacency($mapping, $seq);
	$self->_set_polyA_signals($mapping, $seq);
	my $ori = 0;
	my ($spljnc_ori) = $self->_orient_by_spljnc($mapping, $seq);
	if($spljnc_ori and $mapping->nr_introns > 1) {
	    $ori = $spljnc_ori;
	}
	if($evaluation) {
	    my $ori2 = $self->_orient_by_tail_signal_and_annotation($mapping, 0, {}, {});
	    $self->_update_evaluation($mapping, $ori, $ori2);
	}
        else {
	    $ori = $self->_orient_by_tail_signal_and_annotation($mapping, $spljnc_ori, $minus_splices, $plus_splices)
		unless($ori);
	    # check internal priming; this really has nothing to do with the orientation
	    if(($ori == 1 and $mapping->ends_near_A_stretch and !$mapping->has_plus_strand_polyA_signal) or
	       ($ori == 2 and $mapping->starts_near_T_stretch and !$mapping->has_minus_strand_polyA_signal)) {
		$mapping->internally_primed(1);
		print STDERR $mapping->acc_string, " internally primed\n" if($self->debug());
	    }
	    # push mappings according to ori	
	    print STDERR "ori: $ori\n" if($self->debug);
	    if($ori == 1) {
		$mapping->strand(1);
		$self->_unset_reverse_introns($mapping) if($mapping->nr_introns);
		push @plus, $mapping;
	    }
	    elsif($ori == 2) {
		$mapping->strand(-1);
		$self->_unset_reverse_introns($mapping) if($mapping->nr_introns);
		push @minus, $mapping;
	    }
	    else {
		$mapping->remove_tag('ori_support') if($mapping->has_tag('ori_support'));
		$mapping->add_problem_tag(PROBLEM_TAG_METHOD, 'no_ori_data') if($ori == 0);
		$mapping->strand(0);
		push @ambiv, $mapping;
	    }
	}
    }
    return (\@plus, \@minus, \@ambiv);
}


sub _orient_by_tail_signal_and_annotation
{
    my ($self, $mapping, $spljnc_ori, $minus_splices, $plus_splices) = @_;
    my $ori = $self->_orient_by_tail($mapping);
    unless($ori) {  # if no obvious tail    
	$ori = $self->_orient_by_annot($mapping, $spljnc_ori);
	if($ori == 1 && $self->_find_internal_tail($mapping, -1) or $ori == 2 && $self->_find_internal_tail($mapping, 1)) {
	    $ori = 3;
	}
	elsif(!$spljnc_ori) {
	    if(($ori == 1 and $mapping->starts_near_T_stretch and !$mapping->has_plus_strand_polyA_signal) or
	       ($ori == 2 and $mapping->ends_near_A_stretch and !$mapping->has_minus_strand_polyA_signal)) {
		if($mapping->mRNA_acc_list) {
		    print STDERR "adjunct_T-rich: mRNA=",$mapping->mRNA_acc_string, " EST=", $mapping->EST_acc_string, "; ",
		    join(",",map{$_->mapping_id} $mapping->primary_mapping_list),"\n";
		}
		else {
		    $ori = 3;
		    $mapping->add_problem_tag(PROBLEM_TAG_METHOD, 'adjunct_T-rich');
		}
	    }
	    if($ori == 1) {
		$ori = 3 if ($self->_fuzzy_splice_overlap_check($mapping, $minus_splices));
	    }
	    elsif($ori == 2) {
		$ori = 3 if ($self->_fuzzy_splice_overlap_check($mapping, $plus_splices));
	    }
	}
    }
    if($ori and $spljnc_ori and $ori != 3 and $ori != $spljnc_ori) {
	$ori = 3;
	$mapping->add_problem_tag(PROBLEM_TAG_METHOD, "uncertain_splice_ori");
    }
    return $ori;
}


sub _update_evaluation
{
    my ($self, $mapping, $splice_ori, $test_ori) = @_;
    my $e = $self->evaluation;
    my $qType;
    if($mapping->refSeq_acc_list) { $qType = "RefSeq"; }
    elsif($mapping->mRNA_acc_list) { $qType = "mRNA"; }
    else { $qType = "EST" }

    my $key;
    if($splice_ori == 0 or $splice_ori == 3) {
	if($test_ori == 0 or $test_ori == 3) {
	    $key = "$qType splOri- testOri-";
	} 
	else {
	    $key = "$qType splOri- testOri+";
	}
    }
    elsif($test_ori == 0 or $test_ori == 3) {
	$key = "$qType splOri+ testOri-";
    }
    elsif($splice_ori == $test_ori) {
	$key = "$qType correct";
    }
    else {
	$key = "$qType incorrect";
    }

    $e->{$key}++;
}


# if mrna (refseq or not): trust start, end
# if est:
# 1. if strand is defined (1 or 2)
#  if strand == 1
#   trust start if fp_est
#   trust end if tp_est or tail
#  if strand == 2
#   trust start if tp_est or tail
#   trust end if fp_est
# 2. if strand is ambigous (3)
#  do not trust for now
# 3. if there is no data for the strand
#  do not trust (as there is no data for trust either)


sub _find_introns
{
    my ($self, $mappings, $seq) = @_;
    my %plus_splices;
    my %minus_splices;
    foreach my $mapping (@$mappings) {
	$self->_set_introns($mapping);
        my @hsps = $mapping->HSP_list;
        for my $i (0..@hsps-2) {
	    my $gap = $hsps[$i]->right_gap;
            if($gap->is_intron and $gap->jnc_type eq 'canonical') {
		my $start = $hsps[$i]->end;
		my $end = $hsps[$i+1]->start;
		if($gap->jnc_strand == 1) {
		    $plus_splices{$start}{$end}++;
		}
		else {
		    $minus_splices{$start}{$end}++;
		}
	    }
	}
    }
    
#    # debug
#    print STDERR "plus splices\n";
#    foreach my $start (keys %plus_splices) {
#	foreach my $end (keys %{$plus_splices{$start}}) {
#	    print STDERR " $start-$end\n";
#	}
#    }
#    print STDERR "\n";
#    print STDERR "minus splices\n";
#    foreach my $start (keys %minus_splices) {
#	foreach my $end (keys %{$minus_splices{$start}}) {
#	    print STDERR " $start-$end\n";
#	}
#    }
#    print STDERR "\n";

    return (\%plus_splices, \%minus_splices);
}


sub _set_introns
{
    my ($self, $mapping) = @_;
    my $min_intron_length = $self->{min_intron_length};
    my $nr_introns = 0;
    my $max_qInsert = 0;
    my @hsps = $mapping->HSP_list;
    for my $i (0..@hsps-2) {
	my $gap = $hsps[$i]->right_gap;
	my $jnc_type = $gap->jnc_type;
	my $qInsert = $gap->qInsert_length;
        if(!$gap->is_broken and
	   $qInsert == 0 and
	   $hsps[$i+1]->start - $hsps[$i]->end - 1 >= $min_intron_length and
	   ($jnc_type eq 'canonical' or $jnc_type eq 'minor')) {
	    $gap->is_intron(1);
	    $nr_introns++;
	}
	else {
	    $gap->is_intron(0);
	    $max_qInsert = $qInsert if($qInsert > $max_qInsert);
	}
    }
    $mapping->max_qInsert($max_qInsert);
    $mapping->nr_introns($nr_introns);
    print STDERR $mapping->acc_string, "\tnr_introns=", $nr_introns, "\n" if($self->debug);
}


sub _unset_reverse_introns
{
    my ($self, $mapping) = @_;
    my $strand = $mapping->strand;
    my $false = 0;
    my @hsps = $mapping->HSP_list;
    for my $i (0..@hsps-2) {
	my $gap = $hsps[$i]->right_gap;
	if($gap->is_intron and $gap->jnc_strand != $strand) {
	    $gap->is_intron(0);
	    $false++;
	}
    }
    $mapping->nr_introns($mapping->nr_introns - $false) if($false);
}


sub _orient_by_spljnc
{
    my ($self, $mapping, $seq) = @_;
    my (@plus, @minus, @ambiv);
    my $result;
    #my $confidence = 0;

    # Count splice junctions of different types
    my ($canon_plus, $canon_minus, $minor_plus, $minor_minus) = (0) x 4;
    for my $i (1..$mapping->nr_HSPs-1) {
	my $sep = $mapping->HSP($i)->right_gap;
	next unless($sep->is_intron);
	if($sep->jnc_type eq 'canonical') {
	    if($sep->jnc_strand == 1) { $canon_plus++;
				       }
	    else { $canon_minus++;
		   }
	}
	elsif($sep->jnc_type eq 'minor') {
	    if($sep->jnc_strand == 1) { $minor_plus++;
				        }
	    else { $minor_minus++;
		   }
	}
	#else {
	 #   $other++;
	#}
    }
    my $plus_count = $canon_plus + $minor_plus;
    my $minus_count = $canon_minus + $minor_minus;
    my $tot_count = $plus_count + $minus_count;
    
    # Evaluate the count
    if($tot_count == 0) {
	$result = 0;   # 0 = unspliced / no recognized junctions
    }
    elsif($canon_minus == 0 and $plus_count and
	  $plus_count >= 2*$minus_count) {
	$result = 1;   # 1 = plus strand
	#$confidence = 1 if($plus_count > 1);
	$mapping->add_tag_value('ori_support', "splice($plus_count)");
    }
    elsif($canon_plus == 0 and $minus_count and
	  $minus_count >= 2*$plus_count) {
	$result = 2;	# 2 = minus strand
	#$confidence = 1 if($minus_count > 1);
	$mapping->add_tag_value('ori_support', "splice($minus_count)");
    }
    else {
	$result = 3; 	# 3 = ambigous
	$mapping->add_problem_tag(PROBLEM_TAG_METHOD, 'splice_conflict');
    }

    print STDERR $mapping->qName_string,"\tspljnc\t$result\t$canon_plus+$minor_plus $canon_minus+$minor_minus\n" if($self->debug);
    return ($result);
    #return ($result, $confidence);
}


sub _fuzzy_splice_overlap_check
{
    my ($self,$mapping, $splices) = @_;
    #print STDERR "* FUZZY: checking ", $mapping->qName_string, "\n";
    my @hsps = $mapping->HSP_list;
    for my $i (0..@hsps-2) {
	my $gap = $hsps[$i]->right_gap;
	next if($gap->is_intron or $gap->is_broken);
	my $start = $hsps[$i]->end;
	my $end = $hsps[$i+1]->start;
	for my $s ($start-1..$start+1) {
	    my $some_splices = $splices->{$s} or next;
	    for my $e ($end-1..$end+1) {
		if($some_splices->{$e}) {
		    $mapping->add_problem_tag(PROBLEM_TAG_METHOD, 'opposite_splice');
		    return 1;  
		}
	    }
	}
    }
    return 0;
}


#sub _orient_by_spljnc__strict
#{
#    my ($self, $mapping, $seq) = @_;
#    my (@plus, @minus, @ambiv);
#    my $result;
#
#    # Count splice junctions of different types
#    my @ssCount = (0) x 7;
#    my @intron_coords = $mapping->tInsert_coords
#			(min_tInsert_length => MIN_INTRON_LEN_FOR_SS_CHECK);
#    foreach my $intron_pos (@intron_coords) {
#	my ($start, $end) = @$intron_pos;
#	#print STDERR "$start\t$end\n";
#	$start -= $seq->start - 1;
#	$end -= $seq->start - 1;
#	#print STDERR "$start\t$end\n";
#	my $ss1 = $seq->subseq($start, $start + 1);
#	my $ss2 = $seq->subseq($end - 1, $end);
#	#print STDERR "$ss1 $ss2\n";
#	#if($mapping->strand eq '-') {
#	    #$ss1, $ss2) = ($self->_revcom($ss2), $self->_revcom($ss1));
#	#}	
#	$ssCount[$self->_ss_type($ss1, $ss2)]++;
#    }
#    my $ok_count = $ssCount[1]+$ssCount[2]+$ssCount[3];
#    my $rev_count = $ssCount[4]+$ssCount[5]+$ssCount[6];
#    my $other_count = $ssCount[0];
#    my $tot_count = $ok_count + $rev_count + $other_count;
#
#    # Evaluate the count
#    if($tot_count == 0) {
#	$result = 0;   # 0 = unspliced
#    }
#    elsif($rev_count == 0 and $ok_count >= 2*$other_count) {
#	$result = ($mapping->strand_numeric == 1) ? 1 : 2;
#    }
#    elsif($ok_count == 0 and $rev_count >= 2*$other_count) {
#	$result = ($mapping->strand_numeric == 1) ? 2 : 1;
#    }
#    else {
#	$result = 3; 	# 3 = ambigous
#	$mapping->add_problem_tag(PROBLEM_TAG_METHOD, 'spljnc_conflict');
#    }
##    elsif($ok_count > $rev_count and $ssCount[1] >= 2*$ssCount[4]) {
##	$result = 2;   # 2 = plus
##    }
##    elsif($rev_count > $ok_count and $ssCount[4] >= 2*$ssCount[1]) {
##	$result = 3;	# 3 = minus
##    }
#
#    return $result;
#}
#

sub _orient_by_annot
{
    my ($self, $mapping, $splori) = @_;
    my $result = $self->_orient_by_annot_sub($mapping, $splori);
    print STDERR $mapping->qName_string,"\tannot\t$result\n" if($self->debug);
    $mapping->add_problem_tag(PROBLEM_TAG_METHOD, 'annot_conflict') if($result == 3);
    return $result;
}


sub _orient_by_annot_sub
{
    my ($self, $gf_mapping, $splori) = @_;
    my $result = 0;

    my @primary_mappings = $gf_mapping->primary_mapping_list;

    my $minus_polyA_signal = $gf_mapping->has_minus_strand_polyA_signal;
    my $plus_polyA_signal = $gf_mapping->has_plus_strand_polyA_signal;
    my $trusted_lib = $self->_mapping_is_from_trusted_lib($primary_mappings[0]);
    my $nr_EST_annot = 0;

    foreach my $raw_mapping (@primary_mappings) {
	my $strand = 0;
	if($raw_mapping->qType eq 'mRNA') {
	    $strand = ($raw_mapping->strand_numeric == 1) ? 1 : 2;
	    $gf_mapping->add_tag_value('ori_support', 'mRNA');
	    #print STDERR $raw_mapping->qName, " map_strand=",$raw_mapping->strand_numeric," ori=$strand \n";
	}
        elsif($raw_mapping->qType eq 'EST') {
	    my $readDir = $raw_mapping->qReadDir;
	    if($readDir == 5) {
		if($raw_mapping->strand_numeric == 1) {
		    $strand = 1 if($trusted_lib or $plus_polyA_signal or $splori == 1);    # could change $trusted_lib to ($trusted_lib and !$adjunct_T)
		}
		else {
		    $strand = 2 if($trusted_lib or $minus_polyA_signal or $splori == 2);
		}
	    }
	    elsif($readDir == 3) {
		if($raw_mapping->strand_numeric == 1) {
		    $strand = 2 if($trusted_lib or $minus_polyA_signal or $splori  == 2);
		}
		else {
		    $strand = 1 if($trusted_lib or $plus_polyA_signal or $splori == 1);
		}
	    }
	    elsif($readDir == -1 and $trusted_lib) {
		$strand = 3; #ambigous
	    }
	    else {
	        next; 	# no data
	    }
	    $nr_EST_annot++ if($strand == 1 or $strand == 2);
      	    #print STDERR $raw_mapping->qName, " readDir=$readDir map_strand=",$raw_mapping->strand_numeric," ori=$strand \n";
        }
        else {
	    warn ("Unrecognized qType for mapping");
	    next;
        }
        if($result == 0 or $result == 3) {
	    $result = $strand;
	}
	elsif($strand != 3 and $strand != 0 and $result != $strand) {
	    return 3;
	}
    }
    
    $gf_mapping->add_tag_value('ori_support', "EST-annot($nr_EST_annot)") if($nr_EST_annot);
    if($result == 1) {
	$gf_mapping->add_tag_value('ori_support', 'A-signal') if($plus_polyA_signal);
    }
    elsif($result == 2) {
	$gf_mapping->add_tag_value('ori_support', 'A-signal') if($minus_polyA_signal);
    }
    $gf_mapping->add_tag_value('ori_support', 'trusted_lib') if($trusted_lib);
    
    return $result;
}


sub _mapping_is_from_trusted_lib
{
    my ($self, $mapping) = @_;
    my $map_db = $mapping->db;
    my $id_hash;
    unless($id_hash = $self->{_trusted_lib_lists}{$map_db}) {
	my $ary = $map_db->get_libs_with_accurate_ori_annotation(
		    confidence_level => ANNORI_CONFIDENCE_LEVEL,
		    min_ppt_correct => ANNORI_MIN_PPT_CORRECT);
	$self->{_trusted_lib_lists}{$map_db} = $id_hash = {};
	print STDERR "got ids of ",scalar(@$ary)," trusted libs\n" if($self->debug);
	foreach my $lib_id (@$ary) {
	    $id_hash->{$lib_id} = 1;
	}
    }
    #print STDERR $mapping->qName, "\ttrust\t", ($id_hash->{$mapping->mRNAInfo->library} ? 1 : 0), "\n";
    return $id_hash->{$mapping->mRNAInfo->library} ? 1 : 0;
}


sub _set_EST_read_dir
{
    my ($self, $gf_mapping) = @_;
    my $result = 0;

    foreach my $raw_mapping ($gf_mapping->primary_mapping_list) {
	if($raw_mapping->qType eq 'EST') {
	    next unless($self->_mapping_is_from_trusted_lib($raw_mapping));
	    my $seq = $raw_mapping->query_seq;
	    if($seq) {
	        my $fp = $self->_find_5prime_text($seq->desc);
	        my $tp = $self->_find_3prime_text($seq->desc);
	        if($fp and $tp) {
		    $raw_mapping->qReadDir(-1);	#ambigous
		}
		elsif($fp) {
		    $raw_mapping->qReadDir(5);
		    $gf_mapping->has_5p_EST(1);
		}
		elsif($tp) {
		    $raw_mapping->qReadDir(3);
		    $gf_mapping->has_3p_EST(1);
		}
	    }
	    else {
	        warn ("No query seq for mapping");
	        next;
	    }
        }
	
    }
}


sub _find_5prime_text
{
    my($self, $str) = @_;
    $str =~ s/similar to.*//i;
    return ($str =~ /(?:5|five)(?:'|\-prime|\s+prime)/i) ? 1 : 0;
#    return ($str =~ /(?:5|five)(?:\-|\s*)(?:'|prime)/i) ? 1 : 0;
}


sub _find_3prime_text
{
    my($self, $str) = @_;
    $str =~ s/similar to.*//i;
    return ($str =~ /(?:3|three)(?:'|\-prime|\s+prime)/i) ? 1 : 0;
#    return ($str =~ /(?:3|three)(?:\-|\s*)(?:'|prime)/i) ? 1 : 0;
}


sub _orient_by_tail
{
    my ($self, $mapping, %args) = @_;

    # Gather evidence
    my ($plus, $minus) = (0,0);
    foreach my $m ($mapping->first_primary_mappings) {
	if($m->strand_numeric == 1) {
	    $minus = $self->_find_T_tail_for_mapping($m);
	}
	else {
	    $minus = $self->_find_A_tail_for_mapping($m);
	}
	last if($minus);
    }
    foreach my $m ($mapping->last_primary_mappings) {
	if($m->strand_numeric == 1) {
	    $plus = $self->_find_A_tail_for_mapping($m);
	}
	else {
	    $plus = $self->_find_T_tail_for_mapping($m);
	}
	last if($plus);
    }

    # Evaluate evidence
    my $result = 0;
    if($plus and $minus) {
	$result = 3;
	$mapping->add_problem_tag(PROBLEM_TAG_METHOD, 'tail_conflict');
    }
    elsif($plus) { $result = 1;  }
    elsif($minus) { $result = 2; }

    $mapping->add_tag_value('ori_support', 'A-tail') if($result == 1 or $result == 2);

    # Return result
    print STDERR $mapping->qName_string,"\ttail\t$result\n" if($self->debug);
    return $result;
}


sub _find_internal_tail
{
    my ($self, $mapping, $ori) = @_;
    my $min_tail_size = 10;	 # Min length of contigous A/T stretch
    my $max_tail_distance = 1;   # Max distance between tail and mapping

    my $got_tail = 0;

    if($ori != 1 and $ori != -1) { die "Invalid orientation $ori"; }

    foreach my $m ($mapping->last_primary_mappings) {
        if($m->strand_numeric == $ori) {
	    my $qSize = $m->qSize;
	    my $subseq_start = $m->qEnd+1;
	    next unless ($qSize - $subseq_start + 1 >= $min_tail_size);
	    my $subseq_end = $subseq_start + $min_tail_size + $max_tail_distance - 1;
	    $subseq_end = $qSize if($subseq_end > $qSize);
	    my $subseq = $m->query_seq->subseq($subseq_start, $subseq_end);
	    if($subseq =~ /A{$min_tail_size}/i) { $got_tail = 1; last; }
        }
	else {
	    my $subseq_end = $m->qStart - 1;
	    next unless ($subseq_end >= $min_tail_size);
	    my $subseq_start = $subseq_end - $min_tail_size - $max_tail_distance + 1;
	    $subseq_start = 1 if($subseq_start < 1);
	    my $subseq = $m->query_seq->subseq($subseq_start, $subseq_end);
	    if($subseq =~ /T{$min_tail_size}/i) { $got_tail = 1; last; }
	}
    }
#    elsif($ori == 2) {
#        foreach my $m ($mapping->last_primary_mappings) {
#	    if($m->strand_numeric == 1) {
#		my $subseq_end = $m->qStart-1;
#		next unless ($subseq_end >= $min_tail_size);
#		my $subseq_start = $subseq_end - $min_tail_size - $max_tail_distance + 1;
#		$subseq_start = 1 if($subseq_start < 1);
#		my $subseq = $m->query_seq->subseq($subseq_start, $subseq_end);
#		if($subseq =~ /T{$min_tail_size}/i) { $got_tail = 1; warn "1 1\n"; last; }
#            }
#	    else {
#		my $qSize = $m->qSize;
#		my $subseq_start = $m->qEnd+1;
#		next unless ($qSize - $subseq_start + 1 >= $min_tail_size);
#		my $subseq_end = $subseq_start + $min_tail_size + $max_tail_distance - 1;
#		$subseq_end = $qSize if($subseq_end > $qSize);
#		my $subseq = $m->query_seq->subseq($subseq_start, $subseq_end);
#		if($subseq =~ /A{$min_tail_size}/i) { $got_tail = 1; warn "1 2\n"; last; }
#	    }
#	}
#    }
#    else {
#	die "Invalid orientation $ori";
#    }

    $mapping->add_problem_tag(PROBLEM_TAG_METHOD, 'internal_tail') if($got_tail);

    return $got_tail;
}


sub _find_A_tail_for_mapping
{
    my ($self, $m) = @_;
    my $cutoff = 6;
    my $As = $m->qTermAs;
    my $unmapped = $m->qSize - $m->qEnd;
    return ($As >= $cutoff and
	    $unmapped >= $cutoff and
	    $unmapped - $As <= 1) ? 1 : 0;
}


sub _find_T_tail_for_mapping
{
    my ($self, $m) = @_;
    my $cutoff = 6;
    my $Ts = $m->qInitTs;
    my $unmapped = $m->qStart - 1;
    return ($Ts >= $cutoff and
	    $unmapped >= $cutoff and
	    $unmapped - $Ts <= 1) ? 1 : 0;
}


sub _find_A_tail_in_seq
{
    my ($self, $seq) = @_;
    my $seqstr = $seq->seq;
    return 2 if($seqstr =~ /A{13}$/i);  # last 13 bases are A's: definite A-tail
    my $subseqstr = substr($seqstr, length($seqstr)-20);
    return 1 if($subseqstr =~ /A{13}/i); # 13 contigous A's found withing the last 20: putative A-tail
    return 0;
}


sub _find_T_tail_in_seq
{
    my ($self, $seq) = @_;
    my $seqstr = $seq->seq;
    return 2 if($seqstr =~ /^T{13}/i);  # first 13 bases are T's: definite T-tail
    my $subseqstr = substr($seqstr, 0, 20);
    return 1 if($subseqstr =~ /T{13}/i); # 13 contigous T's found withing the first 20: putative T-tail
    return 0;
}


#sub _orient_by_tail
#{
#    my ($self, $mapping, %args) = @_;
#    my $result;
#
#    # Gather evidence
#    my ($plus, $minus) = (0,0);
#    foreach my $m ($mapping->first_primary_mappings) {
#	next unless($m->query_seq);
#	if($m->strand_numeric == 1) {
#	    $minus += AT::Tools::PolyATailFinder->find_T_tail_in_seq
#			    ($m->query_seq, %args);
#	}
#	else {
#	    $minus += AT::Tools::PolyATailFinder->find_A_tail_in_seq
#			    ($m->query_seq, %args);
#	}
#	last if($minus);
#    }
#    foreach my $m ($mapping->last_primary_mappings) {
#	next unless($m->query_seq);
#	if($m->strand_numeric == 1) {
#	    $plus += AT::Tools::PolyATailFinder->find_A_tail_in_seq
#			    ($m->query_seq, %args);
#	}
#	else {
#	    $plus += AT::Tools::PolyATailFinder->find_T_tail_in_seq
#			    ($m->query_seq, %args);
#	}
#	last if($plus);
#    }
#
#    # Evaluate evidence
#    if($plus and $minus) {
#	$result = 3;
#    }
#    elsif($plus) {
#	$result = 1;
#    }
#    elsif($minus) {
#	$result = 2;
#    }
#    else {
#	$result = 0;
#    }
#
#    # Return result
#    print STDERR $mapping->qName_string,"\ttail\t$result\n" if($self->debug);
#    return ($result);
#}


sub _set_stretch_adjacency
{
    my ($self, $mapping, $seq) = @_;
    $mapping->starts_near_T_stretch($self->_find_T_stretch($mapping, $seq));
    $mapping->ends_near_A_stretch($self->_find_A_stretch($mapping, $seq));
    print STDERR $mapping->qName_string,"\tT=",$mapping->starts_near_T_stretch," A=", $mapping->ends_near_A_stretch, "\n" if($self->debug);
}


sub _find_T_stretch
{
    my ($self, $mapping, $seq) = @_;

    # get genomic seq 12 bp into mapping and 14 bp outside
    my $up_start = $mapping->start - $seq->start + 1 - 14;
    if($up_start < 1) {
	warn "MappingCompass::_find_T_stretch: start of sequence ",$seq->id,":",$seq->start,"-",$seq->end," exceeded\n";
	return 0;
    }
    my $up_seq = $seq->subseq($up_start, $up_start+25);  #+24);

    # search for poly-T/poly-A stretch
    return AT::Tools::PolyNStretchFinder->find_stretch_in_seq
	($up_seq, 'T', win_size => 14, min_content => 10);
}


sub _find_A_stretch
{
    my ($self, $mapping, $seq) = @_;

    # get genomic seq 12 bp into mapping and 14 bp outside
    my $down_start = $mapping->end - $seq->start + 1 - 11; #9;
    if($down_start+25 > $seq->length) {
	warn "MappingCompass::_find_T_stretch: end of sequence ",$seq->id,":",$seq->start,"-",$seq->end," exceeded\n";
	return 0;
    }
    my $down_seq = $seq->subseq($down_start, $down_start+25); #+24);

    # search for poly-T/poly-A stretch
    return AT::Tools::PolyNStretchFinder->find_stretch_in_seq
	($down_seq, 'A', win_size => 14, min_content => 10);
}


sub _orient_by_stretch
{
    my ($self, $mapping, $seq, %args) = @_;
    my $result;

    $args{win_size} = 14 unless $args{win_size};
    $args{min_content} = 10 unless $args{min_content};

    # for both ends of mapping:
    # get genomic seq 10 bp into mapping and 15 bp outside
    my $up_start = $mapping->start - $seq->start + 1 - 15;
    my $down_start = $mapping->end - $seq->start + 1 - 13; #9;
    my $up_seq = $seq->subseq($up_start, $up_start+28);  #+24);
    my $down_seq = $seq->subseq($down_start, $down_start+28); #+24);

    #print "UP: $up_seq\n";
    #print "DOWN: $down_seq\n";

    # search for poly-T/poly-A stretch
    my $upT = AT::Tools::PolyNStretchFinder->find_stretch_in_seq($up_seq, 'T', %args);
    my $downA = AT::Tools::PolyNStretchFinder->find_stretch_in_seq($down_seq, 'A', %args);

    #print "$upT $downA\n";

    if($upT and $downA) {
	$result = 3;
    }
    elsif($downA) {
	$result = ($mapping->strand_numeric == 1) ? 1 : 2;
    }
    elsif($upT) {
	$result = ($mapping->strand_numeric == 1) ? 2 : 1;
    }
    else {
	$result = 0;
    }
    return $result;
}


sub _set_polyA_signals
{
    my ($self, $mapping, $seq) = @_;

    # for both ends of mapping:
    # get genomic seq 37 bp into mapping
    my $map_start = $mapping->start - $seq->start + 1;
    my $map_end = $mapping->end - $seq->start + 1;
    # return if sequence is too short
    if($map_start + 31 > $seq->length or $map_end - 31 < 1) {
	return;
    }
    my $up_seq = $seq->subseq($map_start + 10, $map_start + 31);
    my $down_seq = $seq->subseq($map_end - 31, $map_end - 10);

    my $fwd_re = 'A[AT]TAAA';
    my $rev_re = 'TTTA[AT]T';

    # look for forward in down region
    $mapping->has_plus_strand_polyA_signal( $down_seq =~ /$fwd_re/i ? 1 : 0 );
    # look for a reverse signal in up region
    $mapping->has_minus_strand_polyA_signal( $up_seq =~/$rev_re/i ? 1 : 0 );

    print STDERR $mapping->acc_string, "\tplus_signal=", $mapping->has_plus_strand_polyA_signal, " minus_signal=", $mapping->has_minus_strand_polyA_signal, "\n"
	if($self->debug);
}


sub _orient_by_polyA_signal
{
    my ($self, $mapping, $seq, %args) = @_;
    my $result;

    # for both ends of mapping:
    # get genomic seq 37 bp into mapping
    my $map_start = $mapping->start - $seq->start + 1;
    my $map_end = $mapping->end - $seq->start + 1;
    # return 0 (=not found) if sequence is too short
    if($map_start + 31 > $seq->length or $map_end - 31 < 1) {
	warn ("Sequence too short in polyA signal check for mapping of ".($mapping->qName_string)."\n");
	return 0;
    }
    my $up_seq = $seq->subseq($map_start + 10, $map_start + 31);
    my $down_seq = $seq->subseq($map_end - 31, $map_end - 10);

    ## strip A-stretch of up to 7 bp at the end
    #$up_seq =~ s/^T{0,7}//i;
    #$down_seq =~ s/A{0,7}$//i;
    ## restrict to the region to search
    #$up_seq = substr($up_seq, 9, 30);
    #$down_seq = substr($down_seq, -30, -9);
 
    my $fwd_re = $args{fwd_re} || 'A[AT]TAAA';
    my $rev_re = $args{rev_re} || 'TTTA[AT]T';

    # look for A(A/T)TAAA in down region
    my $downA = ($down_seq =~ /$fwd_re/i);
    #my $downA = ($down_seq =~ /A[ATG]TAAA|[TCG]ATAAA|AATA[TC]A/i);
    #my $downA = ($down_seq =~ /A[AT]TAAA/i);
    #my $downA = ($down_seq =~ /A[AT]TAAA.{0,29}A{0,7}$/i);
    # look for a revcom'd signal in up region
    my $upT = ($up_seq =~/$rev_re/i);
    #my $upT = ($up_seq =~/TTTA[TAC]T|TTTAT[AGC]|A[AG]TATT/i);
    #my $upT = ($up_seq =~/TTTA[AT]T/i);
    #my $upT = ($up_seq =~/^T{0,7}.{0,29}.TTTA[AT]T/i);

    if($upT and $downA) {
	$result = 3;
    }
    elsif($downA) {
	$result = 1;
    }
    elsif($upT) {
	$result = 2;
    }
    else {
	$result = 0;
    }

    print STDERR $mapping->qName_string,"\tsignal\t$result\n" if($self->debug);
    return $result;
    
}


sub __orient_by_polyA_signal
{
    my ($self, $mapping, $seq, %args) = @_;
    my $result;

    # for both ends of mapping:
    # get genomic seq 37 bp into mapping
    my $map_start = $mapping->tStart - $seq->start + 1;
    my $map_end = $mapping->tEnd - $seq->start + 1;
    #my $up_seq = $seq->subseq($map_start + 5, $map_start + 46);
    #my $down_seq = $seq->subseq($map_end - 46, $map_end - 5);
    my $up_seq = $seq->subseq($map_start + 9, $map_start + 30);
    my $down_seq = $seq->subseq($map_end - 30, $map_end - 9);

    my $fwd_re = $args{fwd_re} || 'A[AT]TAAA';
    my $rev_re = $args{rev_re} || 'TTTA[AT]T';

    # look for A(A/T)TAAA in down region
    my $downA = ($down_seq =~ /$fwd_re/i);
    #my $downA = ($down_seq =~ /A[ATG]TAAA|[TCG]ATAAA|AATA[TC]A/i);
    #my $downA = ($down_seq =~ /A[AT]TAAA/i);
    #my $downA = ($down_seq =~ /A[AT]TAAA.{0,29}A{0,7}$/i);
    # look for a revcom'd signal in up region
    my $upT = ($up_seq =~/$rev_re/i);
    #my $upT = ($up_seq =~/TTTA[TAC]T|TTTAT[AGC]|A[AG]TATT/i);
    #my $upT = ($up_seq =~/TTTA[AT]T/i);
    #my $upT = ($up_seq =~/^T{0,7}.{0,29}.TTTA[AT]T/i);

    if($upT and $downA) {
	$result = 3;
    }
    elsif($downA) {
	$result = ($mapping->strand_numeric == 1) ? 1 : 2;
    }
    elsif($upT) {
	$result = ($mapping->strand_numeric == 1) ? 2 : 1;
    }
    else {
	$result = 0;
    }
    return $result;
    
# Note: as a polyA tail left in the transcript may align with some A-residues
#   past the cleavage site in the genome, we allow for a few extra terminal A's/
#   initial T's. We allow a maximum of 7 extra terminal A's / initial T's
#   in the mapped region, as longer genomic polyA/polyT stretches are indicative
#   of internal priming, which would suggest that a putative polyA signal
#   is a false hit.
}


# 10 overrepr motifs: 
#/A[ATGC]TAAA/
#/[TCG]ATAAA/
#/AATA[TCG]A/

#/TTTA[TACG]T/
#/TTTAT[AGC]/
#/A[AGC]TATT/
1;
