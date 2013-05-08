# AT::FT::GFMapping module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::FT::GFMapping

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::FT::GFMapping;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use AT::FT::GFMappingHSP;
use AT::FT::GFMappingGap;
use AT::Tools::SeqHandler;

@ISA = qw/AT::Root/;


sub new  {
    my ($caller, %args) = @_;
    my $self = bless {
	gfmapping_id => ($args{gfmap_id} || undef),
	bin => ($args{bin} || undef),
	cluster_id => ($args{cluster_id} || 0),
	assembly => ($args{assembly} || undef),
	chr => ($args{chr} || undef),
	start => ($args{start} || undef),
	end => ($args{end} || undef),
	strand => ($args{strand} || 0),
	trust_start => ($args{trust_start} || 0),
	trust_end => ($args{trust_end} || 0),
	max_qInsert => ($args{max_qInsert} || 0),
        best_for_query => ($args{best_for_query} || 0),    # can be "T"rue, "F"alse, "A"mbiguous
	query_bp_past_right_border => ($args{query_bp_past_right_border} || 0),
	query_bp_past_left_border => ($args{query_bp_past_left_border} || 0),
	starts_near_T_stretch => ($args{starts_near_T_stretch} || undef),
	ends_near_A_stretch => ($args{ends_near_A_stretch} || undef),
	internally_primed => ($args{internally_primed} || 0),
	has_plus_strand_polyA_signal => ($args{has_plus_strand_polyA_signal} || 0),
	has_minus_strand_polyA_signal => ($args{has_minus_strand_polyA_signal} || 0),
	has_3p_EST => ($args{has_3p_EST} || 0),
	has_5p_EST => ($args{has_5p_EST} || 0),
	nr_introns => ($args{nr_introns} || 0),
	score1 => ($args{score1} || 0),
	score2 => ($args{score2} || 0),
	_tags => {},
	#_primary_mappingL => ($args{primary_mappings} || $args{mappings} || croak 'No primary mappings'),
	_EST_accL => ($args{EST_accs} || []),
	_mRNA_accL => ($args{mRNA_accs} || []),
	_HSPL => []
    }, ref $caller || $caller;
    $self->_init($args{primary_mappings} || $args{mappings});
    return $self;
}


sub _init {
    my ($self, $primary_mappings) = @_;
    if($primary_mappings) {
        # Get these attributes from the first mapping in our list of raw mappings.
        # This is ok since these attributes should be the same for all raw mappings
        # in the list.
	my $m = $primary_mappings->[0];
	my $mrna_info = $m->mRNAInfo;
        $self->{library} = $mrna_info->library unless(defined $self->{library});
        $self->{mrnaClone} = $mrna_info->mrnaClone unless(defined $self->{mrnaClone});
	$self->{assembly} =  $m->target_db unless($self->{assembly});
	$self->{chr} = $m->tName unless($self->{chr});
	# Set mapping and accession lists, keep mRNAs and ESTs separately
	my (@mrna_mappings, @est_mappings, @mrna_acc, @est_acc);
	foreach my $m (@$primary_mappings) {
	    if($m->qType eq 'mRNA') {
		push @mrna_mappings, $m;
		push @mrna_acc, $m->qName;
	    }
	    elsif($m->qType eq 'EST') {
		push @est_mappings, $m;
		push @est_acc, $m->qName;
	    }
	    else {
		warn "GFMapping: skipping primary mapping with unknown qType ".($m->qType);
	    }
	}
	$self->{_primary_mRNA_mappingL} = \@mrna_mappings;
	$self->{_primary_EST_mappingL} = \@est_mappings;
	$self->{_mRNA_accL} = \@mrna_acc;
	$self->{_EST_accL} = \@est_acc;
    }
}


sub HSP { $_[0]->{_HSPL}->[$_[1]-1]; }
sub nr_HSPs { scalar(@{$_[0]->{_HSPL}}); }

sub primary_mapping_list {
    my $self = shift;
    my @list = (@{$self->{_primary_mRNA_mappingL}}, @{$self->{_primary_EST_mappingL}});
    return @list;
}

sub primary_refSeq_mapping {
    my @m = grep {$_->qName =~ /^N[MR]_/} $_[0]->primary_mRNA_mapping_list;
    if(@m > 1) { warn "Several RefSeq mappings in GFMapping (",join(',',map{$_->qName}@m),")"; }
    return $m[0];
}

sub replace_primary_mapping {
    my ($self, $new_m) = @_;
    my $list;
    if($new_m->qType eq 'mRNA') { $list = $self->{_primary_mRNA_mappingL}; }
    elsif($new_m->qType eq 'EST') { $list = $self->{_primary_EST_mappingL}; }
    else{ warn "No qType for mapping of ".($new_m->qName); return; }
    for my $i (0..@$list-1) {
	if($list->[$i]->qName eq $new_m->qName) {
	    $list->[$i] = $new_m;
	    return 1;
       	}
    }
    warn "Could not replace mapping of ",$new_m->qName,"; not found\n";
    return 0;
}

sub qName_list {
    my $self = shift;
    my @list = (@{$self->{_mRNA_accL}}, @{$self->{_EST_accL}});
    return @list;
}

sub qName_string { return join(',',$_[0]->qName_list); }
sub mRNA_acc_string { return join(',',@{$_[0]->{_mRNA_accL}}); }
sub EST_acc_string { return join(',',@{$_[0]->{_EST_accL}}); }


sub acc_list { qName_list(@_) }
sub acc_string { qName_string(@_) }
sub refSeq_acc_list { grep {$_ =~ /^N[MR]_/} shift->mRNA_acc_list }

sub is_stretch_adjacent {
    my ($self) = shift;
    return ($self->{starts_near_T_stretch} or $self->{ends_near_A_stretch}) ? 1 : 0;
}

sub first_primary_mappings {
    my ($self, %args) = @_;
    my $glitch = $args{glitch} || 9;
    my @mappings = $self->primary_mapping_list;
    my $start = $mappings[0]->tStart;
    for my $i (1..@mappings-1) {
	$start = $mappings[$i]->tStart if($start > $mappings[$i]->tStart);
    }
    return grep {$_->tStart <= $start + $glitch} @mappings;
}

sub last_primary_mappings {
    my ($self, %args) = @_;
    my $glitch = $args{glitch} || 9;
    my @mappings = $self->primary_mapping_list;
    my $end = $mappings[0]->tEnd;
    for my $i (1..@mappings-1) {
	$end = $mappings[$i]->tEnd if($end < $mappings[$i]->tEnd);
    }
    return grep {$_->tEnd >= $end - $glitch} @mappings;
}


sub add_right_HSP {
    my ($self, $HSP, $sep) = @_;
    my $list = $self->{_HSPL};
    if(@$list) {
	croak "No gap" unless ($sep);
	$list->[-1]->next($HSP);
	$list->[-1]->right_gap($sep);
	$HSP->left_gap($sep);
	$self->{nr_unbroken_gaps}++ unless($sep->is_broken);
    }
    else {
	$self->start($HSP->start);
    }
    push @$list, $HSP;
    $self->end($HSP->end);
    return $HSP;
}


sub add_left_HSP {
    my ($self, $HSP, $sep) = @_;
    my $list = $self->{_HSPL};
    if(@$list) {
	croak "No gap" unless ($sep);
	$HSP->next($list->[0]);
	$HSP->right_gap($sep);
	$list->[0]->left_gap($sep);
	$self->{nr_unbroken_gaps}++ unless($sep->is_broken);
    }
    else {
	$self->end($HSP->end);
    }
    unshift @$list, $HSP;
    $self->start($HSP->start);
    return $HSP;
}


#sub trust_start { $_[0]->{_trust_start} = 1; }
#sub mistrust_start { $_[0]->{_trust_start} = 0; }
#sub trust_start_unless_mistrusted {
#    $_[0]->{_trust_start} = 1 unless defined($_[0]->{_trust_start}); }
#sub has_trusted_start { $_[0]->{_trust_start} || 0; }
#sub trust_end { $_[0]->{_trust_end} = 1; }
#sub mistrust_end { $_[0]->{_trust_end} = 0; }
#sub trust_end_unless_mistrusted {
#    $_[0]->{_trust_end} = 1 unless defined($_[0]->{_trust_end}); }
#sub has_trusted_end { $_[0]->{_trust_end} || 0; }


sub truncate_at_positions {
    my ($self, $start, $end) = @_;
    my $HSPs = $self->{_HSPL};

    #truncate start
    if($start ne $self->start) {
        while($HSPs->[0]->end < $start) {
	    shift @$HSPs;
	    croak "GFMapping: Error in HSP data structure" unless(@$HSPs);
        }
        $HSPs->[0]->start($start);
        $HSPs->[0]->left_gap(0);
        $self->start($start);
    }

    #truncate end
    if($end ne $self->end) {
        while($HSPs->[-1]->start > $end) {
	    pop @$HSPs;
	    croak "GFMapping: Error in HSP data structure" unless(@$HSPs);
        }
        $HSPs->[-1]->end($end);
        $HSPs->[-1]->right_gap(0);
        $self->end($end);
    }

    return $self;
}


sub qSeqStr_before_tPos {
    my ($self, $tPos) = @_;

    my @raw_mappings = $self->primary_mapping_list;
    croak "Method qSeqStr_before_tPos not implemented for multiple raw mappings "
	if(@raw_mappings > 1);

    my $m = $raw_mappings[0];
    my $qSeq = $m->query_seq;
    return undef unless($qSeq);

    my $qPos = $m->qPos_from_tPos_exact($tPos);    
    return undef unless ($qPos);
    return '' if($qPos == 1);   
    
    if($m->strand_numeric == 1) {
        return $qSeq->subseq(1, $qPos-1);
    }
    else {
	$qPos = $qSeq->length - $qPos + 1;
	return AT::Tools::SeqHandler->revcom($qSeq->subseq($qPos+1,
							   $qSeq->length));
    }
}


sub qSeqStr_after_tPos {
    my ($self, $tPos) = @_;

    my @raw_mappings = $self->primary_mapping_list;
    croak "Method qSeqStr_before_tPos not implemented for multiple raw mappings "
	if(@raw_mappings > 1);

    my $m = $raw_mappings[0];
    my $qSeq = $m->query_seq;
    return undef unless($qSeq);

    my $qPos = $m->qPos_from_tPos_exact($tPos);    
    return undef unless ($qPos);
    return '' if($qPos >= $qSeq->length);    # need >= here in case of version errors in the mapping db
   
    if($m->strand_numeric == 1) {
        return $qSeq->subseq($qPos+1, $qSeq->length);
    }
    else {
	$qPos = $qSeq->length - $qPos + 1;
	return AT::Tools::SeqHandler->revcom($qSeq->subseq(1, $qPos-1));
    }
}


sub merge {
    my ($mapping_a, $mapping_b, %args) = @_;
    my $self = $mapping_a;

    unless($mapping_a->_mRNA_span_check($mapping_b)) { return "Err 'mRNA span'"; }

    my $a = $mapping_a->{_HSPL};
    my $b = $mapping_b->{_HSPL};
    my $glitch = defined($args{'glitch'}) ? $args{'glitch'} : 5;

    my @HSPs;
    my @gaps;
    my($i, $j) = (0, 0);
 
    # Find the first element in a that does not end before b's first element
    while($i < @$a and $a->[$i]->end < $b->[0]->start) {
	push @HSPs, $a->[$i]->clone;
	push @gaps, $a->[$i]->left_gap->clone if($i);
	$i++;
    }  
    
    if($i == @$a) {
	# No overlap between $a and $b. Mark the break in the mapping
	push @gaps, AT::FT::GFMappingGap->new(is_broken => 1);
	push @HSPs, $b->[0]->clone;
	$j = 1;
    }
    else {
	# Some overlap between a and b
	# From this first overlap, all internal borders should agree
	my $prev_ender = $i ? $a->[$i-1] : undef;
	for(; $i < @$a and $j < @$b; $i++, $j++) {

	    my ($a_start, $a_end) = ($a->[$i]->start, $a->[$i]->end);
	    my ($b_start, $b_end) = ($b->[$j]->start, $b->[$j]->end);
	    my ($starter, $ender) = ($a->[$i], $a->[$i]); # default to $a

	    if($i and $a_start > $b_start + $glitch) {
		return "Err 1";
	    }
	    elsif($j and $a_start != $b_start) {
		if($b_start > $a_start + $glitch) {
		    return "Err 2";
		}
		# pick the most reliable of $a_start and $b_start
		$starter = $b->[$j]
		    if ($mapping_a->_gap_confidence($i,$i+1) <
			$mapping_b->_gap_confidence($j,$j+1));
	    }
    
	    if($i != @$a-1 and $a_end < $b_end - $glitch) {
		return "Err 3";
	    }
	    elsif($i == @$a-1 and $j == @$b-1) {
		$ender = $b->[$j] if ($b_end > $a_end);
	    }
	    elsif($j != @$b-1 and $a_end != $b_end) {
		if($b_end < $a_end - $glitch) {
		    return "Err 4";
		}
		# pick the most reliable of $a_end and $b_end
		$ender = $b->[$j]
		    if ($i == @$a-1 or
			$mapping_a->_gap_confidence($i+1,$i+2) <
			$mapping_b->_gap_confidence($j+1,$j+2));
	    }

	    push @HSPs, AT::FT::GFMappingHSP->new
		(start => $starter->start, end => $ender->end);
	    push @gaps, $prev_ender->right_gap->merge($starter->left_gap)
		if($prev_ender);
	    $prev_ender = $ender;
	}
    }

    # add HSPs after the region of overlap
    my ($c, $k) = ($i < @$a) ? ($a, $i) : ($b, $j);
    for(; $k < @$c; $k++) {
	push @HSPs, $c->[$k]->clone;
	push @gaps, $c->[$k]->left_gap->clone;
    }  

    # construct new mapping with the HSPs and gaps made
    my $end_mapping = ($mapping_a->end >= $mapping_b->end) ?
	$mapping_a : $mapping_b;
    my $merged = $self->new(mappings => [$mapping_a->primary_mapping_list,
					 $mapping_b->primary_mapping_list],
			    best_for_query => _merge_best_for_query($mapping_a, $mapping_b),
			    trust_start => $mapping_a->trust_start,
			    trust_end => $end_mapping->trust_end,
			    has_3p_EST => $mapping_a->has_3p_EST || $mapping_b->has_3p_EST,
			    has_5p_EST => $mapping_a->has_5p_EST || $mapping_b->has_5p_EST,	    
			    query_bp_past_left_border => $mapping_a->query_bp_past_left_border,
			    query_bp_past_right_border => $end_mapping->query_bp_past_right_border
			    #,aln_strand => $mapping_a->aln_strand
			    );
    $merged->add_right_HSP($HSPs[0]);
    for $i (1..@HSPs-1) {
	$merged->add_right_HSP($HSPs[$i], $gaps[$i-1]);
    }
    $merged->{_trust_start} = $mapping_a->{_trust_start};
    $merged->{_trust_end} = ($mapping_a->end > $mapping_b->end) ?
	$mapping_a->{_trust_end} : $mapping_b->{_trust_end};

    # return it
    return $merged;

# for if best for query
# if there is an mrna mapping: copy from that
# else
# if one spans all: copy from that
# else take worst: F, A, T

}


sub _mRNA_span_check
{
    my ($a, $b) = @_;
    if(($a->primary_mRNA_mapping_list and
	($b->start < $a->start or $b->end > $a->end))
        or
       ($b->primary_mRNA_mapping_list and
	($a->start < $b->start or $a->end > $b->end)))
    { return 0; }
    else
    { return 1; }
}


sub _merge_best_for_query
{
    my ($a, $b) = @_;
    if($a->primary_mRNA_mapping_list) {
      return $a->best_for_query;
    }
    elsif($b->primary_mRNA_mapping_list) {
      return $b->best_for_query;
    }
    elsif($a->start <= $b->start+5 and $a->end >= $b->end-5) {
      return $a->best_for_query;
    }
    elsif($b->start <= $a->start+5 and $b->end >= $a->end-5) {
      return $b->best_for_query;      
    }
    else {
      my $bfq_a = $a->best_for_query;
      my $bfq_b = $b->best_for_query;
      return "F" if($bfq_a eq 'F' or $bfq_b eq 'F');
      return "A" if($bfq_a eq 'A' or $bfq_b eq 'A');
      return "T";
    }
}


sub _gap_confidence
{
    my ($self, $i, $j) = @_;
    my $ceiling = 10;
    my $a = $self->HSP($i)->length;
    $a = $ceiling if ($a > $ceiling);
    my $b = $self->HSP($j)->length;
    $b = $ceiling if ($b > $ceiling);
    return $a * $b;
}


sub print_HSPs {
    my ($self, $out) = @_;
    $out = \*STDERR unless($out);
    print $out join(',',map {$_->qName} $self->primary_mapping_list),"\t";
    my $list = $self->{_HSPL};
    foreach my $HSP ($self->HSP_list) {
	print $out $HSP->start, '-', $HSP->end;
	if(my $sep = $HSP->right_gap) {
	    print $out ($sep->is_broken) ? '* *' : ' ';
	}
    }
    print $out "\n";
}


sub print_HSPs_and_gaps {
    my ($self, $out) = @_;
    $out = \*STDERR unless($out);
    print $out join(',',map {$_->qName} $self->primary_mapping_list),"\t";
    my $list = $self->{_HSPL};
    foreach my $HSP ($self->HSP_list) {
	print $out $HSP->start, '-', $HSP->end;
	if(my $sep = $HSP->right_gap) {
	    print $out ':',$sep->is_broken,'/',$sep->jnc_str,'/',$sep->jnc_type,':';
	}
    }
    print $out "\n";
}


sub is_spliced
{
    return shift->nr_introns ? 1 : 0;
}


sub loc_str
{
    my ($self) = @_;

    my $str = $self->chr.":".$self->start."-".$self->end.":".
	($self->strand == 1 ? '+' : '-');
    return $str;
}


sub overlaps
{
    my ($self, $other) = @_;
    return ($self->chr eq $other->chr and
	    $self->start <= $other->end and
	    $self->end >= $other->start)
	? 1 : 0;
}


# Tagging interface; should comform to the one in Bio::SeqFeature::Generic

sub has_tag  {

    my ($self, $tag) = @_;
    return exists $self->{'_tags'}->{$tag};
}

sub add_tag_value
{
    my $self = shift;
    my $tag = shift;
    $self->{'_tags'}->{$tag} ||= [];
    push(@{$self->{'_tags'}->{$tag}},@_);

}

sub get_tag_values  {

   my ($self, $tag) = @_;
   if( ! defined $tag ) { return (); }
   if ( ! exists $self->{'_tags'}->{$tag} ) {
       croak("asking for tag value that does not exist $tag");
   }
   return @{$self->{'_tags'}->{$tag}};

}

sub remove_tag  {

   my ($self, $tag) = @_;
   if ( ! exists $self->{'_tags'}->{$tag} ) {
	croak("trying to remove a tag that does not exist: $tag");
   }
   my @vals = @{$self->{'_tags'}->{$tag}};
   delete $self->{'_tags'}->{$tag};
   return @vals;
}

# Wrapper for tag named 'problem' for historical reasons

sub add_problem_tag
{
    my ($self, $method, $tag) = @_;
    $self->add_tag_value('problem', "${method}::${tag}");
}

sub problem_tag_list
{
    my ($self) = @_;
    return $self->has_tag('problem') ? $self->get_tag_values('problem') : ();
}

1;

