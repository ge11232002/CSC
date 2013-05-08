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


# Tests with Sim4 showed that sim4 could find marginal exons where there
# were 18 bp of EST sequence exactly matching the genomic. E.g.: AW468878.
# Sim4 did not find marginal exons if there was only 10 bp EST sequence
# corresponding to the exons (exactly matching). E.g.: C21231, AA884300.
# The program also had problems with longer sequences (up to 33 bp) of low
# quality. E.g.: AA356704.
# Conclusion: Sim4 is better than blat, but there is still a need for this module.


package AT::FT::Factory::MappingEndMover_old;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use AT::FT::GFMappingHSP;
use AT::FT::GFMappingGap;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Bio::AlignIO;

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
	min_HSP_len => ($args{min_HSP_len} || 10),
	word_size => ($args{word_size} || 4),
	range_size => ($args{range_size} || 20),
	match_score => ($args{match_score} || 1),
	mismatch_score => ($args{mismatch_score} || -2),
	aln_score_threshold => ($args{aln_score_threshold} || 5),
	extension_threshold => ($args{extension_threshold} || 9),
	gap_score => ($args{gap_score} || -3),
	verbose => ( $args{verbose} || 0 ),
	_word2start => undef,
	_word2end => undef,
	_starts => undef,
	_ends => undef,
	#_start2word => undef,
	#_end2word => undef,
    }, ref $caller || $caller;
    
    return $self;
}


sub run
{
    my ($self, %args) = @_;
    my $mappings = $args{mappings} || croak "No mappings arg";
    my $target_seq = $args{target_seq} || croak "No target_seq arg";    
 
    $self->_create_indices($mappings, $target_seq);
    $self->_print_indices() if($self->debug);

    my $nr_ends_moved = $self->_move_ends($mappings, $target_seq);
    $self->_verbose_msg("MappingEndMover: moved $nr_ends_moved ends\n");

    $self->_delete_indices();
}


sub _create_indices
{
    my ($self, $mappings, $target_seq) = @_;
    my $min_HSP_len = $self->min_HSP_len();
    my $word_size = $self->word_size();
    my (%start2word, %end2word);
    my (%word2starts, %word2ends);
    foreach my $m (@$mappings) {
	for my $i (1 .. $m->nr_HSPs-1) {
	    if($m->HSP($i)->length >= $min_HSP_len and
	       $m->HSP($i+1)->length >= $min_HSP_len and
	       $m->HSP($i)->right_gap->jnc_type eq 'canonical') {
		my $p1end = $m->HSP($i)->end;
		my $p2start = $m->HSP($i+1)->start;
		my $end_word;
		unless($end_word = $end2word{$p1end}) {
		    $end_word = $target_seq->subseq
			($p1end -$target_seq->start -$word_size +2,
			 $p1end -$target_seq->start +1);
		    $end_word = uc($end_word);
		    $end2word{$p1end} = $end_word;
		    push @{$word2ends{$end_word}}, $p1end;
			# we could extend this %end_idx to contain more info,
			# e.g. the mappings
		}
		# could use $end_word here to add more to the index
		my $start_word;
		unless($start_word = $start2word{$p2start}) {
		    $start_word = $target_seq->subseq
			($p2start -$target_seq->start +1,
			 $p2start -$target_seq->start +$word_size);
		    $start_word = uc($start_word);
		    $start2word{$p2start} = $start_word;
		    push @{$word2starts{$start_word}}, $p2start;
			# likewise, %start_idx could be extended
		}
	    }
	}
    }

    #$self->{_start2word} = \%start2word;   # need this?
    #$self->{_end2word} = \%end2word;       # need this?
    $self->{_starts} = [sort {$a <=> $b} keys(%start2word)];
    $self->{_ends} = [sort {$a <=> $b} keys(%end2word)];
    $self->{_word2starts} = \%word2starts;
    $self->{_word2ends} = \%word2ends;
}


sub _delete_indices
{
    my ($self) = @_;
    #$self->{_start2word} = undef;
    #$self->{_end2word} = undef;
    $self->{_starts} = undef;
    $self->{_ends} = undef;
    $self->{_word2starts} = undef;
    $self->{_word2ends} = undef;
}


sub _move_ends
{
    my ($self, $mappings, $target_seq) = @_;
    my $range_size = $self->range_size();
    my $nr_ends_moved = 0;

    foreach my $m (@$mappings) {

	# Check whether mapping ends at or just after a defined border
	my @r_borders = $self->_find_borders_within_range
	    ($self->{_ends}, $m->end -$range_size +1, $m->end);
	my ($r_hsps, $r_border);
	if(@r_borders) {
           #print STDERR "** several right borders! **\n" if(@r_borders > 1);
	    # Find the best way to break the mapping
            my @move_result = $self->_right_mover($m, $target_seq, @r_borders);
	    if(@move_result == 2) {
		($r_hsps, $r_border) = @move_result;
	    }
	    elsif($move_result[0]) {
		$m->query_bp_past_right_border($move_result[0]);
	    }
	    #print STDERR "MappingEndMover: right HSPS = ",
	    #    join(",", map { join ('-',@$_) } @$r_hsps), "\n--\n" if($r_hsps);
	}

	# Do the same for the start of the mapping
	my @l_borders = $self->_find_borders_within_range
	    ($self->{_starts}, $m->start, $m->start +$range_size -1);
	my ($l_hsps, $l_border);
	if(@l_borders) {
           #print STDERR "** several left borders! **\n" if(@l_borders > 1);
	    # Find the best way to break the mapping
	    my @move_result = $self->_left_mover($m, $target_seq, @l_borders);
	    if(@move_result == 2) {
		($l_hsps, $l_border) = @move_result;
	    }
	    elsif($move_result[0]) {
		$m->query_bp_past_left_border($move_result[0]);
	    }
	    #print STDERR "MappingEndMover: left HSPS = ",
	    #    join(",", map { join ('-',@$_) } @$l_hsps), "\n--\n" if($l_hsps);
	}

	# Default to original borders
	$r_border = $m->end unless($r_border);
	$l_border = $m->start unless($l_border);

	# Check that left and right movements do not conflict
	if($l_border >= $r_border) {
	    warn "MappingEndMover: left - right movement conflict for ".
		join(',', map { $_->qName } $m->primary_mapping_list)."\n";
	    next;
	}

	# Alter mapping
	#$m->print_HSPs(\*STDERR) if($l_hsps or $r_hsps);
	$m->truncate_at_positions($l_border, $r_border);
	if($l_hsps) {
	    my $HSP = AT::FT::GFMappingHSP->new
		(start => $l_hsps->[0]->[0],
		 end => $l_hsps->[-1]->[1]);
	    my $sep = AT::FT::GFMappingGap->new
		(jnc_str => $self->_get_jnc_str($target_seq,
			     		        $l_hsps->[-1]->[1]+1,
						$m->HSP(1)->start-1));
	    $m->add_left_HSP($HSP, $sep);
	}
	if($r_hsps) {
	    my $HSP = AT::FT::GFMappingHSP->new
		(start => $r_hsps->[0]->[0],
		 end => $r_hsps->[-1]->[1]);
	    my $sep = AT::FT::GFMappingGap->new
		(jnc_str => $self->_get_jnc_str($target_seq,
						$m->HSP($m->nr_HSPs)->end+1,
						$r_hsps->[0]->[0]-1));
	    $m->add_right_HSP($HSP, $sep);
	}
	#$m->print_HSPs(\*STDERR) if($l_hsps or $r_hsps);

	$nr_ends_moved++ if($l_hsps or $r_hsps);

    }

    return $nr_ends_moved;
}


sub _get_jnc_str
{
    my ($self, $target_seq, $jnc1start, $jnc2end) = @_;
    $jnc1start -= $target_seq->start - 1;
    $jnc2end -= $target_seq->start - 1;
    my $str = $target_seq->subseq($jnc1start, $jnc1start + 1).
	      $target_seq->subseq($jnc2end - 1, $jnc2end);
    return $str;   
}


# should implement this as a binary search rather than a linear
sub _find_borders_within_range
{
    my ($self, $borders, $start, $end) = @_;
    return (grep { $_ >= $start and $_ <= $end } @$borders);
}


sub _right_mover
# Try to move right end of mapping.
# If there was no movement giving a better score than the original placement:
#  return the number of bases of query seq after the first border.
{
    my ($self, $mapping, $t_seq, @borders) = @_;
    my $word_size = $self->word_size;

    # Find the first border and get query sequence from there
    my $first_border = $borders[0];
    for my $i (1..@borders-1) {
        $first_border = $borders[$i] if($borders[$i] < $first_border);
    }
    my @q_seq_ary = split(//, uc($mapping->qSeqStr_after_tPos($first_border)));

    # Construct array of seqs to align this query against
    # First, take the seq it's aligned against in the original mapping
    my $t_subseq_start = $first_border - $t_seq->start + 2;
    my $t_subseq_end = $mapping->end - $t_seq->start + 11; # add some bases in case blat was more stringent
    $t_subseq_end = $t_seq->length if($t_seq->length < $t_subseq_end);
    my $orig_t_subseq_str = uc($t_seq->subseq($t_subseq_start, $t_subseq_end));
    my @targets = ({d => 0,
		    border2 => 0,
		    seq_ary => [split(//,$orig_t_subseq_str)]}
		   );
    # Add unspliced target sequences
    foreach my $b (@borders) {
        my $q_seq_str_after_border = $mapping->qSeqStr_after_tPos($b) || next;
        #print STDERR "MappingEndMover: ",
         #   join(',',map {$_->qName} $mapping->primary_mapping_list),
         #   " end ",$mapping->end," -> border $b\t",($q_seq_str_after_border||'-'),"\n";
	next if (length($q_seq_str_after_border) < $word_size);
        my $word = substr($q_seq_str_after_border,0,$word_size);
        my $word_hits = $self->{_word2starts}->{$word};
        #print STDERR "MappingEndMover: hits = ",$word_hits?join(",",@$word_hits):'-',"\n";
        foreach my $word_hit (@$word_hits) {
            next if($word_hit <= $b);
	    # get subseq before border
	    my $t_subseq_str = substr($orig_t_subseq_str, 0, $b - $first_border);
	    # get subseq at new border
	    my $t_subseq_start = $word_hit - $t_seq->start + 1;
	    my $t_subseq_end = $t_subseq_start + @q_seq_ary + 14;
	    $t_subseq_end = $t_seq->length if($t_seq->length < $t_subseq_end);
	    $t_subseq_str .= uc($t_seq->subseq($t_subseq_start, $t_subseq_end));
	    # add to array of sequences
	    push @targets, { d => $b - $first_border,
			     border2 => $word_hit,
			     seq_ary => [split(//,$t_subseq_str)] };
        }
    }
    return scalar(@q_seq_ary) if (@targets == 1);

    # Align; find best target
    my ($hsps, $target) = $self->_align_and_select_best(\@q_seq_ary, \@targets);
    return scalar(@q_seq_ary) unless ($hsps);

    # Extract positions in target seq only
    # Offset these positions to absolute
    my $b2_abs_offset = $target->{border2} - $target->{d} - 1;
    my @target_hsps =
	map { [ $b2_abs_offset + $_->[2], $b2_abs_offset + $_->[3] ] }
	@$hsps;

    return (\@target_hsps, $first_border + $target->{d});
}


sub _left_mover
# Try to move left end of mapping.
# See _right_mover for further details.
{
    my ($self, $mapping, $t_seq, @borders) = @_;
    my $word_size = $self->word_size;

    # Find the first border and get query sequence from there
    my $first_border = $borders[0];
    for my $i (1..@borders-1) {
        $first_border = $borders[$i] if($borders[$i] > $first_border);
    }
    my @q_seq_ary = reverse
	split(//, uc($mapping->qSeqStr_before_tPos($first_border)));
    
    # Construct array of seqs to align this query against
    # First, take the seq it's aligned against in the original mapping
    my $t_subseq_start = $mapping->start - $t_seq->start - 9; # add some bases in case blat was more stringent
    $t_subseq_start = 1 if($t_subseq_start < 1);
    my $t_subseq_end = $first_border - $t_seq->start;
    my $orig_t_subseq_str = uc($t_seq->subseq($t_subseq_start, $t_subseq_end));
    my @targets = ({d => 0,         # distance from seq start to junction
		    border2 => 0,   # absolute position of new (distal) border
		    seq_ary => [reverse split(//, $orig_t_subseq_str)] }
		   );
    # Add unspliced target sequences
    foreach my $b (@borders) {
        my $q_seq_str_before_border = $mapping->qSeqStr_before_tPos($b) || '';
        #print STDERR "MappingEndMover: ",
         #   join(',',map {$_->qName} $mapping->primary_mapping_list),
         #   " start ",$mapping->start," -> border $b\t",($q_seq_str_before_border||'-'),"\n";
	next if(length($q_seq_str_before_border) < $word_size);
        my $word = substr($q_seq_str_before_border, -$word_size);
        my $word_hits = $self->{_word2ends}->{$word};
        #print STDERR "MappingEndMover: word $word -> hits ",$word_hits?join(",",@$word_hits):'-',"\n";
        foreach my $word_hit (@$word_hits) {
            next if($word_hit >= $b);
	    # get subseq before border
	    my $t_subseq_str = substr($orig_t_subseq_str,
				      -($first_border-$b),
				      $first_border-$b);
	    # get subseq at new border
	    my $t_subseq_end = $word_hit - $t_seq->start + 1;
	    my $t_subseq_start = $t_subseq_end - @q_seq_ary - 14;
	    $t_subseq_start = 1 if($t_subseq_start < 1);
	    $t_subseq_str = uc($t_seq->subseq($t_subseq_start, $t_subseq_end))
		.$t_subseq_str;
	    # add to array of sequences
	    push @targets, { d => $first_border - $b,
			     border2 => $word_hit,
			     seq_ary => [reverse split(//, $t_subseq_str)] };
        }
    }
    return scalar(@q_seq_ary) if (@targets == 1);

    # Align; find best target
    my ($hsps, $target) = $self->_align_and_select_best(\@q_seq_ary, \@targets);
    return scalar(@q_seq_ary) unless ($hsps);

    # Extract positions in target seq only
    # Offset these positions to absolute
    my $b2_abs_offset = $target->{border2} + $target->{d} + 1;
    my @target_hsps =
	reverse
	map { [ $b2_abs_offset - $_->[3], $b2_abs_offset - $_->[2] ] }
	@$hsps;

    return (\@target_hsps, $first_border - $target->{d});
}


sub _align_and_select_best
{
    my ($self, $q_seq_ary, $targets) = @_;

    # Align query against all target sequences and pick the best aln
    #  (this might seem to be a waste of cpu time since if we make > 2 alignments
    #   the upper-left HSPs of the DP matrices for some aligments may be
    #   identical; however it is very rare that we make > 2 alignments)
    my $best_score = $self->aln_score_threshold() - 1;
    my ($best_target, $best_aln);
    foreach my $t (@$targets) {
	my ($score, $aln) = $self->_dp_align($q_seq_ary, $t->{seq_ary});
	if($score == $best_score) { undef $best_target; }
	elsif($score > $best_score) { $best_score = $score;
				      $best_target = $t;
				      $best_aln = $aln; }
	#print STDERR "score = $score\n\n";
    }
    #print STDERR "best score = $best_score\n";

    # Were there several best alignments or none above the threshold?
    return undef unless($best_target);
    # Was the best aln that of the original mapping?
    return undef if($best_target eq $targets->[0]);

    # Get the aln after the splice junction as HSPs
    my $b1_rel_offset = $best_target->{d} + 1;
    my ($hsps, $subscore) = $self->_dpAlignment_to_hsps($best_aln, $b1_rel_offset);

    # Was the score after the splice junction above the threshold?
    return undef if($subscore < $self->aln_score_threshold);
    # Check that border2 is really covered by the alignment
    return undef if(@$hsps == 0 or $hsps->[0]->[2] != $b1_rel_offset);

    return ($hsps, $best_target);
}


# can be optimized as follows:
# first extend exact matches as far as possible
# then extend by dp until score drops below max_obtained - extension_threshold
# The exact extension step will simply fill in the diagonal of @N.@M
# Before starting the db we will fill in 1 col horiz and vert until below thr
# 
# Q: how to use a partial matrix?

sub _dp_align {
    my($self, $q, $t, $dir) = @_;

    #my @Q = ('-', split(//, $q)); # query sequence to be aligned with
    #my @T = ('-', split(//, $t)); # target sequence
    my @Q = ('-', @$q); # query sequence to be aligned with
    my @T = ('-', @$t); # target sequence
    print STDERR @$q, " vs\n",@$t,"\n";
    my (@M, @N); # score and traceback matrices
    my ($best_i, $best_j) = (0,0); # end position of best alignment

    # Get score parameters
    my $match_score = $self->match_score();
    my $mismatch_score = $self->mismatch_score();
    my $gap_score = $self->gap_score();
    my $extension_threshold = $self->extension_threshold();
    
    # Align (i.e. fill in score and traceback matrices)
    my $best_score = 0;
    foreach my $i (0..@Q-1) { $M[$i][0] = $gap_score * $i; }
    foreach my $j (0..@T-1) { $M[0][$j] = 0; }
    foreach my $i (1..@Q-1) {
	foreach my $j (1..@T-1) {
	    my $p =  ($Q[$i] eq $T[$j]) ? $match_score : $mismatch_score;
	    my $s1 = $M[$i-1][$j-1] + $p;
	    my $s2 = $M[$i-1][$j] + $gap_score;
	    my $s3 = $M[$i][$j-1] + $gap_score;
	    my ($k, $s);
	    if($s1 > $s2) {
		if($s1 > $s3) { $k = 0; $s = $s1; }
		else { $k = 2; $s = $s3; }
	    }
	    else {
		if($s2 > $s3) { $k = 1; $s = $s2; }
		else { $k = 2; $s = $s3; }
	    }
	    $N[$i][$j] = $k;
	    $M[$i][$j] = $s;
	    if($s > $best_score) {
		$best_score = $s;
		($best_i, $best_j) = ($i, $j); 
	    }
	}
    }

    # Output alignment (for DEBUG)
    $self->_debug_print_aln(\@Q, \@T, \@N, $best_i, $best_j);

    return ($M[$best_i][$best_j], [\@M, \@N, $best_i, $best_j]);
}


sub _align {
    my($self, $q, $t, $dir) = @_;

    my @Q = ('-', @$q); # query sequence to be aligned with
    my @T = ('-', @$t); # target sequence
    #print STDERR @Q, " vs\n",@T,"\n";
    my (@M, @N); # score and traceback matrices
    my $best_score = 0;

    # Get score parameters
    my $match_score = $self->match_score();
    my $mismatch_score = $self->mismatch_score();
    my $gap_score = $self->gap_score();
    my $extension_threshold = $self->extension_threshold();
    my $minus_infinity = -1000;

    # Extend as far as there are exact matches
    $M[0][0] = 0; $N[0][0] = 0;
    my $I = 1;
    while($I < @Q and $I < @T and $Q[$I] eq $T[$I]) {
	$best_score += $match_score;
	$M[$I][$I] = $best_score;
	$N[$I][$I] = 1;
	$I++;
    }
    my ($best_i, $best_j) = ($I-1,$I-1); # end position of best alignment so far
    
    # Compute scores for exact alignment followed by only gaps
    my ($i, $j) = ($I-1, $I);
    my $min_score = $best_score - $extension_threshold;
    for(my $s = $best_score-$gap_score; $s >= $min_score; $s -= $gap_score) {
	$M[$j][$i] = $M[$i][$j] = $s;
	$N[$j][$i] = 1;
	$N[$i][$j] = 2;
	$j++;
    }
    $M[$j][$i] = $M[$i][$j] = $minus_infinity;
       
    # DP align
    my ($next_start_col, $next_end_col) = ($I, @T-1);
    for my $i ($I..@Q-1) {
	my ($start_col, $end_col) = ($next_start_col, $next_end_col);
	$next_start_col = undef;
	for my $j ($start_col..$end_col) {
	    my ($k, $s);
	    my $p =  ($Q[$i] eq $T[$j]) ? $match_score : $mismatch_score;
	    my $s1 = $M[$i-1][$j-1] + $p;
	    my $s2 = $M[$i-1][$j] + $gap_score;
	    my $s3 = $M[$i][$j-1] + $gap_score;
	    if($s1 > $s2) {
		if($s1 > $s3) { $k = 0; $s = $s1; }
		else { $k = 2; $s = $s3; }
	    }
	    else {
		if($s2 > $s3) { $k = 1; $s = $s2; }
		else { $k = 2; $s = $s3; }
	    }
	    if($s < $min_score) {
		$M[$i][$j] = $minus_infinity;
	    }
	    else {
		$next_start_col = $j unless(defined($next_start_col));
		$next_end_col = $j;
		$N[$i][$j] = $k;
		$M[$i][$j] = $s;
		if($s > $best_score) {
		    $best_score = $s;
		    ($best_i, $best_j) = ($i, $j); 
		}
	    }
	}
	last unless($next_start_col);
	$next_end_col++;
    }

    # Output alignment (for DEBUG)
    #$self->_debug_print_aln(\@Q, \@T, \@N, $best_i, $best_j);

    return ($M[$best_i][$best_j], [\@M, \@N, $best_i, $best_j]);
}




sub _dpAlignment_to_hsps
{
    my ($self, $aln, $tStart) = @_;

    # Reconstruct alignment as a set of HSPs
    my ($M, $N, $last_i, $last_j) = @$aln;
    my ($i, $j) = ($last_i, $last_j);
    my ($hsp_qEnd, $hsp_tEnd);
    my @hsps;
    while($i >= 1 and $j >= $tStart) {
	if($N->[$i][$j] == 0) {		# 0 = match
	    ($hsp_qEnd, $hsp_tEnd) = ($i, $j) unless($hsp_qEnd);
	    if($i == 1 or $j == $tStart or $N->[$i-1][$j-1]) {
		push @hsps, [$i, $hsp_qEnd, $j, $hsp_tEnd];
		$hsp_qEnd = 0;
	    }
	    $i--; $j--;
	}
	elsif($N->[$i][$j] == 1) {	# 1 = tGap
	    $i--;
        }
        else {   			# 2 = qGap
	    $j--;
        }
    }
    @hsps = reverse(@hsps);
    my $subscore = $M->[$last_i][$last_j] - $M->[$i][$j];
    print STDERR "Qry HSPs: ", join(", ", map {$_->[0].'-'.$_->[1]} @hsps), "\n";
    print STDERR "Tgt HSPs: ", join(", ", map {$_->[2].'-'.$_->[3]} @hsps), "\n";
    print STDERR "Subscore = $subscore; Totscore = ",$M->[$last_i][$last_j],"\n";

    # return hsps
    return (\@hsps, $subscore);
}


sub _debug_print_aln
{
    my ($self, $Q, $T, $N, $last_i, $last_j) = @_;

    if($last_i == 0 or $last_j == 0) {
	print STDERR "(no alignment)\n";
	return;
    }

    # Reconstruct alignment    
    my @qAli;# = ($Q->[$i]);
    my @tAli;# = ($T->[$j]);
    my ($i, $j) = ($last_i, $last_j);
    while($i > 0 and $j > 0) {
	if($N->[$i][$j] == 0) {		# 0 = match
	    push @qAli, $Q->[$i--];
	    push @tAli, $T->[$j--];
	}
	elsif($N->[$i][$j] == 1) {	# 1 = tGap
	    push @qAli, $Q->[$i--];
	    push @tAli, '-';
	}
	else {   			# 2 = qGap
	    push @qAli, '-';
	    push @tAli, $T->[$j--];
        }
    }
    
    # Output to STDERR in ClustalW format
    my $qSeq = Bio::LocatableSeq->new(-id => 'query',
				      -seq => join('',reverse @qAli),
				      -start => 1,
				      -end => $last_i);
    my $tSeq = Bio::LocatableSeq->new(-id => 'target',
				      -seq => join('',reverse @tAli),
				      -start => 1,
				      -end => $last_j);
    my $aln = Bio::SimpleAlign->new();
    $aln->add_seq($qSeq); $aln->add_seq($tSeq);
    my $aln_out = Bio::AlignIO->new(-fh => \*STDERR, -format => 'clustalw');
    $aln_out->write_aln($aln);
#    print STDERR "Query:\t", reverse(@qAli), "\nTarget:\t", reverse(@tAli), "\n";
} 


sub _max_element
{
    my ($self, @A) = @_;
    my $max = 0;
    for my $i (1..@A-1) {
	$max = $i if($A[$max] < $A[$i]);
    }
    return $max;
}


sub _print_indices
{
    my ($self) = @_;
    my $idx;
    print STDERR "---\nword2starts index\n";
    $idx = $self->{_word2starts};
    foreach my $word (sort {$a cmp $b} keys %$idx) {
	print STDERR $word, "\t", join(",",@{$idx->{$word}}), "\n";
    }
    print STDERR "---\nword2ends index\n";
    $idx = $self->{_word2ends};
    foreach my $word (sort {$a cmp $b} keys %$idx) {
	print STDERR $word, "\t", join(",",@{$idx->{$word}}), "\n";
    }
#    print STDERR "---\nstart2word index\n";
#    $idx = $self->{_start2word};
#    foreach my $coord (sort {$a <=> $b} keys %$idx) {
#	print STDERR $coord, "\t", $idx->{$coord}, "\n";
#    }
#    print STDERR "---\nend2word index\n";
#    $idx = $self->{_end2word};
#    foreach my $coord (sort {$a <=> $b} keys %$idx) {
#	print STDERR $coord, "\t", $idx->{$coord}, "\n";
#    }
}


1;
