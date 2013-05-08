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


package AT::FT::Factory::MappingEndMover;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use AT::FT::GFMappingHSP;
use AT::FT::GFMappingGap;
use AT::MappingFactory;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Bio::AlignIO;

use constant FWD_IDX => scalar 0;
use constant REV_IDX => scalar 1;

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
	max_search_intron_size => ($args{max_intron_size} || 120000),
	max_intron_size => ($args{max_intron_size} || 7500),
	word_size => ($args{word_size} || 4),
	range_size => ($args{range_size} || 20),
	match_score => ($args{match_score} || 1),
	mismatch_score => ($args{mismatch_score} || -2),
	aln_score_threshold => ($args{aln_score_threshold} || 10),
	extension_threshold => ($args{extension_threshold} || 9),
	gap_score => ($args{gap_score} || -3),
	verbose => ( $args{verbose} || 0 ),
	move_mRNA_mappings_only => ( $args{move_mRNA_mappings_only} || 0),
	replace_primary_mappings => ($args{replace_primary_mappings} || 0),
	primap_outstream => $args{primap_outstream},
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
    my @start2word_ary = ({}, {});
    my @end2word_ary = ({}, {});
    my @word2starts_ary = ({}, {});
    my @word2ends_ary = ({}, {});
    foreach my $m (@$mappings) {
	#print STDERR "strand_idx = $strand_idx\t", $m->strand, "\n";
	#my $start2word = $start2word_ary[$strand_idx];
	#my $end2word = $end2word_ary[$strand_idx];
	#my $word2starts = $word2starts_ary[$strand_idx];
	#my $word2ends = $word2ends_ary[$strand_idx];
	for my $i (1 .. $m->nr_HSPs-1) {
	    my $gap = $m->HSP($i)->right_gap;
	    if($m->HSP($i)->length >= $min_HSP_len and
	       $m->HSP($i+1)->length >= $min_HSP_len and
	       $gap->jnc_type eq 'canonical') {
	        my $strand_idx = $gap->jnc_strand == 1 ? FWD_IDX : REV_IDX;
		my $p1end = $m->HSP($i)->end;
		my $p2start = $m->HSP($i+1)->start;
		my $end_word;
		unless($end_word = $end2word_ary[$strand_idx]{$p1end}) {
		    $end_word = $target_seq->subseq
			($p1end -$target_seq->start -$word_size +2,
			 $p1end -$target_seq->start +1);
		    $end_word = uc($end_word);
		    $end2word_ary[$strand_idx]{$p1end} = $end_word;
		    push @{$word2ends_ary[$strand_idx]{$end_word}}, $p1end;
			# we could extend this %end_idx to contain more info,
			# e.g. the mappings
		}
		# could use $end_word here to add more to the index
		my $start_word;
		unless($start_word = $start2word_ary[$strand_idx]{$p2start}) {
		    $start_word = $target_seq->subseq
			($p2start -$target_seq->start +1,
			 $p2start -$target_seq->start +$word_size);
		    $start_word = uc($start_word);
		    $start2word_ary[$strand_idx]{$p2start} = $start_word;
		    push @{$word2starts_ary[$strand_idx]{$start_word}}, $p2start;
			# likewise, %start_idx could be extended
		}
	    }
	}
    }

    #$self->{_start2word} = \%start2word;   # need this?
    #$self->{_end2word} = \%end2word;       # need this?
    $self->{_starts} = [ map { [sort {$a <=> $b} keys(%$_) ] } @start2word_ary ];
    $self->{_ends} = [ map { [ sort {$a <=> $b} keys(%$_) ] } @end2word_ary ];
    $self->{_word2starts} = \@word2starts_ary;
    $self->{_word2ends} = \@word2ends_ary;
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
    my $primap_outstream = $self->primap_outstream;
    my $replace_primary_mappings = $self->replace_primary_mappings;

    foreach my $m (@$mappings) {

	next if($self->{move_mRNA_mappings_only} and $m->mRNA_acc_list == 0);

	# Check whether mapping ends at or just after a defined border
	my @r_borders = $self->_find_borders_within_range
	    ($self->{_ends}, $m->end -$range_size +1, $m->end);
	my ($r_hsps, $r_border, $r_score);
	if(@{$r_borders[FWD_IDX]} or @{$r_borders[REV_IDX]}) {
           #print STDERR "** several right borders! **\n" if(@r_borders > 1);
	    # Find the best way to break the mapping
            my @move_result = $self->_right_mover($m, $target_seq, @r_borders);
	    if(@move_result == 3) {
		($r_hsps, $r_border, $r_score) = @move_result;
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
	my ($l_hsps, $l_border, $l_score);
	if(@{$l_borders[FWD_IDX]} or @{$l_borders[REV_IDX]}) {
           #print STDERR "** several left borders! **\n" if(@l_borders > 1);
	    # Find the best way to break the mapping
	    my @move_result = $self->_left_mover($m, $target_seq, @l_borders);
	    if(@move_result == 3) {
		($l_hsps, $l_border, $l_score) = @move_result;
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
		(start => $l_hsps->[0]->[2],
		 end => $l_hsps->[-1]->[3]);
	    my $sep = AT::FT::GFMappingGap->new
		(jnc_str => $self->_get_jnc_str($target_seq,
			     		        $l_hsps->[-1]->[3]+1,
						$m->HSP(1)->start-1));
	    $m->add_left_HSP($HSP, $sep);
	    print STDERR "MappingEndMover: moved left end of ",$m->qName_string," from ", $l_border, " to ", $l_hsps->[0][2], " (", $HSP->length, "nt)\n" if($self->debug);
	}
	if($r_hsps) {
	    my $HSP = AT::FT::GFMappingHSP->new
		(start => $r_hsps->[0]->[2],
		 end => $r_hsps->[-1]->[3]);
	    my $sep = AT::FT::GFMappingGap->new
		(jnc_str => $self->_get_jnc_str($target_seq,
						$m->HSP($m->nr_HSPs)->end+1,
						$r_hsps->[0]->[2]-1));
	    $m->add_right_HSP($HSP, $sep);
	    print STDERR "MappingEndMover: moved right end of ",$m->qName_string," from ", $r_border, " to ", $r_hsps->[-1][3], " (", $HSP->length, "nt)\n" if($self->debug);
	}
	#$m->print_HSPs(\*STDERR) if($l_hsps or $r_hsps);

	if($l_hsps or $r_hsps) {
	    $nr_ends_moved++;
	    if($primap_outstream or $replace_primary_mappings) {
		my $new_primap = $self->_create_new_primap($m->primary_mapping_list, $l_hsps, $r_hsps, $target_seq);
		_hack_in_score($new_primap, $m, $l_score, $r_score);
		$primap_outstream->write_mapping($new_primap) if($primap_outstream);
		$m->replace_primary_mapping($new_primap) if ($replace_primary_mappings);
	    }
	}

    }

    return $nr_ends_moved;
}


sub _hack_in_score
{
    my ($m, $gfm, $l_score, $r_score) = @_;
    my $string = "";
    if($l_score) {
      my $d = $gfm->HSP(2)->start - $gfm->HSP(1)->end - 1;
      $string .= "L_${l_score}_$d";
    }
    if($r_score) {
      my $nr_hsps = $gfm->nr_HSPs;
      my $d = $gfm->HSP($nr_hsps)->start - $gfm->HSP($nr_hsps-1)->end - 1;
      $string .= "R_${r_score}_$d";
    }
    $m->bin($string);
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
    return (map { [ grep { $_ >= $start and $_ <= $end } @$_ ] } @$borders);
}


sub _right_mover
# Try to move right end of mapping.
# If there was no movement giving a better score than the original placement:
#  return the number of bases of query seq after the first border.
{
    my ($self, $mapping, $t_seq, @borders_by_strand) = @_;
    my $word_size = $self->word_size;
    my $max_intron_size = $self->max_search_intron_size;

    # Find the first border and get query sequence from there
    my $first_border;
    foreach my $border_list (@borders_by_strand) {
	foreach my $border (@$border_list) {
	    $first_border = $border if(!defined($first_border) or $border < $first_border);
	}
    }
    my @q_seq_ary = split(//, uc($mapping->qSeqStr_after_tPos($first_border)));

    # Construct array of seqs to align this query against
    # First, take the seq it's aligned against in the original mapping
    my $t_subseq_start = $first_border - $t_seq->start + 2;
    my $t_subseq_end = $mapping->end - $t_seq->start + 11; # add some bases in case blat was more stringent
    $t_subseq_end = $t_seq->length if($t_seq->length < $t_subseq_end);
    my $orig_t_subseq_str = uc($t_seq->subseq($t_subseq_start, $t_subseq_end));
    my @targets = ({d => 0,
		    #border2 => 0,
		    intron_size => 0,
		    seq_ary => [split(//,$orig_t_subseq_str)]}
		   );
    # Add unspliced target sequences
    for my $strand_idx (FWD_IDX, REV_IDX) {
	foreach my $b (@{$borders_by_strand[$strand_idx]}) {
            my $q_seq_str_after_border = $mapping->qSeqStr_after_tPos($b) || next;
            #print STDERR "MappingEndMover: ",
             #   join(',',map {$_->qName} $mapping->primary_mapping_list),
             #   " end ",$mapping->end," -> border $b\t",($q_seq_str_after_border||'-'),"\n";
            next if (length($q_seq_str_after_border) < $word_size);
            my $word = substr($q_seq_str_after_border,0,$word_size);
            my $word_hits = $self->{_word2starts}[$strand_idx]{$word};
            #print STDERR "MappingEndMover: hits = ",$word_hits?join(",",@$word_hits):'-',"\n";
	    foreach my $word_hit (@$word_hits) {
	        next if($word_hit <= $b or $word_hit > $b+$max_intron_size);
	        # get subseq before border
	        my $t_subseq_str = substr($orig_t_subseq_str, 0, $b - $first_border);
	        # get subseq at new border
	        my $t_subseq_start = $word_hit - $t_seq->start + 1;
	        my $t_subseq_end = $t_subseq_start + @q_seq_ary + 14;
		$t_subseq_end = $t_seq->length if($t_seq->length < $t_subseq_end);
		$t_subseq_str .= uc($t_seq->subseq($t_subseq_start, $t_subseq_end));
		# add to array of sequences
		push @targets, { d => $b - $first_border,
				 #border2 => $word_hit,
				 intron_size => $word_hit - $b - 1,
				 seq_ary => [split(//,$t_subseq_str)] };
	    }
        }
    }
    return scalar(@q_seq_ary) if (@targets == 1);

    # Align; find best target
    my ($hsps, $target, $score) = $self->_align_and_select_best(\@q_seq_ary, \@targets);
    return scalar(@q_seq_ary) unless ($hsps);

    # Offset hsp positions to absolute
    #my $b2_abs_offset = $target->{border2} - $target->{d} - 1;
    my $b2_abs_offset = $first_border + $target->{intron_size};
    my $primap = ($mapping->primary_mapping_list)[0];
    my $q_offset = $primap->qPos_from_tPos_exact($first_border);
    my @hsps =
	    map { [ $q_offset + $_->[0], $q_offset + $_->[1],
	        $   b2_abs_offset + $_->[2], $b2_abs_offset + $_->[3] ] }
	    @$hsps;
    return (\@hsps, $first_border + $target->{d}, $score);
}


sub _left_mover
# Try to move left end of mapping.
# See _right_mover for further details.
{
    my ($self, $mapping, $t_seq, @borders_by_strand) = @_;
    my $word_size = $self->word_size;
    my $max_intron_size = $self->max_search_intron_size;

    # Find the first border and get query sequence from there
    my $first_border = 0;
    foreach my $border_list (@borders_by_strand) {
	for my $border (@$border_list) {
	    $first_border = $border if($border > $first_border);
	}
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
		    #border2 => 0,   # absolute position of new (distal) border
		    intron_size => 0,
		    seq_ary => [reverse split(//, $orig_t_subseq_str)] }
		   );
    # Add unspliced target sequences
    for my $strand_idx (FWD_IDX, REV_IDX) {
	foreach my $b (@{$borders_by_strand[$strand_idx]}) {
            my $q_seq_str_before_border = $mapping->qSeqStr_before_tPos($b) || '';
            #print STDERR "MappingEndMover: ",
             #   join(',',map {$_->qName} $mapping->primary_mapping_list),
             #   " start ",$mapping->start," -> border $b\t",($q_seq_str_before_border||'-'),"\n";
	    next if(length($q_seq_str_before_border) < $word_size);
	    my $word = substr($q_seq_str_before_border, -$word_size);
	    my $word_hits = $self->{_word2ends}[$strand_idx]{$word};
	    #print STDERR "MappingEndMover: word $word -> hits ",$word_hits?join(",",@$word_hits):'-',"\n";
	    foreach my $word_hit (@$word_hits) {
                next if($word_hit >= $b or $word_hit < $b-$max_intron_size);
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
			     #border2 => $word_hit,
				 intron_size => $b - $word_hit - 1,
			     seq_ary => [reverse split(//, $t_subseq_str)] };
	    }
        }
    }
    return scalar(@q_seq_ary) if (@targets == 1);

    # Align; find best target
    my ($hsps, $target, $score) = $self->_align_and_select_best(\@q_seq_ary, \@targets);
    return scalar(@q_seq_ary) unless ($hsps);

    # Offset hsp positions to absolute
    #my $b2_abs_offset = $target->{border2} + $target->{d} + 1;
    my $b2_abs_offset = $first_border - $target->{intron_size};
    my $primap = ($mapping->primary_mapping_list)[0];
    my $q_offset = $primap->qPos_from_tPos_exact($first_border) - 1;
    my @hsps =
	    reverse
	    map { [ $q_offset - $_->[1] + 1, $q_offset - $_->[0] + 1,
	            $b2_abs_offset - $_->[3], $b2_abs_offset - $_->[2] ] }
	    @$hsps;

    return (\@hsps, $first_border - $target->{d}, $score);
}

# if - strand:
# 
#
# query coords are in reversed seq; reverse if + strand
# 

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
	my ($score, $aln) = $self->_align($q_seq_ary, $t->{seq_ary});
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
    # Check that we don't exceed max intron size
    return undef if($best_target->{intron_size} > $self->max_intron_size);

    # Get the aln after the splice junction as HSPs
    my $b1_rel_offset = $best_target->{d} + 1;
    my ($hsps, $subscore) = $self->_dpAlignment_to_hsps($best_aln, $b1_rel_offset);

    # Was the score after the splice junction above the threshold?
    return undef if($subscore < $self->aln_score_threshold);
    # Check that border2 is really covered by the alignment
    return undef if(@$hsps == 0 or $hsps->[0]->[2] != $b1_rel_offset);

    return ($hsps, $best_target, $subscore);
}


sub _dp_align {
    my($self, $q, $t, $dir) = @_;

    #my @Q = ('-', split(//, $q)); # query sequence to be aligned with
    #my @T = ('-', split(//, $t)); # target sequence
    my @Q = ('-', @$q); # query sequence to be aligned with
    my @T = ('-', @$t); # target sequence
    print STDERR @Q, " vs\n",@T,"\n";
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
	    if($s1 >= $s2) {
		if($s1 >= $s3) { $k = 0; $s = $s1; }
		else { $k = 2; $s = $s3; }
	    }
	    else {
		if($s2 >= $s3) { $k = 1; $s = $s2; }
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
    $self->_debug_print_aln(\@Q, \@T, \@N, \@M, $best_i, $best_j);

    return ($M[$best_i][$best_j], [\@M, \@N, $best_i, $best_j]);
}

# do: find out why fewer ends moved
#       probably bcs subscore calculation affected
# do: test if hsps generated are the same

sub _align {
    my($self, $q, $t) = @_;

    #print STDERR @$q, " vs\n",@$t,"\n";
    my (@M, @N); # score and traceback matrices

    # Get sequence lengths
    my $q_len = @$q;
    my $t_len = @$t;
    my $min_seq_len = $q_len < $t_len ? $q_len : $t_len;

    # Get score parameters
    my $match_score = $self->match_score();
    my $mismatch_score = $self->mismatch_score();
    my $gap_score = $self->gap_score();
    my $extension_threshold = $self->extension_threshold();
    my $minus_infinity = -1000;

    # Extend as far as there are exact matches
    my $offset = 0;
    $offset++ while($offset < $min_seq_len and $q->[$offset] eq $t->[$offset]);
    my $best_score = $offset * $match_score;
    my ($best_i, $best_j) = (0,0);
    $M[0][0] = $best_score;
    $N[0][0] = $offset;
    $offset--;

    #print STDERR " Exact match up to $best_i; score = $best_score\n";

    # Compute scores for exact alignment followed by only gaps
    my $min_score = $best_score - $extension_threshold;
    my $i = 1;
    for(my $s = $best_score+$gap_score;	$s >= $min_score; $s += $gap_score) {
	$M[$i][0] = $M[0][$i] = $s;
	$N[$i][0] = 1;
	$N[0][$i] = 2;
	#print STDERR " gap $i; score = $s\n";
	$i++;
    }
    my $first_infinity = $i;
    for my $j ($first_infinity..$t_len) { $M[0][$j] = $minus_infinity;  }
      
    # DP align
    my $end_row = $q_len-$offset-1;
    my ($next_start_col, $next_end_col) = (1, $first_infinity);
    for my $i (1..$end_row) {
	my $q_nt = $q->[$i+$offset];
	my ($start_col, $end_col) = ($next_start_col, $next_end_col);
	$next_start_col = undef;
	$M[$i][$start_col-1]=$minus_infinity unless(defined $M[$i][$start_col-1]);
	for my $j ($start_col..$t_len-$offset-1) {
	    my ($k, $s);
	    my $p =  ($q_nt eq $t->[$j+$offset]) ? $match_score : $mismatch_score;
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
	    if($s >= $min_score) {
		$next_start_col = $j unless(defined $next_start_col);
		$next_end_col = $j;
		$N[$i][$j] = $k;
		$M[$i][$j] = $s;
		if($s > $best_score) {
		    $best_score = $s;
		    ($best_i, $best_j) = ($i, $j); 
		}
	    }
	    else {
		if($j < $end_col) {
		    $M[$i][$j]=$minus_infinity;
		}
		else {
		    for my $k ($j..$t_len) { $M[$i][$k]=$minus_infinity; }
		    last;
		}
	    }
	}
	last unless($next_start_col);
	$next_end_col++;
    }

    # Output alignment (for DEBUG)
    #$self->_debug_print_matrix($q, $t, \@N, \@M);
    #$self->_debug_print_aln($q, $t, \@N, \@M, $best_i, $best_j) if($self->debug);

    return ($M[$best_i][$best_j], [\@M, \@N, $best_i, $best_j]);
}


sub _dpAlignment_to_hsps
{
    my ($self, $aln, $tStart) = @_;

    #print STDERR "to_hsps: $tStart\n";

    # Reconstruct alignment as a set of HSPs
    my ($M, $N, $last_i, $last_j) = @$aln;
    my ($i, $j) = ($last_i, $last_j);
    my $offset = $N->[0][0];
    my ($qs,$qe,$ts,$te);
    my @hsps;
    my $min_j = $tStart > $offset ? $tStart - $offset : 1;
    while($i >= 1 and $j >= $min_j) {
	if($N->[$i][$j] == 0) {		# 0 = match
	    ($qs, $ts) = ($i, $j);
	    ($qe, $te) = ($i, $j) unless($qe);
	    $i--; $j--;
	}
	else {
	    if($N->[$i][$j] == 1) {	# 1 = tGap
		$i--;
            }
            else {			# 2 = qGap
		$j--;
	    }
	    if($qe) {			# add hsp if we come from a match
		push @hsps, [$qs+$offset, $qe+$offset, $ts+$offset, $te+$offset];
		$qe = 0;
	    }
	}
    }
    if($qe and $qs == 1 and $ts == 1) {
	push @hsps, [$tStart, $qe+$offset, $tStart, $te+$offset];
    }
    else {
	push @hsps, [$qs+$offset, $qe+$offset, $ts+$offset, $te+$offset] if($qe);
        push @hsps, [$tStart, $offset, $tStart, $offset] if($offset >= $tStart);
    }
    @hsps = reverse(@hsps);
    my $subscore = $offset >= $tStart ?
	$M->[$last_i][$last_j] - ($tStart-1) * $self->match_score :
	$M->[$last_i][$last_j] - $M->[$i][$j];
    #print STDERR "Qry HSPs: ", join(", ", map {$_->[0].'-'.$_->[1]} @hsps), "\n";
    #print STDERR "Tgt HSPs: ", join(", ", map {$_->[2].'-'.$_->[3]} @hsps), "\n";
    #print STDERR "Subscore = $subscore; Totscore = ",$M->[$last_i][$last_j],"\n";
    
    # return hsps
    return (\@hsps, $subscore);
}


sub _create_new_primap
{
    my ($self, $m0, $l_hsps, $r_hsps, $t_seq) = @_;
    
    my @hsps0 = map { [ $_->qStart, $_->qEnd, $_->tStart, $_->tEnd ] } $m0->all_HSPs;

    #print STDERR "\nInitial  HSPs:\n";
    #$self->_debug_print_hsp_coords(\@hsps0);
    #print STDERR "\nLeft HSPs:\n";
    #$self->_debug_print_hsp_coords($l_hsps) if($l_hsps);
    #print STDERR "\nRight HSPs:\n";
    #$self->_debug_print_hsp_coords($r_hsps) if($r_hsps);

    my $hsps = $l_hsps ? $self->_change_left_hsps($l_hsps, \@hsps0) : \@hsps0;
    $hsps = $self->_change_right_hsps($r_hsps, $hsps) if($r_hsps);

    #print STDERR "\nChanged  HSPs:\n";
    #$self->_debug_print_hsp_coords($hsps);
    
    my $m = AT::MappingFactory->create_mapping(HSP_coords => $hsps,
					       target_seq => $t_seq,
					       query_seq => $m0->query_seq,
					       mapping_id => $m0->mapping_id,
					       strand => $m0->strand,
					       target_db => $m0->target_db,
					       tName => $m0->tName,
					       tSize => $m0->tSize,
					       qType => $m0->qType,
					       mRNAInfo => $m0->mRNAInfo);
					       

    return $m;
}


sub _debug_print_hsp_coords
{
    my ($self, $coords) = @_;
    foreach my $hsp (@$coords) {
	print STDERR join("\t", @$hsp), "\n";
    }
}


sub _change_left_hsps
{
    my ($self, $l_hsps_ref, $m_hsps_ref) = @_;    
    
    my @hsps;
    foreach my $hsp (@$l_hsps_ref) {
	push @hsps, [@$hsp];
    }
    my $h = $hsps[-1];

    my $i;
    for ($i=0; $i < @$m_hsps_ref; $i++) {
	last if($m_hsps_ref->[$i][1] > $h->[1]);
    }

    my $first_mid = $m_hsps_ref->[$i];
    if($first_mid and $first_mid->[0] <= $h->[1]) {   		# overlap?
        if($first_mid->[3]-$h->[2] == $first_mid->[1]-$h->[0]) {	# continuation?
	    $h->[2] = $first_mid->[2];
	    $h->[3] = $first_mid->[3];
	}
	else {
	    my $d = $h->[1] - $first_mid->[0] + 1;
	    $h = [$h->[1]+1, $first_mid->[1], $first_mid->[2]+$d, $first_mid->[3]];
	    push @hsps, $h;
	}
	$i++;
    }

    for (; $i < @$m_hsps_ref; $i++) {
        push @hsps, [@{$m_hsps_ref->[$i]}];
    }

    return \@hsps;
}


sub _change_right_hsps
{
    my ($self, $r_hsps_ref, $m_hsps_ref) = @_;    
    
    my @hsps;
    my @h = @{$r_hsps_ref->[0]};

    my $i;
    for ($i = 0; $i < @$m_hsps_ref; $i++) {
	last if($m_hsps_ref->[$i][1] >= $h[0]);
        push @hsps, [@{$m_hsps_ref->[$i]}];
    }

    my $last_mid = $m_hsps_ref->[$i];
    if($last_mid and $last_mid->[0] < $h[0]) {   		# starts before?
        if($h[3]-$last_mid->[2] == $h[1]-$last_mid->[0]) {	# continuation?
	    $h[0] = $last_mid->[0];
	    $h[1] = $last_mid->[1];
	}
	else {
	    my $d = $h[0] - $last_mid->[0] - 1;
	    push @hsps, [$last_mid->[0], $h[0]-1, $last_mid->[2], $last_mid->[2]+$d];
	}
    }

    push @hsps, \@h;

    for my $j (1.. @$r_hsps_ref-1) {
        push @hsps, [@{$r_hsps_ref->[$j]}];
    }

    return \@hsps;
}


sub _debug_print_aln
{
    my ($self, $Q, $T, $N, $M, $last_i, $last_j) = @_;

    # Reconstruct alignment    
    my @qAli;# = ($Q->[$i]);
    my @tAli;# = ($T->[$j]);
    my $offset = $N->[0][0]-1;
    my ($i, $j) = ($last_i, $last_j);
    while($i > 0 or $j > 0) {
	if($N->[$i][$j] == 0) {		# 0 = match
	    push @qAli, $Q->[$i+$offset];
	    push @tAli, $T->[$j+$offset];
	    $i--; $j--;
	}
	elsif($N->[$i][$j] == 1) {	# 1 = tGap
	    push @qAli, $Q->[$i+$offset];
	    push @tAli, '-';
	    $i--;
	}
	else {   			# 2 = qGap
	    push @qAli, '-';
	    push @tAli, $T->[$j+$offset];
	    $j--;
        }
    }
    while($offset >= 0) {
	push @qAli, $Q->[$offset];
	push @tAli, $T->[$offset];
	$offset--;
    }

    if(@qAli == 0) {
	print STDERR "(no alignment)\n";
	return;
    }
    
    # Output to STDERR in ClustalW format
    my $qSeq = Bio::LocatableSeq->new(-id => 'query',
				      -seq => join('',reverse @qAli),
				      -start => 1,
				      -end => $last_i+$N->[0][0]);
    my $tSeq = Bio::LocatableSeq->new(-id => 'target',
				      -seq => join('',reverse @tAli),
				      -start => 1,
				      -end => $last_j+$N->[0][0]);
    my $aln = Bio::SimpleAlign->new();
    $aln->add_seq($qSeq); $aln->add_seq($tSeq);
    my $aln_out = Bio::AlignIO->new(-fh => \*STDERR, -format => 'clustalw');
    $aln_out->write_aln($aln);
    #print STDERR "Score = ",$M->[$last_i][$last_j],"\n";
#    print STDERR "Query:\t", reverse(@qAli), "\nTarget:\t", reverse(@tAli), "\n";
} 


sub _debug_print_matrix
{
    my ($self, $Q, $T, $N, $M) = @_;
    my $offset = $N->[0][0];
    print STDERR "Q/T   *";
    for(my $i = $offset; $i < @$T; $i++) { print STDERR "     ",$T->[$i]; }
    print STDERR "\n";
    for my $i (0..@$Q-$offset) {
	if($i) { print STDERR $Q->[$i+$offset-1]; }
	else { print STDERR "*"; }
	for my $j (0..@$T-$offset) {
	    my $score = $M->[$i][$j];
	    if(defined $score) { print STDERR sprintf("%6d",$score); }
	    else { print STDERR "     x"; }
	}
	print STDERR "\n\n";
    }
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
    foreach my $i (FWD_IDX, REV_IDX) {
      print STDERR "---\nword2starts index $i\n";
      $idx = $self->{_word2starts}[$i];
      foreach my $word (sort {$a cmp $b} keys %$idx) {
	print STDERR $word, "\t", join(",",@{$idx->{$word}}), "\n";
      }
      print STDERR "---\nword2ends index $i\n";
      $idx = $self->{_word2ends}[$i];
      foreach my $word (sort {$a cmp $b} keys %$idx) {
	print STDERR $word, "\t", join(",",@{$idx->{$word}}), "\n";
      }
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


# do: introduce offset for length of exact alignment; put score in M[0][0]; put offset in N[0][0]

#sub _align {
#    my($self, $q, $t) = @_;
#
#    my @Q = ('-', @$q); # query sequence to be aligned with
#    my @T = ('-', @$t); # target sequence
#    print STDERR @Q, " vs\n",@T,"\n";
#    my (@M, @N); # score and traceback matrices
#    my $best_score = 0;
#
#    # Get score parameters
#    my $match_score = $self->match_score();
#    my $mismatch_score = $self->mismatch_score();
#    my $gap_score = $self->gap_score();
#    my $extension_threshold = $self->extension_threshold();
#    my $minus_infinity = -1000;
#
#    # Extend as far as there are exact matches
#    $M[0][0] = 0; $N[0][0] = 0;
#    my $i = 1;
#    while($i < @Q and $i < @T and $Q[$i] eq $T[$i]) {
#	$best_score += $match_score;
#	$M[$i][$i] = $best_score;
#	$N[$i][$i] = 0;
#	$i++;
#    }
#    my $first_mismatch = $i;  # index of first_mismatch;
#    my ($best_i, $best_j) = ($i-1,$i-1); # end position of best alignment so far
#
#    #print STDERR " Exact match up to $best_i; score = $best_score\n";
#
#    # Compute scores for exact alignment followed by only gaps
#    my $min_score = $best_score - $extension_threshold;
#    for(my $s = $best_score+$gap_score;	$s >= $min_score; $s += $gap_score) {
#	$M[$i][$best_i] = $M[$best_i][$i] = $s;
#	$N[$i][$best_i] = 1;
#	$N[$best_i][$i] = 2;
#	#print STDERR " gap $i; score = $s\n";
#	$i++;
#    }
#    my $first_infinity = $i;
#
#    # Fill rest of matrix with minus_infinity
#    for my $j ($first_infinity..@T-1) { $M[$best_i][$j] = $minus_infinity;  }
#      
#    # DP align
#    my ($next_start_col, $next_end_col) = ($first_mismatch, $first_infinity);
#    for my $i ($first_mismatch..@Q-1) {
#	my ($start_col, $end_col) = ($next_start_col, $next_end_col);
#	$next_start_col = undef;
#	$M[$i][$start_col-1]=$minus_infinity unless(defined $M[$i][$start_col-1]);
#	for my $j ($start_col..@T-1) {
#	    my ($k, $s);
#	    my $p =  ($Q[$i] eq $T[$j]) ? $match_score : $mismatch_score;
#	    my $s1 = $M[$i-1][$j-1] + $p;
#	    my $s2 = $M[$i-1][$j] + $gap_score;
#	    my $s3 = $M[$i][$j-1] + $gap_score;
#	    if($s1 >= $s2) {
#		if($s1 >= $s3) { $k = 0; $s = $s1; }
#		else { $k = 2; $s = $s3; }
#	    }
#	    else {
#		if($s2 >= $s3) { $k = 1; $s = $s2; }
#		else { $k = 2; $s = $s3; }
#	    }
#	    if($s >= $min_score) {
#		$next_start_col = $j unless(defined $next_start_col);
#		$next_end_col = $j;
#		$N[$i][$j] = $k;
#		$M[$i][$j] = $s;
#		if($s > $best_score) {
#		    $best_score = $s;
#		    ($best_i, $best_j) = ($i, $j); 
#		}
#	    }
#	    else {
#		if($j < $end_col) {
#		    $M[$i][$j]=$minus_infinity;
#		}
#		else {
#		    for my $k ($j..@T-1) { $M[$i][$k]=$minus_infinity; }
#		    last;
#		}
#	    }
#	}
#	last unless($next_start_col);
#	$next_end_col++;
#    }
#
#    # Output alignment (for DEBUG)
#    #$self->_debug_print_matrix(\@Q, \@T, \@M);
#    $self->_debug_print_aln(\@Q, \@T, \@N, \@M, $best_i, $best_j);
#
#    return ($M[$best_i][$best_j], [\@M, \@N, $best_i, $best_j]);
#}



1;
