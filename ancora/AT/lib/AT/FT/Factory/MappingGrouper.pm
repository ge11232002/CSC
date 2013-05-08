# AT::FT::Factory::MappingGrouper module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

package AT::FT::Factory::MappingGrouper;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;

@ISA = qw/AT::Root/;


sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
sing_mul_maxJncGlitch_bp => ($args{sing_mul_maxJncGlitch_bp} || 10),
mul_mul_maxJncGlitch_bp => ($args{mul_mul_maxJncGlitch_bp} || 0),
sing_sing_minOl_bp => ($args{sing_sing_minOl_bp} || 50),
sing_sing_minOl_pc => ($args{sing_sing_minOl_pc} || 50),
sing_mulExt_minOl_bp => ($args{sing_mulExt_minOl_bp} || 50),
sing_mulExt_minOl_pc => ($args{sing_mulExt_minOl_pc} || 50),
mul_mul_minJncMatches_abs => ($args{mul_mul_minJncMatches_abs} || 1),
mul_mul_minJncMatches_pc => ($args{mul_mul_minJncMatches_pc} || 50),
verbose => ( $args{verbose} || 0 )
		      
    }, ref $caller || $caller;
    
    return $self;
}


sub run
{
    my ($self, %args) = @_;
    my $mappings = $args{mappings} || croak "No mappings arg";
    #my $target_seq = $args{target_seq} || croak "No target_seq arg";
    my $strand = $args{strand} || croak "No strand arg";

    my $groups = $self->_group_mappings($mappings, $strand);

    return $groups;
}


sub _group_mappings {
    my ($self, $mappings, $strand) = @_;

    # Sort in incr exon order to make sure that seed groups
    # with unspliced before seeding groups with spliced.
    # It is probably also an optimization to seed groups with
    # mappings with a large number of exons.
    $mappings = [ sort
		  { $a->nr_unbroken_gaps <=> $b->nr_unbroken_gaps }
		  @$mappings ];
    my $nr_mappings = scalar(@$mappings);

    if($self->debug) {
        print STDERR "\nGrouping mappings: ";
        foreach my $i (0..@$mappings-1) {
            print STDERR $i.'.'.$mappings->[$i]->qName_string.'; ';
        }
        print STDERR "\n\n";
    }

    # strategy: seed group with a mapping, try to include others
    # We don't need to reassess the similarity between two
    # mappings that we have already compared. Use @similarity to store
    # this info.
    my @group_lists;
    my @grouped = (0) x $nr_mappings;
    my $ungrouped = @$mappings-1; # index of some mapping yet ungrouped
    my @similarity; $#similarity = $nr_mappings-1;
    while($ungrouped != -1) {	# loop while there are ungrouped mappings
	#print STDERR "group seed: $ungrouped\n";
	my @group_flg = (0) x scalar(@$mappings);
	my @group_lst = ($mappings->[$ungrouped]);
	$group_flg[$ungrouped] = 1; # start group with the ungrouped
	$grouped[$ungrouped] = 1;
	my $grouped_one = 1;
	while($grouped_one and $ungrouped != -1) {
	    $grouped_one = 0;
	    $ungrouped = -1;    # assume we are out of mappings to group
	    for(my $i = 0; $i < @$mappings; $i++) {
		if($group_flg[$i]) {
		    # mapping is already in this group - do nothing
		    #print STDERR "$i: present\n";
		}
		elsif(($grouped[$i] == 0
			#or @{$mappings->[$i]->exons} == 1
			) and
		      $self->_fits_in_group($mappings->[$i], \@group_lst, $strand,
					    $i, \@similarity)) {
		    # mapping is ungrouped and fits in this group
		    push @group_lst, $mappings->[$i];
		    $group_flg[$i] = 1;
		    $grouped[$i] = 1;
		    $grouped_one = 1;
		    #print STDERR "$i: added\n";
		}
		elsif($grouped[$i] == 0) {
		    # mapping is ungrouped but does not fit in this group
		    $ungrouped = $i;
		    #print STDERR "$i: ungrouped\n";
		}
	    }
	}
	#print STDERR (join ' ', @group_flg), " ",
	#    (join ",", map {$_->mapping->qName} @group_lst),"\n";
	push @group_lists, \@group_lst;
    }

    return \@group_lists;
}


sub _fits_in_group {
    my ($self, $candidate, $group, $strand, $i, $sim) = @_;
    foreach my $member (@$group) {
	$sim->[$i] = {} unless(defined $sim->[$i]);
	my $result;
	unless(defined($result = $sim->[$i]->{$member})) {
	    $result = $sim->[$i]->{$member} = 
		$self->_mappings_compatible($candidate, $member, $strand);
	}
        return 1 if ($result);
    }
    return 0;
}


sub _mappings_compatible {
    my ($self, $a, $b, $strand) = @_;

    # note: a and b are treated differently
    # a should be the candidate, b the member

    my $hspsA = $a->HSP_ref;
    my $hspsB = $b->HSP_ref;

    # return false if they don't overlap at all
    if($hspsA->[0]->start > $hspsB->[-1]->end or $hspsA->[-1]->end < $hspsB->[0]->start)
    { return 0 };

    my $a_is_spliced = $a->nr_unbroken_gaps;
    my $b_is_spliced = $b->nr_unbroken_gaps;

    if(!$a_is_spliced and !$b_is_spliced) {
	# both single-exon
	# compatible if overlap is at least 50 bp or 50% of the smaller
	foreach my $hspA (@$hspsA) {
	    foreach my $hspB (@$hspsB) {
		return 1 if $self->_overlap_size_check
		    ($hspA->start, $hspA->end,
		     $hspB->start, $hspB->end,
		     $self->sing_sing_minOl_bp,
		     $self->sing_sing_minOl_pc);
	    }
	}
    }
    elsif(!$a_is_spliced) {
	# a is single-exon, b is multi-exon
	# Compatible if
	# (a) a overlaps with an exteral exon of b and does not extend into
	#     the gene past that exon. The overlap should be at least x bp
	#     or y% of the smallest of the two exons.
	# (b) a overlaps with an internal exon of b and does not extend
	#     past the borders of that exon
        my $conf_glitch = $self->sing_mul_maxJncGlitch_bp;
	foreach my $hspA (@$hspsA) {
	    my($s1, $e1) = ($hspA->start, $hspA->end);
	    if($e1 <= $hspsB->[0]->end+$conf_glitch) {
		return $self->_overlap_size_check($s1, $e1,
					      $hspsB->[0]->start, $hspsB->[0]->end,
					      $self->sing_mulExt_minOl_bp,
					      $self->sing_mulExt_minOl_pc);
	    }
	    elsif($s1 >= $hspsB->[-1]->start-$conf_glitch) {
		return $self->_overlap_size_check($s1, $e1,
					      $hspsB->[-1]->start, $hspsB->[-1]->end,
					      $self->sing_mulExt_minOl_bp,
					      $self->sing_mulExt_minOl_pc);
	    }
	    else {
	        foreach my $i (1 .. @$hspsB-1) {
		    if($s1 >= $hspsB->[$i]->start-$conf_glitch and
		        $e1 <= $hspsB->[$i]->end+$conf_glitch) {
		        return 1;
		    }
		}
	    }
	}
    }
    elsif(!$b_is_spliced) {
	# a is multi-exon, b is single-exon
	# Never compatible.
	# We don't allow single-exon genes to bring in multi-exon genes.
	# Otherwise large unspliced transcripts can cause unrelated spliced
	# transcripts to be grouped together.
	return 0;
    }
    else {
	# both multi-exon
	# Compatible if they share at least x or n/y splice junctions,
	# where n is the number of splice junctions for the transcript with
	# the fewest exons.
	#print STDERR "comp-check: ",$a->{mapping}->qName," ", $b->{mapping}->qName, "\n";
        my $spljnc_glitch = $self->mul_mul_maxJncGlitch_bp;
	my $shared_spljnc = 0;
	my($i, $j) = (0,0);
	while($i < @$hspsA and $j < @$hspsB) {
	    my($s1, $e1) = ($hspsA->[$i]->start, $hspsA->[$i]->end);
	    my($s2, $e2) = ($hspsB->[$j]->start, $hspsB->[$j]->end);
	    #print STDERR $i."[$s1-$e1]:".$j."[$s2-$e2]  ";
	    if ($i and $j and
		!$hspsA->[$i]->left_gap->is_broken and
		!$hspsB->[$j]->left_gap->is_broken and
		$s1 >= $s2 - $spljnc_glitch and $s1 <= $s2 + $spljnc_glitch) {
		$shared_spljnc++;
	    }
	    if($i < @$hspsA-1 and $j < @$hspsB-1 and
		!$hspsA->[$i]->right_gap->is_broken and
		!$hspsB->[$j]->right_gap->is_broken and
		$e1 >= $e2 - $spljnc_glitch and $e1 <= $e2 + $spljnc_glitch) {
		$shared_spljnc++;
	    }
	    if ($e1 < $e2) {
		$i++;
	    }
	    else {
		$j++;
	    }
	}
	return 1 if ($shared_spljnc > $self->mul_mul_minJncMatches_abs);
	my $min_exons = (@$hspsA < @$hspsB) ? @$hspsA : @$hspsB;
	my $min_spljnc =
	    $self->mul_mul_minJncMatches_pc * ($min_exons * 2 - 2) / 100;
	return ($shared_spljnc >= $min_spljnc) ? 1 : 0;
    }
}


sub _overlap_size_check
{
    my ($self, $s1, $e1, $s2, $e2, $min_bp, $min_pc) = @_;
    my $ol_start = ($s1 > $s2) ? $s1 : $s2;
    my $ol_end = ($e1 < $e2) ? $e1 : $e2;
    my $ol = $ol_end - $ol_start;
    return 1 if ($ol >= $min_bp);
    my $size1 = $e1 - $s1 + 1;
    my $size2 = $e2 - $s2 + 1;
    my $min_size = ($size1 < $size2) ? $size1 : $size2;
    return ($ol >= $min_size * $min_pc / 100) ? 1 : 0;
}


1;
