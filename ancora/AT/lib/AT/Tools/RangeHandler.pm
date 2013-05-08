# AT::Tools::RangeHandler module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::Tools::RangeHandler - various methods for coordinate range manipulation

=head1 SYNOPSIS

 # Create two sets of ranges to play with
 my @r1 = ([100, 200], [300,400], [500,600]);
 my @r2 = ([150, 350], [550,700]);

 # Compute intersection
 my $i = AT::Tools::RangeHandler->compute_intersection(\@r1, \@r2);
 print join(', ', map {$_->[0].'-'.$_->[1]} @$i), "\n";

 # Compute union
 my $u = AT::Tools::RangeHandler->compute_union(\@r1, \@r2);
 print join(', ', map {$_->[0].'-'.$_->[1]} @$u), "\n";

=head1 OVERVIEW

This class should be made into a class AT::RangeSet that would be an array.
All standart array methods (push,pop) should be allowed, but there should
also be the extra methods in this class (and more?).

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::Tools::RangeHandler;

use strict;
use vars '@ISA';
use AT::Root;
use Carp;

@ISA = qw(AT::Root);


=head2 sum

 Title     : sum
 Usage     : my $sum = AT::Tools::RangeHandler->
		sum(\@ranges);
 Function  : Sum the sizes of all ranges in a set.
 Returns   : Integer
 Args      : Reference to two 2D arrays of coordinate pairs.

=cut


sub sum {
    my ($self, $ranges) = @_;
    my $sum = 0;
    foreach my $range (@$ranges) {
	$sum += $range->[1]-$range->[0]+1;
    }
    return $sum;
}


=head2 span

 Title     : span
 Usage     : my $span = AT::Tools::RangeHandler->span(\@ranges);
 Function  : Compute the span of a range set, i.e.
             the distance between the start of the first range
             and the end of the last range.
 Returns   : Integer
 Args      : Reference to 2D array of coordinate pairs

=cut


sub span {
    my ($self, $ranges) = @_;
    return $ranges->[-1][1] - $ranges->[0][0] + 1;
}


=head2 merge_adjacent

 Title     : merge_adjacent
 Usage     : my $merged = AT::Tools::RangeHandler->
		merge_adjacent(\@original);
 Function  : Merge adjacent ranges in a range set.
             Ranges a and b are adjacent if end(a)+1 == start(b)
 Returns   : Reference to 2D array of coordinate pairs
 Args      : Reference to 2D array of coordinate pairs

=cut


sub merge_adjacent {
    my ($self, $r_in) = @_;
    my @r_out;
    my ($start1,$end1) = @{$r_in->[0]};
    for my $i (1..@$r_in-1) {
	my ($start2, $end2) = @{$r_in->[$i]};
	if($start2 == $end1+1) {
	   $end1 = $end2;
	}
	else {
	    push @r_out, [$start1, $end1];
	    ($start1, $end1) = ($start2, $end2);
	}
    }
    push @r_out, [$start1,$end1];
    return \@r_out;
}




=head2 compute_inverse

 Title     : compute_inverse
 Usage     : my $intergenic = AT::Tools::RangeHandler->
		compute_inverse(\@gene_bounds, 1, $chr_size);
 Function  : Compute inverse of a set of ranges
 Returns   : Reference to 2D array of coordinate pairs
 Args      : 1. References to 2D array of coordinate pairs
             2. Start bound (optional)
             3. End bound (optional)

=cut


sub compute_inverse {
    my ($self, $r_ref, $start_bound, $end_bound) = @_;
    my @r = sort { $a->[0] <=> $b->[0] } @$r_ref;
    my @inverse;

    $start_bound = $r[0][1]+1 unless($start_bound);
    $end_bound = $r[-1][0]-1 unless($end_bound);
    my $prev_end = 0;
    foreach my $range (@r) {
	my ($start,$end) = @$range;
	my $inv_start = $prev_end >= $start_bound ? $prev_end+1 : $start_bound;
	my $inv_end = $start <= $end_bound ? $start-1 : $end_bound;
	push @inverse, [$inv_start,$inv_end] if($inv_start <= $inv_end);
	$prev_end = $end;
	last if($prev_end >= $end_bound);
    }

    if($prev_end <= $end_bound) {
	my $inv_start = $prev_end >= $start_bound ? $prev_end+1 : $start_bound;
	push @inverse, [$inv_start, $end_bound];
    }

    return \@inverse;
}


=head2 compute_intersection

 Title     : compute_intersection
 Usage     : my $int = AT::Tools::RangeHandler->
		compute_intersection(\@ranges1, \@ranges2);
 Function  : Compute intersection of two sets of ranges.
 Returns   : Reference to 2D array of coordinate pairs
 Args      : References to two 2D arrays of coordinate pairs.

=cut


sub compute_intersection {
    my ($self, $r1_ref, $r2_ref) = @_;
    my @r1 = sort { $a->[0] <=> $b->[0] } @$r1_ref;
    my @r2 = sort { $a->[0] <=> $b->[0] } @$r2_ref;
    my @intersections;

    my($i, $j) = (0, 0);

    while($i < @r1 and $j < @r2) {
        my $max_start = ($r1[$i][0] > $r2[$j][0]) ?
	    $r1[$i][0] : $r2[$j][0];
	my $min_end = ($r1[$i][1] < $r2[$j][1]) ?
	    $r1[$i][1] : $r2[$j][1];
	push @intersections, [$max_start, $min_end] if($min_end >= $max_start);
	if ($r1[$i][1] < $r2[$j][1]) {
	    $i++;
	}
	else {
	    $j++;
	}
    }
    return \@intersections;
}


=head2 compute_union

 Title     : compute_union
 Usage     : my $union = AT::Tools::RangeHandler->
		compute_union(\@ranges1, \@ranges2, ...);
 Function  : Compute union of two or more sets of ranges.
 Returns   : Reference to 2D array of coordinate pairs
 Args      : References to two 2D arrays of coordinate pairs.

=cut


sub compute_union {
    my ($self, @range_refs) = @_;
    my @r = sort { $a->[0] <=> $b->[0] } map { @$_ } @range_refs;
    return [] unless(@r);
    my @u = ([@{$r[0]}]);
    foreach my $range (@r[1..@r-1]) {
	if($range->[0] <= $u[-1][1]) {
	    $u[-1][1] = $range->[1] if($range->[1] > $u[-1][1]);
	}
	else {
	    push @u, [@$range];  
	}
    }
    return \@u;
}


=head2 partition_by_overlap

 Title     : partition_by_overlap
 Usage     : my ($overlapping, $nonoverlapping)
              = AT::Tools::RangeHandler->
		partition_by_overlap(\@ranges1, \@ranges2);
 Function  : Partition the ranges in one set based on whether or
             not they overlap with ranges in a second set.
             Note that range datastructures are not duplicated,
             i.e. the returned range sets will combined contain
             the same references as the first input range set.
             TODO: Generalize so that we can give an overlap
             threshold.
 Returns   : References to two 2D arrays of coordinate pairs
 Args      : References to two 2D arrays of coordinate pairs.

=cut


sub partition_by_overlap {
    my ($self, $r1_ref, $r2_ref) = @_;
    my @r1 = sort { $a->[0] <=> $b->[0] } @$r1_ref;
    my @r2 = sort { $a->[0] <=> $b->[0] } @$r2_ref;
    my (@overlapping, @nonoverlapping);

    my($i, $j) = (0, 0);

    while($i < @r1 and $j < @r2) {
        # If a ends before b starts, put a in nonoverlapping, get next a
	if($r1[$i][1] < $r2[$j][0]) {
	    push @nonoverlapping, $r1[$i];
	    $i++;
	}
        # elsif b ends before a starts, get next b
	elsif($r2[$j][1] < $r1[$i][0]) {
	    $j++;
	}
        # else put a in overlapping, get next a
	else {
	    push @overlapping, $r1[$i];
	    $i++;
	}
    }
    push @nonoverlapping, @r1[$i..@r1-1];
    return (\@overlapping, \@nonoverlapping);
}


sub partition_by_overlap2 {
# a twist on the above method
    my ($self, $r1_ref, $r2_ref) = @_;
    my @r1 = sort { $a->[0] <=> $b->[0] } @$r1_ref;
    my @r2 = sort { $a->[0] <=> $b->[0] } @$r2_ref;
    my @r3 = map { [] } (1..@r2);

    my($i, $j) = (0, 0);

    while($i < @r1 and $j < @r2) {
        # If a ends before b starts, get next a
	if($r1[$i][1] < $r2[$j][0]) {
	    $i++;
	}
        # elsif b ends before a starts, get next b
	elsif($r2[$j][1] < $r1[$i][0]) {
	    $j++;
	}
        # else (overlapping): put a in overlapping,
	# get next for that which ends first
	else {
	    push @{$r3[$j]}, [@{$r1[$i]}];
	    if($r1[$i][1] < $r2[$j][1]) {
		$i++;
	    }
	    else {
		$j++;
	    }
	}
    }
    return \@r3;
}




#
# Got this sub from the antisense project
# may want to add it as a method
#

sub get_unique_regions
{
    my ($regions1_ref, $regions2_ref) = @_;   
    my @regions1 = @$regions1_ref;
    my @regions2 = @$regions2_ref;
    my @uni1;
    my @uni2;

    my ($start1, $end1) = @{shift @regions1 || []};
    my ($start2, $end2) = @{shift @regions2 || []};

    while($start1 and $start2) {

	my $overlap = 1;

	# Get a unique region if starts differ
	if($start1 != $start2) {
	    if($start1 < $start2) {
		if($end1 < $start2) {
		    push @uni1, [$start1, $end1];
		    ($start1, $end1) = @{shift @regions1 || []}; # shift both
		    $overlap = 0;
		}
		else {
		    push @uni1, [$start1, $start2-1];
		}
	    }
	    else {
		if($end2 < $start1) {
		    push @uni2, [$start2, $end2];
		    ($start2, $end2) = @{shift @regions2 || []};
		    $overlap = 0;
		}
		else {
		    push @uni2, [$start2, $start1-1];
		}
	    }
	}

	# Clip and shift regions if we found some overlap
	if($overlap) {
	    if($end1 == $end2) {		#   if ends equal
	        ($start1, $end1) = @{shift @regions1 || []}; # shift both
	        ($start2, $end2) = @{shift @regions2 || []};
	    }
	    elsif($end1 < $end2) {		#   elsif 1 ends before 2
		$start2 = $end1+1;		#     start2 = end1+1
	        ($start1, $end1) = @{shift @regions1 || []}; # shift 1
	    }
	    else {				#   else (1 ends after 2)
		$start1 = $end2+1;		#     start1 = end2+1
	        ($start2, $end2) = @{shift @regions2 || []}; # shift 2
	    }
	}
    }

    if($start1) {
	push @uni1, [$start1, $end1];
	foreach my $r (@regions1) {
	    push @uni1, [@$r];
	}
    }
    elsif($start2) {
	push @uni2, [$start2, $end2];
	foreach my $r (@regions2) {
	    push @uni2, [@$r];
	}
    }
   
    return (\@uni1, \@uni2);
}










;

