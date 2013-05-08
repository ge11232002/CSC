# AT::GFX::AlignedLoci module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::GFX::Locus - module for drawing loci with features in AT::DB::GenomeFeature

=head1 SYNOPSIS


=cut


package AT::GFX::AlignedLoci;

use strict;
use vars '@ISA';
use Carp;
use Class::Struct;
use AT::Root;
use Bio::Graphics::Panel;
use Bio::SeqFeature::Generic;
use Benchmark;


@ISA = qw(AT::Root);


struct(MyGene => [ id => '$',
		   chr => '$',
		   start => '$',
		   end => '$',
		   strand => '$',
		   exons => '@' ]);

struct(MyLoc => [ asm_param => '$',
		  chr => '$',
		  start => '$',
		  end => '$',
		  strand => '$',
		  aln_regions => '@',
		  offset => '$',
		  image => '$' ]);


# In:
# 1-3 organism, db, tracks sets
# 
# To draw_locus:
# Chr,s,e for 1-3 organisms
# 
# If 0 : die
# If 1/2 : as now
# If 1/3 : middle->out
# If 2/3 : sandwich
# If 3/3 : in order, all fixed

sub new
{
    my ($caller, %args) = @_;

    my $svg = $args{svg};
    my $gd_pkg = $svg ? 'GD::SVG' : 'GD';
    eval "use $gd_pkg";

    my $self = bless {
		      asm_param => ($args{assembly_param} || die "No assembly parameters"),
		      width => ($args{width} || 1250),
		      pad_left => ($args{pad_left} || 30),
		      pad_right => ($args{pad_right} || 220),
		      pad_top => ($args{pad_top} || 10),
		      pad_bottom => ($args{pad_bottom} || 10),
		      max_nr_secondary_panels => defined($args{max_nr_secondary_panels}) ? $args{max_nr_secondary_panels} : -1,
		      connector_height => ($args{connector_height} || 100),
		      min_context => ($args{min_context} || 0),
		      fill_with_context => defined($args{fill_with_context}) ? $args{fill_with_context} : 1,
		      max_space_factor => $args{max_space_factor} || 1.5, 
		      panel => undef,
		      svg => $svg,
		      _image_pkg => $gd_pkg.'::Image',
		      _polygon_pkg => $gd_pkg.'::Polygon',
		      _font_pkg => $gd_pkg.'::Font',
		      _errstr => ''
		     }, ref $caller || $caller;
    
    return $self;   
}


sub errstr {
    return shift->{_errstr};
}


sub draw_locus
{
    my ($self, @args) = @_;
    my $coord_sets;

    if($args[0] and ref($args[0])) {
	$coord_sets = shift @args;
    }
    elsif(@args >= 3) {
	my $chr = shift @args;
	my $start = shift @args;
	my $end = shift @args;
	$coord_sets = [[0, $chr, $start, $end]];
    }
    else {
	croak("Missing arguments to method draw_locus");
    }
      
    # if we want to add some options we can here do:
    # my %options = @args;

    my $nr_coord_sets = @$coord_sets;
    my $nr_asm = @{$self->{asm_param}};

    my $draw_method;

    if($nr_coord_sets == 1) {
	if($nr_asm == 2) { 
	    $draw_method = \&draw_pair;
	}
	elsif($nr_asm == 3) {
	    $draw_method = \&draw_outward_triplet;
	}
    }
    elsif($nr_coord_sets == 2) {
	if($nr_asm == 3) {
	    $draw_method = \&draw_inward_triplet;
	}
    }

    if($draw_method) {
	return $draw_method->($self, $coord_sets);
    }
    else {
	$self->{_errstr} = "Drawing with constraints on $nr_coord_sets out of $nr_asm not implemented";
	return 0;
    }
}


# do:
# modify draw_pair to work in two directions
# decouple the part below _get_aln_data so that it can be invoked by others as "escape route" if alignments are missing
# or decoupled

sub draw_pair
{
    my ($self, $coord_sets) = @_;

    my $rgn_width = $self->{width};
    my $conn_height = $self->{connector_height};

    # Get specs for the reference assembly
    my ($ref_index, $chr, $rgn_start, $rgn_end) = @{$coord_sets->[0]};
    my $ref_loc = MyLoc->new(asm_param => $self->asm_param->[$ref_index],
			     chr => $chr, start => $rgn_start, end => $rgn_end, strand => '+',
			     offset => 0);
    
    # Check that index values are valid
    croak "Undefined reference index\n" unless(defined($ref_index));
    croak "Invalid reference index $ref_index\n" unless($ref_index == 0 or $ref_index == 1);

    # Determine which assembly is the reference
    my $sec_index = 1 - $ref_index;

    # Set scale
    my $scale = ($rgn_end - $rgn_start + 1) / $rgn_width;
    print STDERR "using scale: $scale bp/pixel\n";

    # Get ranked alignment groups, aligned ranges and offsets for each reference
    my ($aln_groups, $sec_loc) = $self->_get_aln_data($sec_index, $ref_loc);

    # Increment offsets so the lowest is zero, and scale the offsets
    my @all_loc = ($ref_loc, @$sec_loc);
    $self->_adjust_offsets($scale, \@all_loc);

    # Add context secondary panels if there is space
    $self->_fill_with_context($ref_loc, $sec_loc, $scale) if($self->{fill_with_context});

    # Draw panels
    $self->_draw_panels(\@all_loc, $scale);

    # Make empty image to draw on
    my ($img, $col) = $self->_create_blank_background_image(\@all_loc, scalar(@$sec_loc) ? 1 : 0);
    my @conn_color = ($col->{conn1}, $col->{conn_inv1});

    # Draw everything
    my $y = 0;
    if($ref_index == 1 and @$sec_loc) {
	my $mid_loc = shift @$sec_loc;
	$sec_loc = [reverse @$sec_loc];
	$y = $self->_draw_backgrounds_and_add_panels($img, $sec_loc, $y, $scale, [$col->{bg1}]);
	$y = $self->_draw_backgrounds_and_add_panels($img, [$mid_loc], $y, $scale, [$col->{bg1}]);
	$y = $self->_draw_relations($img, @conn_color, $mid_loc, $ref_loc, 0, 0, $y, $scale);
    }
    $y = $self->_draw_backgrounds_and_add_panels($img, [$ref_loc], $y, $scale, [$col->{bg1}, $col->{bg2}]);
    if($ref_index == 0 and @$sec_loc) {
	my $mid_loc = shift @$sec_loc;
	$y = $self->_draw_relations($img, @conn_color, $ref_loc, $mid_loc, 0, 0, $y, $scale);
	$y = $self->_draw_backgrounds_and_add_panels($img, [$mid_loc], $y, $scale, [$col->{bg1}]);
	$y = $self->_draw_backgrounds_and_add_panels($img, $sec_loc, $y, $scale, [$col->{bg2}]);
    }

    $self->{_gd_image} = $img;
    return 1;
}


sub draw_outward_triplet
{
    my ($self, $coord_sets) = @_;

    my $rgn_width = $self->{width};
    my $conn_height = $self->{connector_height};

    # Get specs for the reference assembly
    my ($ref_index, $chr, $rgn_start, $rgn_end) = @{$coord_sets->[0]};
    my $ref_loc = MyLoc->new(asm_param => $self->asm_param->[$ref_index],
			     chr => $chr, start => $rgn_start, end => $rgn_end, strand => '+',
			     offset => 0);
    
    # Check that index value is valid
    unless($ref_index == 1) {
	$self->{_errstr} = "Invalid reference index $ref_index\n";
	return 0;
    }

    # Set scale to the highest among the refs
    my $scale = ($rgn_end - $rgn_start + 1) / $rgn_width;
    print STDERR "using scale: $scale bp/pixel\n";

    # Get ranked alignment groups, aligned ranges and offsets for each reference
    my ($aln_groups1, $sec_loc1) = $self->_get_aln_data(0, $ref_loc);
    my ($aln_groups2, $sec_loc2) = $self->_get_aln_data(2, $ref_loc);

    # Increment offsets so the lowest is zero, and scale the offsets
    my @all_loc = ($ref_loc, @$sec_loc1, @$sec_loc2);
    $self->_adjust_offsets($scale, \@all_loc);

    # Add context secondary panels if there is space
    if($self->{fill_with_context}) {
	$self->_fill_with_context($ref_loc, $sec_loc1, $scale);
	$self->_fill_with_context($ref_loc, $sec_loc2, $scale);
    }

    # Draw panels
    $self->_draw_panels(\@all_loc, $scale);

    # Make empty image to draw on
    my $nr_connectors = (scalar(@$sec_loc1) ? 1 : 0) + (scalar(@$sec_loc2) ? 1 : 0);
    my ($img,$col) = $self->_create_blank_background_image(\@all_loc, $nr_connectors);
    my @conn_color = ($col->{conn1}, $col->{conn_inv1});
    my $y = 0;

    my @bg_color_list;
    my $nr_upper_panels = @$sec_loc1;
    # Add first secondary panels 
    if(@$sec_loc1) {
	my $mid_loc = shift @$sec_loc1;
	$sec_loc1 = [reverse @$sec_loc1];
	$y = $self->_draw_backgrounds_and_add_panels($img, $sec_loc1, $y, $scale, [$col->{bg2}]);
	$y = $self->_draw_backgrounds_and_add_panels($img, [$mid_loc], $y, $scale, [$col->{bg1}]);
	$y = $self->_draw_relations($img, @conn_color, $mid_loc, $ref_loc, 0, 0, $y, $scale);
	@bg_color_list = ($col->{bg1}, ($col->{bg2}) x scalar(@$sec_loc1), $col->{bg1}, $col->{bg2});
    }
    else {
	@bg_color_list = ($col->{bg1}, $col->{bg2});
    }
    # Add ref
    $y = $self->_draw_backgrounds_and_add_panels($img, [$ref_loc], $y, $scale, \@bg_color_list);
    if(@$sec_loc2) {
	my $mid_loc = shift @$sec_loc2;
	$y = $self->_draw_relations($img, @conn_color, $ref_loc, $mid_loc, $nr_upper_panels, 0, $y, $scale);
	$y = $self->_draw_backgrounds_and_add_panels($img, [$mid_loc], $y, $scale, [$col->{bg1}]);
	$y = $self->_draw_backgrounds_and_add_panels($img, $sec_loc2, $y, $scale, [$col->{bg2}]);
    }

    $self->{_gd_image} = $img;
    return 1;
}


sub draw_inward_triplet
{
    my ($self, $coord_sets) = @_;

    my $rgn_width = $self->{width};
    my $conn_height = $self->{connector_height};

    # Get specs for the two reference assemblies
    my ($ref_index1, $chr1, $rgn_start1, $rgn_end1) = @{$coord_sets->[0]};
    my ($ref_index2, $chr2, $rgn_start2, $rgn_end2) = @{$coord_sets->[1]};
    my $ref_loc1 = MyLoc->new(asm_param => $self->asm_param->[$ref_index1],
			      chr => $chr1, start => $rgn_start1, end => $rgn_end1, strand => '+',
			      offset => 0);
    my $ref_loc2 = MyLoc->new(asm_param => $self->asm_param->[$ref_index2],
			      chr => $chr2, start => $rgn_start2, end => $rgn_end2, strand => '+',
			      offset => 0);
    
    # Check that index values are valid
    unless($ref_index1 != $ref_index2 and $ref_index1 >= 0 and $ref_index1 <= 2 and $ref_index2 >= 0 and $ref_index2 <= 2) {
	croak "Invalid reference indices $ref_index1, $ref_index2\n";
    }

    # Determine which assembly is to be shown in the middle (the non-reference assembly)
    my $mid_index = 3 - $ref_index1 - $ref_index2;

    # Set scale to the highest among the refs
    my $span1 = $rgn_end1 - $rgn_start1 + 1;
    my $span2 = $rgn_end2 - $rgn_start2 + 1;
    my $scale = ($span1 > $span2 ? $span1 : $span2) / $rgn_width;
    print STDERR "using scale: $scale bp/pixel\n";

    # Get ranked alignment groups, aligned ranges and offsets for each reference
    my ($aln_groups1, $sec_loc1) = $self->_get_aln_data($mid_index, $ref_loc1);
    my ($aln_groups2, $sec_loc2) = $self->_get_aln_data($mid_index, $ref_loc2);

    my $mid_loc; 
  
    if(@$aln_groups1 and @$aln_groups2) {
	# Get the best target region for the top locus, and find an overlapping target region for
	# the bottom locus
	my $mid_loc1 = $sec_loc1->[0];
	my $mid_loc2;
	for my $i (0..@$sec_loc2-1) {
	    my $l = $sec_loc2->[$i];
	    if($mid_loc1->chr eq $l->chr and
	       $mid_loc1->start <= $l->end and $mid_loc1->end >= $l->start) {
		# Remove the middle region for the list of target regions for the first locus
		shift @$sec_loc1;
		# Remove the middle region from the list of target regions for the second locus
		$mid_loc2 = $l;
		splice(@$sec_loc2, $i, 1);
		# Calculate combined mid bounds ( min(s1,s2), max(e1,e2) ) 
		$mid_loc = $mid_loc1;
		if($mid_loc->start > $mid_loc2->start) {
		    $mid_loc->start($mid_loc2->start);
		    $mid_loc->offset($mid_loc->offset - $mid_loc2->start + $mid_loc->start);
		}
		$mid_loc->end($mid_loc2->end) if($mid_loc->end < $mid_loc2->end);
		# Fix aligned ranges
		$mid_loc->aln_regions->[1] = $mid_loc2->aln_regions->[0];
		my $ref_loc2_aln_ranges = $ref_loc2->aln_regions;
		unshift @$ref_loc2_aln_ranges, splice(@$ref_loc2_aln_ranges, $i, 1);
		# If strand is not the same, we need to flip all regions and offsets for the second group
		#  i.e. invert all strands (index 3) and do something :-) with the offsets
		if($mid_loc->strand ne $mid_loc2->strand) {
		    foreach my $m ($ref_loc2, @$sec_loc2) {
			$m->strand($m->strand eq '+' ? '-' : '+');
			$m->offset(($mid_loc2->end-$mid_loc2->start+1) - ($m->end-$m->start+1) - $m->offset + $mid_loc2->offset);
		    }
		    $mid_loc2->offset(0);
		}
		# Reconciliate offsets (based on mid region, start bound)
		#   d = s2 - s  + o - o2
		#   add d to all offsets from bottom panel 
		my $d;
		if($mid_loc->strand eq '+') {
		    $d = $mid_loc2->start - $mid_loc->start + $mid_loc->offset - $mid_loc2->offset;
		}
		else {
		    $d = $mid_loc->end - $mid_loc2->end + $mid_loc->offset - $mid_loc2->offset;
		}
		foreach my $m ($ref_loc2, @$sec_loc2) {
		    $m->offset($m->offset + $d);
		}
		last;
	    }
	}
    }

    # Increment offsets so the lowest is zero, and scale the offsets
    my @all_loc = ($ref_loc1, $ref_loc2, $mid_loc?$mid_loc:(), @$sec_loc1, @$sec_loc2);
    $self->_adjust_offsets($scale, \@all_loc);

    # Add context secondary panels if there is space
    if($self->{fill_with_context}) {
	$self->_fill_with_context($ref_loc1, $sec_loc1, $scale);
	$self->_fill_with_context($ref_loc2, $sec_loc2, $scale);
	if($mid_loc) {
	    $self->_fill_with_context($ref_loc1, [$mid_loc], $scale);
	    $self->_fill_with_context($ref_loc2, [$mid_loc], $scale);
	}
    }

    # Draw panels
    $self->_draw_panels(\@all_loc, $scale);

    # Unless there is a shared middle locus, get two middle loci to draw
    my ($mid_loc1, $mid_loc2, $nr_connectors);
    if($mid_loc) {
	$nr_connectors = 2 + (scalar(@$sec_loc1)?1:0) + (scalar(@$sec_loc2)?1:0);
    }
    else {
	$mid_loc1 = shift @$sec_loc1;
	$mid_loc2 = shift @$sec_loc2;
	$nr_connectors = ($mid_loc1?1:0) + ($mid_loc2?1:0) + (scalar(@$sec_loc1)?1:0) + (scalar(@$sec_loc2)?1:0);
    }

    # Make empty image to draw on
    my ($img, $col) = $self->_create_blank_background_image(\@all_loc, $nr_connectors);
    my @conn_color1 = ($col->{conn1}, $col->{conn_inv1});
    my @conn_color2 = ($col->{conn2}, $col->{conn_inv2});
    my $y = 0;

    # Add secondary panels for ref1 
    if(@$sec_loc1) {
	$sec_loc1 = [reverse @$sec_loc1];
	$y = $self->_draw_backgrounds_and_add_panels($img, $sec_loc1, $y, $scale, [$col->{bg2}]);
	$y = $self->_draw_relations($img, @conn_color2, $sec_loc1->[-1], $ref_loc1, 0, 1, $y, $scale);
    }
    # Add ref1
    $y = $self->_draw_backgrounds_and_add_panels($img, [$ref_loc1], $y, $scale, [$col->{bg1}, $col->{bg2}]);
    # Add middle
    if($mid_loc) {
	$y = $self->_draw_relations($img, @conn_color1, $ref_loc1, $mid_loc, 0, 0, $y, $scale);
	$y = $self->_draw_backgrounds_and_add_panels($img, [$mid_loc], $y, $scale, [$col->{bg1}]);
	$y = $self->_draw_relations($img, @conn_color1, $mid_loc, $ref_loc2, 1, 0, $y, $scale);
    }
    else {
	if($mid_loc1) {
	    $y = $self->_draw_relations($img, @conn_color1, $ref_loc1, $mid_loc1, 0, 0, $y, $scale);
	    $y = $self->_draw_backgrounds_and_add_panels($img, [$mid_loc1], $y, $scale, [$col->{bg1}]);
	}
	if($mid_loc2) {
	    $y = $self->_draw_backgrounds_and_add_panels($img, [$mid_loc2], $y, $scale, [$col->{bg1}]);
	    $y = $self->_draw_relations($img, @conn_color1, $mid_loc2, $ref_loc2, 0, 0, $y, $scale);
	}
    }
    # Add ref2
    $y = $self->_draw_backgrounds_and_add_panels($img, [$ref_loc2], $y, $scale, [$col->{bg1}, $col->{bg2}]);
    # Add secondary panels for ref2
    if(@$sec_loc2) {
	$y = $self->_draw_relations($img, @conn_color2, $ref_loc2, $sec_loc2->[0], 1, 0, $y, $scale);
	$y = $self->_draw_backgrounds_and_add_panels($img, $sec_loc2, $y, $scale, [$col->{bg2}]);
    }

    $self->{_gd_image} = $img;
    return 1;
}


sub _create_blank_background_image
{
    my ($self, $loc_list, $nr_connectors) = @_;
    my ($width, $height) = (0,0);
    foreach my $loc (@$loc_list) {
	my ($w,$h) = $loc->image->getBounds();
	$w += $loc->offset;
	$width = $w if($w > $width);
	$height += $h;
    }
    $height += $nr_connectors * ($self->{connector_height} + 10);
    $height += 2; # this should not be needed, but just in case
    my $img = $self->{_image_pkg}->new($width,$height);

    my $colors = {canvas => $img->colorAllocate(255,255,255),
		  bg1 => $img->colorAllocate(220,220,220),
		  bg2 => $img->colorAllocate(200,255,200),
		  bg3 => $img->colorAllocate(200,200,255),
		  conn1 => $img->colorAllocate(0,0,0),
		  conn2 => $img->colorAllocate(0,180,0),
		  conn_inv1 => $img->colorAllocate(180,0,0),
		  conn_inv2 => $img->colorAllocate(180,180,0)};

    return ($img, $colors);
}


sub _get_aln_data
{
    my ($self, $sec_index, $ref_loc) = @_;

    my $sec_asm_param = $self->asm_param->[$sec_index];
    my $alndb = $ref_loc->asm_param->{alignment_db};
    my $gendb1 = $ref_loc->asm_param->{assembly_db};
    my $gendb2 = $sec_asm_param->{assembly_db};
    my $max_nr_sec_panels = $self->max_nr_secondary_panels;
    my $context = $self->min_context;
    my $t0; # For benchmarking

    # Get alignments over region;
    $t0 = Benchmark->new;
    my $alnset = $alndb->get_alignments_for_region(chr => $ref_loc->chr,
						   start => $ref_loc->start,
						   end => $ref_loc->end,
						   seq_dbs => [$gendb1, $gendb2],
						   aln_type => 'net',
						   confine => 1);
    my $alns = $alnset->get_alignments || [];
# HERE: if just one region specified
    my $aln_groups = get_alignment_groups($alnset,$self->max_space_factor*($ref_loc->end-$ref_loc->start));
    $aln_groups = sort_aln_groups_by_span($aln_groups);
# HERE: else: find all alignments in 2nd region, first then make alignment groups and sort the rest
    print STDERR "Got ".scalar(@$aln_groups)." alignment groups ".timediff_str($t0,Benchmark->new)."\n";

    # Limit number of secondary panels if requested
    if($max_nr_sec_panels >= 0) {
	$#$aln_groups = scalar(@$aln_groups) < $max_nr_sec_panels ? @$aln_groups-1 : $max_nr_sec_panels-1;
    }

    # Get coordinate ranges:
    #   aligned regions in the primary panel
    #   x offsets for secondary and primary panels
    $t0 = Benchmark->new;
    my @sec_loc;
    foreach my $alng (@$aln_groups) {
        my (undef, $chr2, $aln_strand, $r1, $r2) = get_aligned_ranges($alng, $gendb1->id, $gendb2->id);
	# Find bounds of sets of aligned ranges
        my ($start1, $end1); # For primary panel
	my ($start2, $end2); # For secondary panel
        foreach my $r (@$r1) {
	    $start1 = $r->[0] if(!$start1 or $start1 > $r->[0]);
	    $end1 = $r->[1] if(!$end1 or $end1 < $r->[1]);
	}
	foreach my $r (@$r2) {
	    $start2 = $r->[0] if(!$start2 or $start2 > $r->[0]);
	    $end2 = $r->[1] if(!$end2 or $end2 < $r->[1]);
	}
	$start2 = $start2 - $context;
	$end2 = $end2 + $context;
	# Determine offset for secondary panel
	my $o = ((($end1 - $start1 + 1) - ($end2 - $start2 + 1)) / 2 + $start1-$ref_loc->start);
	# Create new loc object for secondary panel
	my $sec_loc = MyLoc->new(asm_param => $sec_asm_param,
				 chr => $chr2, start => $start2, end => $end2, strand => $aln_strand,
				 aln_regions => [$r2],
				 offset => $o);
	push @sec_loc, $sec_loc;
	# Record the aligned ranges for the reference panel
	push @{$ref_loc->aln_regions}, $r1;
    }
    print STDERR "got coord ranges ".timediff_str($t0, Benchmark->new)."\n";

    return ($aln_groups, \@sec_loc);
}    


sub _adjust_offsets
{
    my ($self, $scale, $loc_list) = @_;
    # find lowest offset
    my $min = $loc_list->[0]->offset;
    foreach my $loc (@$loc_list) {
	$min = $loc->offset if($min > $loc->offset);
    }
    # adjust all offsets
    foreach my $loc (@$loc_list) {
	$loc->offset(int(($loc->offset - $min) / $scale))
    }
}


sub _fill_with_context
{
    my ($self, $ref_loc, $sec_loc_list, $scale) = @_;

    # expand each secondary loc to the left so it fills down to the offset for the reference panel
    # and to the right so it fills up to the end of the reference panel

    my $left_bound = $ref_loc->offset;
    my $right_bound = $ref_loc->offset + ($ref_loc->end - $ref_loc->start + 1) / $scale;

    foreach my $sec_loc (@$sec_loc_list) {
	my $chr_size = $sec_loc->asm_param->{assembly_db}->get_chr_size($sec_loc->chr);
	if($sec_loc->offset > $left_bound) {
	    if($sec_loc->strand eq '+') {
		my $new_start = $sec_loc->start - int(($sec_loc->offset - $left_bound) * $scale + .5);
		if($new_start >= 1) {
		    $sec_loc->offset($left_bound);
		    $sec_loc->start($new_start);
		}
		else {
		    $sec_loc->offset(int($sec_loc->offset - $sec_loc->start / $scale + .5));
		    $sec_loc->start(1);
		}
	    }
	    else {
		my $new_end = $sec_loc->end + int(($sec_loc->offset - $left_bound) * $scale + .5);
		if($new_end <= $chr_size) {
		    $sec_loc->offset($left_bound);
		    $sec_loc->end($new_end);
		}
		else {
		    $sec_loc->offset(int($sec_loc->offset - ($chr_size - $sec_loc->end) / $scale + .5));
		    $sec_loc->end($chr_size);
		}
	    }
	}
	my $right_end = $sec_loc->offset + ($sec_loc->end - $sec_loc->start + 1) / $scale;
	if($right_end < $right_bound) {
	    if($sec_loc->strand eq '+') {
		my $new_end = $sec_loc->end + int(($right_bound - $right_end) * $scale + .5);
		$new_end = $chr_size if($new_end > $chr_size);
		$sec_loc->end($new_end);
	    }
	    else {
		my $new_start = $sec_loc->start - int(($right_bound - $right_end) * $scale + .5);
		$new_start = 1 if($new_start < 1);
		$sec_loc->start($new_start);
	    }
	}

    }

}


sub timediff_str
{
    my ($b1, $b2) = @_;
    my $t = timestr(timediff($b2, $b1));
    return "[$t sec]";
}

sub gd_image
{
    return shift->{_gd_image};
}



sub get_aligned_ranges
{
    my ($alns, $prefix1, $prefix2) = @_;
    my (@r1, @r2);
    my ($chr1, $chr2);
    my $plus_strand_span = 0;
    my $minus_strand_span = 0;
    foreach my $caln (@$alns) {
        my $seq1 = $caln->get_seq1;
        my $seq2 = $caln->get_seq2;
	my $strand = $caln->get_strand;
	unless($chr1) {
	    $chr1 = $prefix1 ? substr($seq1->id, length($prefix1)+1) : $seq1->id;
	    $chr2 = $prefix2 ? substr($seq2->id, length($prefix2)+1) : $seq2->id;
	    #($chr1) = ($seq1->id =~ /(chr.+)/);
	    #($chr2) = ($seq2->id =~ /(chr.+)/);
	}
        foreach my $a (@{$caln->get_subalignments}) {
	    # get absolute coords for subalignment
            my ($rel_start1, $rel_end1) = @{$a->get_pos_seq1};
            my ($rel_start2, $rel_end2) = @{$a->get_pos_seq2};
            my $abs_start1 = $seq1->start + $rel_start1 - 1;
            my $abs_end1 = $seq1->start + $rel_end1 - 1;
            my ($abs_start2, $abs_end2);
            if($strand eq '+') {
        	$abs_start2 = $seq2->start + $rel_start2 - 1;
        	$abs_end2 = $seq2->start + $rel_end2 - 1;
		$plus_strand_span += $abs_end2 - $abs_start2 + 1;
            }
	    else {
	        $abs_start2 = $seq2->end - $rel_end2 + 1;
	        $abs_end2 = $seq2->end - $rel_start2 + 1;
		$minus_strand_span += $abs_end2 - $abs_start2 + 1;
	    }
	    push @r1, [$abs_start1, $abs_end1, 1];
	    push @r2, [$abs_start2, $abs_end2, $strand eq '+' ? 1 : 0];
	}
    }
    my $strand;
    if($plus_strand_span >= $minus_strand_span) {
	$strand = '+';
    }
    else {
	$strand = '-';
	foreach my $r (@r2) {
	    $r->[2] = !($r->[2]);
	}
    }
    #my $strand = $plus_strand_span >= $minus_strand_span ? '+' : '-';
    return ($chr1, $chr2, $strand, \@r1, \@r2);
}


sub get_alignment_groups
{
    my ($alnset, $max_d, $sep_strands) = @_;
    my $alns_ref = $alnset->get_alignments || [];
    my @alns;
    if($sep_strands) {
	@alns = sort { $a->get_seq2->id cmp $b->get_seq2->id or
			   $a->get_strand cmp $b->get_strand or
			   $a->get_seq2->end <=> $b->get_seq2->end } @$alns_ref;
    }
    else {
	@alns = sort { $a->get_seq2->id cmp $b->get_seq2->id or
			   $a->get_seq2->end <=> $b->get_seq2->end } @$alns_ref;
    }

    my @groups;      
    my (@group, $group_strand, $group_end);
    my $group_id = "";
    foreach my $aln (@alns) {
	my $seq = $aln->get_seq2;
	if($seq->id eq $group_id and
	   (!$sep_strands or $aln->get_strand eq $group_strand)
	   and $seq->start < $group_end+$max_d) {
	    push @group, $aln;
	    $group_end = $seq->end if($group_end < $seq->end);
	}
	else {
	    push @groups, [@group] if(@group);
	    $group_id = $seq->id;
	    $group_strand = $aln->get_strand;
	    $group_end = $seq->end;
	    @group = ($aln);
        }
    }
    push @groups, [@group] if(@group);

#    foreach my $g (@groups) {
#	print STDERR "Group\n";
#	foreach my $aln (@$g) {
#	    my $seq2 = $aln->get_seq2;
#            print STDERR $seq2->id, " ", $aln->get_strand, " ", $seq2->start, "-", $seq2->end, "\n";
#	}
#    }
    
    return \@groups;
}


sub sort_aln_groups_by_span
{
    my $groups = shift;

    my %span;    
    foreach my $group (@$groups) {
	my $s = 0;
	foreach my $aln (@$group) {
	    foreach my $a (@{$aln->get_subalignments}) {
		my ($rel_start1, $rel_end1) = @{$a->get_pos_seq1};    
		$s += ($rel_end1 - $rel_start1 + 1);
	    }
	}
	$span{$group} = $s;
    }

    my @sorted_groups = sort { $span{$b} <=> $span{$a} } @$groups;

    return \@sorted_groups;
}


sub get_largest_aln
{
    my $alnset = shift;
    my $alns = $alnset->get_alignments || return;
    my $largest_aln;
    my $max_span = 0;
    foreach my $aln (@$alns) {
	my $span = 0;
	foreach my $a (@{$aln->get_subalignments}) {
            my ($rel_start1, $rel_end1) = @{$a->get_pos_seq1};    
	    $span += ($rel_end1 - $rel_start1 + 1);
	}
	if($span > $max_span) {
	    $max_span = $span;
	    $largest_aln = $aln;
	}
    }
    return $largest_aln;
}


sub print_aln
{
    my ($a, $stream) = @_;
    $stream = \*STDOUT unless (defined $stream);
    my $alnout = Bio::AlignIO->new(-format => 'clustalw', -fh => $stream);

    print $stream "Genomic alignment  [",$a->get_strand,' ',$a->get_score,"]\n";
    my @sa = @{$a->get_subalignments};
    foreach my $aln (@sa) {
	my $alnobj = $aln->get_alignment_obj();
	$alnout->write_aln($alnobj);
    }
}


sub loc_str
{
    my ($chr,$start,$end,$strand) = @_;
    if($strand eq '1') { $strand = '+' }
    elsif($strand eq '-1') { $strand = '-' }
    return "$chr:$start-$end:$strand";
}


#####################
# GFX routines
#####################


sub _draw_backgrounds_and_add_panels
{
    my ($self, $img, $loc_list, $y, $scale, $bg_colors) = @_;

    my $l_margin = $self->{pad_left};
    my $t_margin = $self->{pad_top};
    my $max_color_idx = @$bg_colors-1;

    foreach my $loc (@$loc_list) {

	# Draw background at aligned segments
	my $sub_img = $loc->image;
	my $aln_rgn_sets = $loc->aln_regions;
	my ($w,$h) = $sub_img->getBounds();
	my $color_idx = 0;
	foreach my $regions (@$aln_rgn_sets) {
	    my $c = $bg_colors->[$color_idx];
	    $self->_draw_aln_regions($img, $c,
				     $l_margin+$loc->offset, $y+$t_margin+35, $h-$t_margin-35,
				     $loc->start, $loc->end, $loc->strand,
				     $scale, $regions);
	    $color_idx++ if($color_idx < $max_color_idx);
	}

	# Copy panel
        #$sub_img->transparent($sub_img->colorClosest(255,255,255));
	$sub_img->transparent($sub_img->getPixel(0,0));
        $img->copy($sub_img, $loc->offset, $y, 0, 0, $w, $h);

	$y += $h;
    }

    return $y;
}


sub _draw_aln_regions
{
    my ($self, $img, $color, $x_offset, $y_offset, $h, $start, $end, $strand, $scale, $regions) = @_;
    foreach my $r (@$regions) {
	my ($s, $e) = @$r;
	if($strand eq '+') {
	    $s = ($s - $start ) / $scale + $x_offset;
	    $e = ($e - $start ) / $scale + $x_offset;
	}
	else {
	    ($s,$e) = (($end - $e) / $scale + $x_offset,
		       ($end - $s) / $scale + $x_offset);
	}
	$s = int($s+.5);
	$e = int($e+.5);
        $img->filledRectangle($s, $y_offset,
			      $e, $y_offset+$h-1,
			      $color);
    }
}


sub _draw_relations
{
    my ($self, $img, $def_color, $inv_color, $loc1, $loc2, $aln_idx1, $aln_idx2, $y, $scale) = @_;

    my $r1 = $loc1->aln_regions->[$aln_idx1]; 
    my $r2 = $loc2->aln_regions->[$aln_idx2]; 
    my $pad_left1 = $self->{pad_left} + $loc1->offset;
    my $pad_left2 = $self->{pad_left} + $loc2->offset;
    my $start1 = $loc1->start;
    my $start2 = $loc2->start;
    my $end1 = $loc1->end;
    my $end2 = $loc2->end;
    my $strand1 = $loc1->strand;
    my $strand2 = $loc2->strand;
    my $flip1 = $strand1 eq '+' ? 0 : 1;
    my $flip2 = $strand2 eq '+' ? 0 : 1;
    my $height = $self->{connector_height};
    $height = 15 if($height < 15);
    $y += 5;
    
    for my $i (0..@$r1-1) {
	my ($s1,$e1,$straight1) = @{$r1->[$i]};
	my ($s2,$e2,$straight2) = @{$r2->[$i]};
	my $color = ($straight1 and $straight2) ? $def_color : $inv_color;
        #$s1 = ($s1 - $start1) / $scale + $pad_left1;
        #$e1 = ($e1 - $start1) / $scale + $pad_left1;
	if($flip1) {
	    ($s1,$e1) = (($end1 - $e1) / $scale + $pad_left1,
			 ($end1 - $s1) / $scale + $pad_left1);
	}
	else {
	    $s1 = ($s1 - $start1) / $scale + $pad_left1;
	    $e1 = ($e1 - $start1) / $scale + $pad_left1;
	}
	if($flip2) {
	    ($s2,$e2) = (($end2 - $e2) / $scale + $pad_left2,
			 ($end2 - $s2) / $scale + $pad_left2);
	}
	else {
	    $s2 = ($s2 - $start2) / $scale + $pad_left2;
	    $e2 = ($e2 - $start2) / $scale + $pad_left2;
	}
	my $p = $self->{_polygon_pkg}->new();
	$s1 = int($s1+.5);
	$e1 = int($e1+.5);
	$s2 = int($s2+.5);
	$e2 = int($e2+.5);
	$p->addPt($s1, $y);
	$p->addPt($s1, $y+5);
	$p->addPt($s2, $y+$height-6);
	$p->addPt($s2, $y+$height-1);
	$p->addPt($e2, $y+$height-1);
	$p->addPt($e2, $y+$height-6);
	$p->addPt($e1, $y+5);
	$p->addPt($e1, $y);
	$img->filledPolygon($p,$color);
	#print STDERR "$loc1 $loc2 $s1-$e1 $s2-$e2\n";
    }

    return $y + $height + 5;
}


sub _draw_panels
{
    my ($self, $loc_list, $scale) = @_;
    my $t0 = Benchmark->new;
    foreach my $loc (@$loc_list) {
	$loc->image($self->_gen_panel($loc, $scale)->gd);
    }
    print STDERR "generated panels ".timediff_str($t0,Benchmark->new)."\n";
}


sub _gen_panel
{
    my ($self, $loc, $scale) = @_;

    my $params = $loc->asm_param;
    croak "No parameters for assembly" unless($params);
    my $width = int(($loc->end-$loc->start+1)/$scale);
    $width = 1 if($width < 1);

    my $drawer = AT::GFX::Locus->new(%$params,
				     width => $width, #+$self->{pad_left}+$self->{pad_right},
				     pad_left => $self->{pad_left},
				     pad_right => $self->{pad_right},
				     pad_top => $self->{pad_top},
				     pad_bottom => $self->{pad_bottom},
				     flip => $loc->strand eq '+' ? 0 : 1,
				     svg => $self->{svg});

    $drawer->draw_locus($loc->chr,$loc->start,$loc->end);

    return $drawer->panel;
}



1;
