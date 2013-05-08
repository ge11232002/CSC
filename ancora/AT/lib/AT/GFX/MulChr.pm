# AT::GFX::MulChr module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::GFX::MulChr - draw multiple chromosomes with features

=head1 SYNOPSIS

 use AT::GFX::MulChr;
 use AT::DB::GenomeAssemblyNibs;

 my $gendb = AT::DB::GenomeAssemblyNibs->new(assembly_name => 'HS_APR03',
        dir => '/net/mordor/data/goldenpath/HS_APR03/nib');

 # read in positions of features for first track
 my @track1_feats;
 open POINTS, "points.gff" || die;
 while (my $line = <POINTS>) {
    chomp $line;
    my ($chr, undef, undef, $start, $end) = (split /\t/, $line);
    $chr =~ s/^chr//;
    push @track1_feats, {
        chr => $chr,
        start => $start,
        end => $end,
        color => 'red'};
 }
 close POINTS;

 # read in positions of features for second track
 # (we'll make this a histogram track)
 my @track2_feats;
 open DENS, "densities.gff" || die;
 while (my $line = <DENS>) {
    chomp $line;
    my ($chr, undef, undef, $start, $end, $score) = split /\t/, $line;
    $chr =~ s/^chr//;
    push @track2_feats, {
        chr => $chr,
        start => $start,
        end => $end,
        score => $score }; # histogram features have a score
 }
 close DENS;

 # Create the tracks
 my $track1 = {
    side => 0,
    features => \@track1_feats
 };
 my $track2 = {
    side => 1,
    features => \@track2_feats,
    type => 'histogram',
    color => 'green',
    feat_width => 50, # height (in pixels) of a bar with score max_score
    max_score => 100000,
    min_score => 0
 };

 # Draw the image and display it
 my $img = AT::GFX::MulChr->new(x_step => 100);
 $img->draw($gendb, [(1..22, 'X', 'Y')], [$track1, $track2] );
 $img->display;

=head1 DESCRIPTION

Draw a set of chromosomes with features. Features may be points or histogram
bars. Histogram bars are useful for displaying e.g. gene density data.
The module is highly configurable through the various mysterious parameters.
Any number of tracks can be added. Each track can be on the left or right
side of chromosomes.

Point tracks may contain any number of features and the
features can be colored differently. Point features in the same track are
bumped if they collide.

Histogram tracks must contain non-overlapping features and they will all
have the same color, specified in the track definition. There can be weird
effects if they overlap even by 1 bp (yes, we should implement a check for this).
In addition, if a very tall image is made, gaps can appear between adjacent bars
(due to rounding errors?).
The image can be drawn directly on screen, as a png bitmap or a svg vector format.
Alternatively, png:s and svg:s cant be returned as strings (binary resp. XML) for piping or redirection.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::GFX::MulChr;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use AT::GFX::ColorHash;
#use GD;
#use GD::SVG;
use Data::Dumper;

@ISA = qw(AT::Root);


=head2 new

 Title     : new
 Usage     :
 Function  : Constructor
 Returns   : AT::GFX::MulChr
 Args      : ARGUMENT			DEFAULT VALUE
	     height			2000
	     left_margin		15
	     right_margin		15
	     top_margin			10
	     bottom_margin		10
	     label_height		20
	     x_step			38
	     chr_width			7
	     chr_feat_spacer		2
	     feat_width			2
	     feat_spacer		2
	     feat_min_span		2
	     feat_bump_threshold	2
             track_spacer               2
	     svg			undef
	    Some explanations: 
	    <height> is the height of the image.
	    <x_step> is the distance between chromosomes
	    <chr_feat_spacer> is the distance between
	     chromosomes and the nearest features
	    <feat_spacer> is the x distance between features
	    <feat_min_span> is the minimum extent of a feature
	    <feat_bump_threshold> is the minimum y distance between
	     features; set to -1 to disable bumping (collapse
	     features)
            <track_spacer> ...
	    if <svg> is true, image is drawn as vectors, not as bitmap (see below)
	    For all arguments except <svg>, the unit is pixels.
	    Some arguments only apply to point features.

=cut

sub new
{
    my ($caller, %args) = @_;

    my $svg = $args{svg};
    my $gd_pkg = $svg ? 'GD::SVG' : 'GD';
    eval "use $gd_pkg";

    my $self = bless {
	height => $args{height} || 2000,
	left_margin => $args{left_margin} || 15,
	right_margin => $args{right_margin} || 15,
	top_margin => $args{top_margin} || 10,
	bottom_margin => $args{bottom_margin} || 10,
	label_height => $args{label_height} || 20,
	x_step => $args{x_step} || 38,
	chr_width => $args{chr_width} || 7,
	chr_feat_spacer => $args{chr_feat_spacer} || 2,
	feat_width => $args{feat_width} || 2,
	feat_spacer => $args{feat_spacer} || 2,
	track_spacer => $args{track_spacer} || 2,
	feat_min_span => $args{feat_min_span} || 2,
	feat_bump_threshold => $args{feat_bump_threshold} || 2,
	color_xlat => AT::GFX::ColorHash->new,
	chr_pattern => $args{chr_pattern} || [],
	chr_thin => $args{chr_thin} || [],
	default_chr_color => $args{default_chr_color} || 'black',
	colors => {},
	svg => $svg,
	_image_pkg => $gd_pkg.'::Image',
	_polygon_pkg => $gd_pkg.'::Polygon',
	_font_pkg => $gd_pkg.'::Font'
    }, ref $caller || $caller;

    return $self;
}


sub _new_img {
    my ($self, $num_chr) = @_;

    # calculate required width
    $self->{width} = $self->left_margin + $self->right_margin +
		     $self->x_step * ($num_chr+1);

    # create image object
    $self->{img} = $self->{_image_pkg}->new($self->width, $self->height);
    
    # allocate some colors
    for my $col_name ('white','black', $self->color_xlat->names) {
	$self->_add_color($col_name,
			 @{$self->color_xlat->rgb_from_name($col_name)});
    }
} 


sub _add_color {
    my ($self, $name, @rgb) = @_;
    $self->colors->{$name} = $self->img->colorAllocate(@rgb);
}


=head2 draw

 Title     : draw
 Usage     : See sypnosis
 Function  : Create and draw a new image
 Returns   : -
 Args      : See sypnosis
 
=cut

sub draw {
    my ($self, $gendb, $chr_list, $tracks) = @_;

    $self->_new_img(scalar @$chr_list);
    $self->_draw($gendb, $chr_list, $tracks);
}


sub _draw {
    my ($self, $gendb, $chr_list, $tracks) = @_;

    my $img = $self->{img};

    my $y_offset = $self->top_margin + $self->label_height;
    my $x_offset = $self->left_margin;
    my $x_step = $self->x_step;
    my $max_chrlen_pix = $self->height - $y_offset - $self->bottom_margin;

    my $max_chrlen_bp = 0;
    my %chr_data;
    for my $i (0..@$chr_list-1) {
	my %data;
	my $name = $chr_list->[$i];
	$data{length} = $gendb->get_chr_size($name);
	print STDERR "$name\t",$data{length},"\n";
	$data{column} = $i;		# Column the chromosome is drawn in
	$data{fcol_maxends} = [];
	$data{fcol_bases} = [0,0];
	$chr_data{$name} = \%data;
	$max_chrlen_bp = $data{length} if($max_chrlen_bp < $data{length});
   }
   foreach my $chr_region (@{$self->{chr_pattern}}) {
	push @{$chr_data{$chr_region->{chr}}{pattern}}, $chr_region;
   }
   foreach my $chr_region (@{$self->{chr_thin}}) {
	push @{$chr_data{$chr_region->{chr}}{thin}}, $chr_region;
   }

    my $pix_per_bp = $max_chrlen_pix / $max_chrlen_bp;
    #print STDERR "using $pix_per_bp pixels/basepair\n";

    # Draw chromosomes
    my $x = $x_offset + $x_step;
    for my $i (0..@$chr_list-1) {
	$img->stringUp($self->{_font_pkg}->MediumBold, $x-7, $y_offset - 7, $chr_list->[$i],
		     $self->{colors}->{black});
	my $cd = $chr_data{$chr_list->[$i]};
	$self->_draw_chromosome($x, $y_offset, $cd->{length}, $cd->{pattern}, $cd->{thin}, $pix_per_bp);
	$x += $x_step;
    }

    # Draw tracks
    foreach my $track (@$tracks) {
	my $side = $track->{side} || 0;
	my $feat_spacer = $track->{feat_spacer} || $self->feat_spacer;
	my $track_spacer = $track->{track_spacer} || $self->track_spacer;
	my $feat_width = $track->{feat_width} || $self->feat_width;
	my $feat_min_span = $track->{feat_min_span} || $self->feat_min_span;
	my $feat_bump_threshold = $track->{feat_bump_threshold} ||
	    $self->feat_bump_threshold;
	my $feat_step = $feat_width + $feat_spacer;
	my $track_type = $track->{type} || '';

	if($track_type eq 'histogram') {
	    my @h_bars;
	    my $min_score = $track->{min_score} || 0;
	    my $max_score = $track->{max_score} || 100;
	    my $h_scale = ($max_score > $min_score) ?
			($feat_width / ($max_score - $min_score)) : 1;
	    my $prev_chr = "";
	    my $max_height = 0;
	    my $x1;

	    # determine color
	    my $color = $self->{colors}->{$track->{color} || 'black'};

	    # Go through features sorted by 1. chromosome, 2. start position
	    foreach my $f (sort {$a->{chr} cmp $b->{chr} or
			         $a->{start} <=> $b->{start} }
			    @{$track->{features}})
	    {
		my $chr = $f->{chr};
		if($prev_chr ne $chr) {
		    croak("Chromosome $chr not in chromosome list") unless($chr_data{$chr});
 		    $chr_data{$prev_chr}->{fcol_bases}->[$side] += ($max_height+$feat_spacer);  
		    $self->_draw_filled_bars(\@h_bars, $color) if (@h_bars);
   
		    # reset stuff
		    @h_bars = ();
		    $prev_chr = $chr;
		    $max_height = 0;

		    # determine which column we go in
		    my $fcol_maxends = $chr_data{$chr}->{fcol_maxends};
		    #my $i = 0;
		    #$i = $chr_data{$chr}->{fcol_bases}->[$side];
		    my $fcol_offset = $chr_data{$chr}->{fcol_bases}->[$side];
		    $x1 = $x_offset +($chr_data{$chr}->{column}+1)*$x_step
		      + ($side ?
			(+$self->chr_width-1 +$self->chr_feat_spacer +$fcol_offset) :
			(-$self->chr_feat_spacer -$fcol_offset));
		}

		# set coordinates
		my $y1 = $y_offset + $f->{start} * $pix_per_bp;
		my $y2 = $y_offset + $f->{end} * $pix_per_bp;
		my $height = $f->{score} * $h_scale;
		my $x2 = $side ?
		    ($x1 + $height) :
		    ($x1 - $height);
		push @h_bars, [$x1, $y1, $x2, $y2];
		$max_height = $height if($max_height < $height);
	    }
	    $chr_data{$prev_chr}->{fcol_bases}->[$side] += ($max_height+$feat_spacer);  
	    $self->_draw_filled_bars(\@h_bars, $color) if (@h_bars);

	}
	else {
	    # Go through features sorted by 1. chromosome, 2. start position
	    foreach my $f (sort { $a->{chr} cmp $b->{chr} or
				$a->{start} <=> $b->{start} }
			@{$track->{features}})
	    {
		# determine color
		my $color = $self->{colors}->{$f->{color} || 'black'};
		unless(defined($color)) {
		    warn "Unrecognized color ".$f->{color}."; skipping feature";
		    next;
		}

		# determine start and end pos
		my $start = $y_offset + $f->{start} * $pix_per_bp;
		my $end = $y_offset + $f->{end} * $pix_per_bp;
		$end = $start + $feat_min_span-1 if
		    ($end < $start + $feat_min_span-1);

		# determine which column we go in
		my $chr = $f->{chr};
		croak("Chromosome $chr not in chromosome list") unless($chr_data{$chr});
		my $fcol_maxends = $chr_data{$chr}->{fcol_maxends};
		my $i = 0;
		if($feat_bump_threshold > 0) {
		    while(defined($fcol_maxends->[$i]) and
		        $fcol_maxends->[$i] >= $start-$feat_bump_threshold)
		    {
		        $i++;
		    }
		}
		$fcol_maxends->[$i] = $end;
		#$i += $chr_data{$chr}->{fcol_bases}->[$side];
		my $fcol_offset = $chr_data{$chr}->{fcol_bases}->[$side];
		my $x = $x_offset +($chr_data{$chr}->{column}+1)*$x_step
		+ ($side ?
		    (+$self->chr_width-1 +$self->chr_feat_spacer +$i*$feat_step + $fcol_offset) :
		    (-$feat_width -$self->chr_feat_spacer -$i*$feat_step - $fcol_offset));

		# draw
		$img->filledRectangle($x, $start, $x+$feat_width-1, $end, $color);
	    }
	    foreach my $chr_data (values %chr_data) {
		$chr_data->{fcol_bases}->[$side] += @{$chr_data->{fcol_maxends}}*$feat_step + $track_spacer;
		$chr_data->{fcol_maxends} = [];
	    }
	}
    }

}


=head2 display

 Title     : display
 Usage     : See sypnosis
 Function  : Display the most recently drawn image
 Returns   : -
 Args      : -

 This should probably be rewritten so that draw() returns an image object
 that has the capability to display itself.

=cut

sub display
{
    my ($self) = @_;
    my $pid;
    if($pid = open (CHILD, "|-")) {
	CHILD->autoflush(1);
	#CHILD->binmode();
	if ($self->{svg}){
	  print CHILD $self->{img}->svg;   
	}
	else{
	    print CHILD $self->{img}->png;
	    
	}
    }
    else {
	die "cannot fork: $!" unless defined($pid);
	system("display -");
	exit;
    }
}


=head2 display

 Title     : display
 Usage     : See sypnosis
 Function  : Opens the most recently drawn image as png or svg on screen
 Returns   : -
 Args      : -

 This should probably be rewritten so that draw() returns an image object
 that has the capability to return a png.

=cut

=head2 as_png

 Title     : as_png
 Usage     : $image->as_png
 Function  : Draws a pixel version of the image and returns a binary png string
 Returns   : -
 Args      : -
 Note that the argument svg in the constructor must not be true for this to work.
 This should probably be rewritten so that draw() returns an image object
 that has the capability to return a png.

=cut

sub as_png
{
    my $self = shift;
    croak 'Cannot create PNG image in SVG mode' if($self->svg);
    return $self->{img}->png;
}


=head2 as_svg

 Title     : as_svg
 Usage     : $image->as_svg
 Function  : Returns a SVG string: a vector version of the drawn image
 Returns   : -
 Args      : -
 Note that the argument svg in the constructor must be true for this to work.
 This should probably be rewritten so that draw() returns an image object
 that has the capability to return a svg string.

=cut


sub as_svg
{
    my $self = shift;
    croak 'Cannot create SVG string in bitmap mode' unless($self->svg);
    return $self->{img}->svg;
}


sub add_legend
{
    my ($self, $x, $y, $text_list, $color_list) = @_;
    
    my $img = $self->img;
    my $i = 0;
    my $n = @$text_list;
    die "sizes of text_list and color_list do not match in _draw_legend" unless($n == @$color_list);
    for my $i (0..@$text_list-1) {
	my $text = $text_list->[$i];
	my $col = $self->colors->{$color_list->[$i]};
	$img->filledRectangle($x, $y, $x+10, $y+10, $col);
	$img->string($self->{_font_pkg}->MediumBold, $x+15, $y, $text, $self->colors->{black});
	$y += 15;
    }
}


sub _draw_chromosome {
    my ($self, $x1, $y_offset, $length, $features, $thin, $pix_per_bp) = @_;
    my $img = $self->{img};
    my $chrlen_pix = $length * $pix_per_bp;
    my $x2 = $x1 + $self->{chr_width} -1;

    $img->filledRectangle($x1, $y_offset-1,
   		         $x2, $y_offset + $chrlen_pix+1,
		         $self->{colors}->{$self->{default_chr_color}});

    if($features) {
	print STDERR "drawing ", scalar(@$features), " features on chr\n";
	my $fx1 = $x1;
	my $fx2 = $x2;
	if($fx2 - $fx1 > 1) {
	    $fx1++;
	    $fx2--;
	}
	foreach my $f (sort { $a->{start} <=> $b->{start} } @$features) {
	    my $y1 = $f->{start} * $pix_per_bp + $y_offset; 
	    my $y2 = $f->{end} * $pix_per_bp + $y_offset;
	    my $col = $f->{color};
	    $img->filledRectangle($fx1, $y1, $fx2, $y2,
				  $self->{colors}->{$col});
	    
	}
    }

    # hack: redraw at thin regions
    my $tx = int(($x1+$x2)/2);
    foreach my $t (@$thin) {
	my $y1 = $t->{start} * $pix_per_bp + $y_offset; 
	my $y2 = $t->{end} * $pix_per_bp + $y_offset;
	$img->filledRectangle($x1, $y1, $x2, $y2,
			      $self->{colors}->{'white'});
	$img->filledRectangle($tx, $y1, $tx, $y2,
		   $self->{colors}->{$self->{default_chr_color}});

    }

}


sub _draw_filled_bars {
  my ($self, $parts, $fgcolor) = @_;
  my $img = $self->{img};

  print STDERR "drawing ",scalar(@$parts)," bars\n";

  # draw each of the component lines of the histogram surface
  for (my $i = 0; $i < @$parts; $i++) {
    my $part = $parts->[$i];
    my ($x1,$y1,$x2,$y2) = @$part;
    for my $y ($y1..$y2) {
        $img->line($x1,$y,$x2,$y,$fgcolor);
    }
  }
}


# based on Bio::Graphics::Glyph::xyplot::_draw_histogram
sub _draw_open_bars {
  my ($self, $parts, $fgcolor) = @_;
  my $img = $self->{img};

  print STDERR "drawing ",scalar(@$parts)," bars\n";

  # draw each of the component lines of the histogram surface
  for (my $i = 0; $i < @$parts; $i++) {
    my $part = $parts->[$i];
    my $next = $parts->[$i+1];
    my ($x1,$y1,$x2,$y2) = @$part;
    $img->line($x2,$y1,$x2,$y2,$fgcolor);
    next unless $next;
    my ($x3,$y3,$x4,$y4) = @$next;
    if ($y2 == $y3) {# connect vertically to next level
      $img->line($x2,$y2,$x4,$y2,$fgcolor); 
    } else {
      $img->line($x2,$y2,$x1,$y2,$fgcolor); # to bottom
      $img->line($x1,$y2,$x1,$y3,$fgcolor);                    # to right
      $img->line($x3,$y3,$x4,$y3,$fgcolor);   # up
    }
  }

  # end points: from bottom to first
  my ($x1,$y1,$x2,$y2) = @{$parts->[0]};
  $img->line($x1,$y1,$x2,$y1,$fgcolor);
  # from last to bottom
  my ($x3,$y3,$x4,$y4) = @{$parts->[-1]};
  $img->line($x4,$y4,$x3,$y4,$fgcolor);
}



1;
