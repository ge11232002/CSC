package Bio::Graphics::Glyph::fast_xyplot;

use strict;
use Benchmark;
#use GD 'gdTinyFont';

use base qw(Bio::Graphics::Glyph::xyplot);

sub minmax {
  my ($self, $scores_list) = @_;

  return $self->SUPER::minmax(@_) if (@$scores_list == 0 or ref($scores_list->[0]) != 'ARRAY');

  my $max_score = $self->option('max_score');
  my $min_score = $self->option('min_score');

  my $do_min = !defined $min_score;
  my $do_max = !defined $max_score;

  if ($do_min or $do_max) {
      for my $scores (@$scores_list) {
	  for my $s (@$scores) {
	      next unless defined $s;
	      $max_score = $s if $do_max && (!defined $max_score or $s > $max_score);
	      $min_score = $s if $do_min && (!defined $min_score or $s < $min_score);
	  }
      }
  }

  ($min_score,$max_score);
}


sub draw {
  my $self = shift;

  my ($gd,$dx,$dy) = @_;
  my ($left,$top,$right,$bottom) = $self->calculate_boundaries($dx,$dy);

  my $feature = $self->feature;
  
  my $score_start = $feature->attributes('score_start');
  my $score_span = $feature->attributes('score_span');
  my @scores_list = $feature->attributes('score_list');
  my @color_list = $feature->attributes('color_list');

  my ($min_score,$max_score) = $self->minmax(\@scores_list);

  my $side = $self->_determine_side();

  my $height = $self->height;
  my $scale  = $max_score > $min_score ? $height/($max_score-$min_score)
                                       : 1;
  my $x = $left;
  my $y = $top + $self->pad_top;

  # position of "0" on the scale
  my $y_origin = $min_score <= 0 ? $bottom - (0 - $min_score) * $scale : $bottom;
  $y_origin    = $top if $max_score < 0;

  my $clip_ok = $self->option('clip');
  $self->{_clip_ok}   = $clip_ok;
  $self->{_scale}     = $scale;
  $self->{_min_score} = $min_score;
  $self->{_max_score} = $max_score;
  $self->{_top}       = $top;
  $self->{_bottom}    = $bottom;

  # Draw each set of scores
  if(@scores_list) {
      foreach my $scores (@scores_list) {
	  # compute the y position corresponding to each score
	  my @y_pos;
	  foreach my $s (@$scores) {
	      next unless defined $s;
	      push @y_pos, $self->score2position($s);
	  }
	  
	  # call a polygon drawing method
	  my $color = shift @color_list;
	  die "colors must be references to three-elements arrays" unless($color and @$color == 3);
	  my $color_idx = $gd->colorResolve(@$color);
	  $color_idx = $gd->colorResolve(0,0,0) unless(defined $color_idx);
	  $self->_draw_filled($gd, $x, $y, $right, $y_origin, $color_idx,
			      $score_start, $score_span, \@y_pos);
      }
  }

  # draw scale, label and description
  $self->_draw_scale($gd,$scale,$min_score,$max_score,$dx,$dy,$y_origin);
  $self->draw_label(@_)       if $self->option('label');
  $self->draw_description(@_) if $self->option('description');
}


sub score2position {
  my $self  = shift;
  my $score = shift;

  return unless defined $score;

  if ($self->{_clip_ok} && $score < $self->{_min_score}) {
    return $self->{_bottom};
  }

  elsif ($self->{_clip_ok} && $score > $self->{_max_score}) {
    return $self->{_top};
  }

  else {
    my $position      = int( ($score-$self->{_min_score}) * $self->{_scale} +.5);
    return $self->{_bottom} - $position;
  }
}


sub log10 { log(shift)/log(10) }

sub max10 {
  my $a = shift;
  return 0 if $a==0;
  return -min10(-$a) if $a<0;
  return max10($a*10)/10 if $a < 1;
  
  my $l=int(log10($a));
  $l = 10**$l; 
  my $r = $a/$l;
  return $r*$l if int($r) == $r;
  return $l*int(($a+$l)/$l);
}

sub min10 {
  my $a = shift;
  return 0 if $a==0;
  return -max10(-$a) if $a<0;
  return min10($a*10)/10 if $a < 1;
  
  my $l=int(log10($a));
  $l = 10**$l; 
  my $r = $a/$l; 
  return $r*$l if int($r) == $r;
  return $l*int($a/$l);
}


sub timediff_str
{
    my ($b1, $b2) = @_;
    my $t = timestr(timediff($b2, $b1));
    return "[$t sec]";
}


sub _draw_lines {
  my $self = shift;
  my ($gd, $left, $top, $right, $y_origin, $first_start, $span, $y_values) = @_;

  my $fgcolor = $self->fgcolor;
  my $bgcolor = $self->bgcolor;
  my $scale = $self->scale;
  my $step = $span * $scale;

  # we should check at some point that first_start is larger or equal to the segment start
  # we may also want to deal with the flip case

  my $x1 = ($first_start - $self->start) * $scale + $left;
  my $x2 = $x1 + $step;
  my $current_x = ($x1+$x2)/2;
  my $current_y = $y_values->[0];
  my $n = @$y_values;

  for my $next_y ((@$y_values)[1..$n-1]) {
    my $next_x = $current_x + $step;
    $gd->line($current_x, $current_y, $next_x, $next_y, $fgcolor);
    ($current_x,$current_y) = ($next_x,$next_y);
  }

}


sub _draw_filled {
  my $self = shift;
  my ($gd, $left, $top, $right, $y_origin, $color, $first_start, $span, $y_values) = @_;

  my $scale = $self->scale;
  my $step = $span * $scale;
  my $extension = $step > 1 ? $step - 1 : 0;
  my $size = $right - $left + 1;

  #my $t0 = Benchmark->new;
  my @y_at_x = ($y_origin) x $size;
  my $x1 = ($first_start - $self->feature->start) * $scale;
  foreach my $y (@$y_values) {
      my $x2 = $x1 + $extension;
      for my $x (int($x1+.5)..int($x2+.5)) {
	  next if($x < 0 or $x >= $size); # out of bounds check
	  $y_at_x[$x] = $y if($y < $y_at_x[$x]);
      }
      $x1 += $step;
  }

#  print STDERR "Computed pixel array ".timediff_str($t0,Benchmark->new)."\n";

  #$t0 = Benchmark->new;
  if($self->{flip}) {
      my $x = $right;
      foreach my $y (@y_at_x) {
	  $gd->line($x, $y_origin, $x, $y, $color);
	  $x--;
      }  
  }
  else {
      my $x = $left;
      foreach my $y (@y_at_x) {
	  $gd->line($x, $y_origin, $x, $y, $color);
	  $x++;
      }  
  }
#  print STDERR "Drew lines ".timediff_str($t0,Benchmark->new)."\n";
}


# Main drawing code from jk:

#x1 = ((wi->start-seqStart) + (dataOffset * usingDataSpan)) * pixelsPerBase;
#x2 = x1 + (usingDataSpan * pixelsPerBase);
#for (i = x1; i <= x2; ++i)
#{
#    int xCoord = preDrawZero + i;
#    if ((xCoord >= 0) && (xCoord < preDrawSize))
#    {
#	double dataValue =
#	    BIN_TO_VALUE(datum,wi->lowerLimit,wi->dataRange);
#	
#	++preDraw[xCoord].count;
#	if (dataValue > preDraw[xCoord].max)
#	    preDraw[xCoord].max = dataValue;
#	if (dataValue < preDraw[xCoord].min)
#	    preDraw[xCoord].min = dataValue;
#	preDraw[xCoord].sumData += dataValue;
#	preDraw[xCoord].sumSquares += dataValue * dataValue;
#    }
#}
#}



1;

__END__



