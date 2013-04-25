package Bio::Graphics::Glyph::my_text;

use strict;
use base qw(Bio::Graphics::Glyph::generic);

sub default_text
{
  return "Advertising space";  
}

sub default_text_pad
{
  return 3;  
}

sub title { "none" }


sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);

  $self->configure(-link => "hello",
		   -title => "hellooo");
  
  my $fg = $self->fgcolor;
  
  my $font = $self->option('labelfont') || $self->font;
  
  my $text = defined $self->option('text') ? $self->option('text') : $self->default_text();
  my $text_pad = defined $self->option('text_pad') ? $self->option('text_pad') : $self->default_text_pad();
  
  my $height = $font->height;

  my $midY = ($y2+$y1) / 2;

  $gd->string($font, $x1+$text_pad, $midY-$height/2, $text, $self->fontcolor);
}


1;

