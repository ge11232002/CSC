package CONSITE::Alignment;


use vars qw(@ISA $AUTOLOAD);
use CGI qw(:standard);
use TFBS::Matrix::PWM;
use TFBS::Matrix::PFM;
use TFBS::Matrix;
use Bio::Seq;
use Bio::SimpleAlign;
use DBI;
use GD;
use PDL;
use IO::String;
use Class::MethodMaker
    get_set => [qw(cutoff alignseq1 alignseq2 fsiteset1 fsiteset2 
		   seq1name seq2name conservation1 start_at end_at window)],
    list    => [qw(alignseq1array alignseq2array)],
    hash    => [qw(seq1siteset seq2siteset)],
    new_with_init => ['new'];
#use TFBS::Matrix::ICM;

use strict;

use constant WEBSERVER =>'http://forkhead.cgr.ki.se';
use constant TMPDIR =>'/home/httpd/TEMP/';
use constant IMAGE_WIDTH =>600;
use constant IMAGE_MARGIN =>50;

sub init  {

    # this is ugly; it simply cries for refactoring 
    my ($self, %args) = @_;
    # print @_;
    $self->{job} = $args{-job};
    $self->{seqobj1} = $args{-seqobj1};
    $self->window($args{-window} or 50);

    $self->_parse_alignment($args{-alignstring} );
    $self->_calculate_conservation($self->window());
    $self->_do_sitesearch(($args{-MatrixSet}), "80%", 
			  $self->cutoff($args{-cutoff} or 60));
    $self->start_at($args{-start_at} or 1);
    $self->end_at($args{-end_at} or length(_strip_gaps($self->alignseq1())));
    return $self;
    
}

sub DESTROY {
    my $self = shift;
    delete $self->{job};
}

sub HTML_table {
    my ($self) = @_;
    my $table = Tr (th({-rowspan=>2},"Transcription factor").
		    th({-colspan=>5}, [$self->seq1name, $self->seq2name]));
    $table .= Tr(td( [qw(Sequence From To Score Strand 
			    Sequence From To Score Strand)]))."\n";
    $self->fsiteset1->reset;
    $self->fsiteset2->reset;
    my $site1 = $self->fsiteset1->next_site();
    my $site2 = $self->fsiteset2->next_site();
    while (defined $site1 and defined $site2)  {
	my $color = "white";
	
	if ($self->pdlindex($site1->start(), 1=>3)) { #exon
	    $color = "#ffff7f";
	}
	my $Name = $site1->Matrix->{name};
	my $Url = WEBSERVER."/cgi-bin/jaspartf?ID=".$site1->Matrix->{ID}.
	    "&score1=".$site1->score().
		"&score2=".$site2->score().
		    "&pos1=".$site1->start().
			"&pos2=".$site2->start().
			    "&seq1=".$self->seq1name().
				"&seq2=".$self->seq2name();
	$table .= Tr(td({-bgcolor=>$color},[
			 a({ -onClick=> 
"window.open('$Url', 'TFprofile', 'height=400,width=500,toolbar=no,menubar=no,scrollbars=yes,resizable=yes'); return false;", 
			     -href=>$Url, -target=>'_blank'
			   },$site1->Matrix->{name}),
			 $site1->siteseq(),
			 $site1->start(),
			 $site1->end(),
			 $site1->score(),
			 $site1->strand(),
			 $site2->siteseq(),
			 $site2->start(),
			 $site2->end(),
			 $site2->score(),
			 $site2->strand()
			]
			))
	    ."\n";
	$site1 = $self->fsiteset1->next_site();
	$site2 = $self->fsiteset2->next_site();
    }
    return table({-border=>1}, $table);
    
}


sub HTML_alignment {
    my ($self) = @_;
    
    my $TITLE_AREA  = (length($self->seq1name) >length($self->seq2name)) ?
	length($self->seq1name) : length($self->seq2name);
    my $TICKS_EVERY = 10;
    my $BLOCK_WIDTH = 60;
    my $SITE_SHIFT  = 3;
    my ($lastindex1, $lastindex2, );
    my $blockstart=1;
    my $outstring ="";
    $self->fsiteset1->reset();
    $self->fsiteset2->reset();
    my $site1 = $self->fsiteset1->next_site();
    my $site2 = $self->fsiteset2->next_site();
    
    do {
	my ($seq1slice, $seq2slice, $identity_bars, $ruler_bars, $ruler);
	my $seq1_ruler = " "x($TITLE_AREA +4);
	my $seq2_ruler = " "x($TITLE_AREA +4);
	my $seq1_ruler_bars;
	my $seq2_ruler_bars;
	my $span_open=0;
	for my $i ($blockstart..($blockstart+$BLOCK_WIDTH-1)) {
	    last unless defined (($self->alignseq1array)[$i]);
	    print STDERR "POINT 1 I=$i\n";
	    # record the actual sequences and conservation bar
	    
	    $seq1slice .= ($self->alignseq1array)[$i]; 
	    $seq2slice .= ($self->alignseq2array)[$i]; 

	    # ruler and ticks
	    
            # alignment ruler - currently not used
	    #if (($i % $TICKS_EVERY) == 0)  {
	    #	$ruler_bars .= sprintf("%".$TICKS_EVERY."s", '|');
	    #	$ruler      .= sprintf("%".$TICKS_EVERY."d", $i);
	    #}
	    
	    # seq1 ruler
	    my $pos1 = $self->pdlindex($i, 0=>1);
	    if (($pos1  % $TICKS_EVERY) == 0
		and $pos1 !=0)  {
		$seq1_ruler_bars .= '|';
		substr ($seq1_ruler, 1-length($pos1), length($pos1)) = $pos1;
		
	    } else  {
		$seq1_ruler_bars .= " ";
		$seq1_ruler .= " ";
	    }


	    # identity bars with cutoff shading

	    if ($self->conservation1->[$pos1] >= $self->cutoff 
		and ($i == $blockstart 
		     or $self->conservation1->[$pos1-1] < $self->cutoff)) 
	    {
		$identity_bars .= '<SPAN style="background-color:#dfdfdf">';
		$span_open= 1;
	    }

	    $identity_bars .= (
		(($self->alignseq1array)[$i] eq ($self->alignseq2array)[$i])
		    ? '|' : ' ');

	    if ($self->conservation1->[$pos1] >= $self->cutoff 
		and ($i == $blockstart+$BLOCK_WIDTH-1
		     or $self->conservation1->[$pos1+1] < $self->cutoff)) 
	    {
		$identity_bars .= '</SPAN>'; $span_open = 0;
	    }
	
	    # note the sequence number of the most recent nucleotide

	    #if (my $x = $self->pdlindex($i, 0 => 1)) { $lastindex1 = $x ;}
	    #if (my $x = $self->pdlindex($i, 0 => 2)) { $lastindex2 = $x ;}


	    # seq2 ruler
	    my $pos2 = $self->pdlindex($i, 0=>2);
	    if (($pos2  % $TICKS_EVERY) == 0
		and $pos2 !=0)  {
		$seq2_ruler_bars .= '|';
		substr ($seq2_ruler, 1-length($pos2), length($pos2)) = $pos2;
		
	    } else  {
		$seq2_ruler_bars .= " ";
		$seq2_ruler .= " ";
	    }
	    
	}
	if ($span_open) {$identity_bars .= '</SPAN>';} #dirty fix
	# the sites
	
	my (@site1strings, @site2strings, %url1, %url2);
	while ($site1 
	       and 
	       ((my $sitepos = $self->pdlindex($site1->start(), 1=>0)) 
		< 
	       ( $blockstart + $BLOCK_WIDTH - $SITE_SHIFT)))  {
	    print STDERR "POINT 2 site=".$site1->Matrix->{name}."\n";

	    my $strlen; my $is_free = 0; my $offset=0;
	    while (($strlen =(length($site1strings[$is_free]) or 0))
		    > 
		    ( $offset = ($sitepos+$SITE_SHIFT) % $BLOCK_WIDTH ) 
		    )
	    {  	    print STDERR "POINT 2.5 isfree= $is_free strlen=$strlen offset=$offset $site1strings[$is_free]\n";
		    $is_free++; last if $offset <0;}
	    my $site_id; 

	    do 
	    {$site_id = sprintf("%03d", rand(1000))} 
	    while defined $url1{$site_id};

	    print STDERR "OFFSET: $offset   STRLEN: $strlen \n" ;
	    $site1strings[$is_free] .= (" "x($offset-$strlen-1)).
		$site_id.$site1->strand().$site1->siteseq().
		    ":".$site1->Matrix->{name};
	    $site2strings[$is_free] .= (" "x($offset-$strlen-1)).
		$site_id.$site2->strand().$site2->siteseq().
		    ":".$site2->Matrix->{name};
	    $url1{$site_id} = $self->site_url($site1, $site2, $site1);
	    $url2{$site_id} = $self->site_url($site1, $site2, $site2);
	    

	    ($site1 = $self->fsiteset1->next_site()) || last;
	    ($site2 = $self->fsiteset2->next_site()) || last;
	    print STDERR "POINT 3 site2=".$site2->Matrix->{name}."\n";
		
	}
	foreach (0..$#site1strings)  {
	    $site1strings[$_] =~ s/(\d{3})(\+|\-)\S+/$url1{$1} /g;
	    $site2strings[$_] =~ s/(\d{3})(\+|\-)\S+/$url2{$1} /g;
	    print STDERR "POINT 4 site2strings=".$site2->Matrix->{name}."\n";
	}
	

	$outstring .= join 
	    ("", hr,br,"&nbsp;\n",
	     (map {(" "x($TITLE_AREA + 4-$SITE_SHIFT). $_."\n")} 
	          reverse(@site1strings)),
	     "$seq1_ruler\n",
	     sprintf("%".$TITLE_AREA."s    %-".$BLOCK_WIDTH."s  %s\n",
		     " ", $seq1_ruler_bars, " "),
	     sprintf("<I>%".$TITLE_AREA."s</I>    %-".$BLOCK_WIDTH."s  %s\n",
		     $self->seq1name(), $seq1slice, " "),
	     sprintf("%".$TITLE_AREA."s    %-".$BLOCK_WIDTH."s  %s\n",
		     "", $identity_bars, ""),
	     sprintf("<I>%".$TITLE_AREA."s</I>    %-".$BLOCK_WIDTH."s  %s\n",
		     $self->seq2name(), $seq2slice, " "),
	     sprintf("%".$TITLE_AREA."s    %-".$BLOCK_WIDTH."s  %s\n",
		     " ", $seq2_ruler_bars, " "),
	     "$seq2_ruler\n", 
	     (map {(" "x($TITLE_AREA + 4-$SITE_SHIFT). $_."\n")} @site2strings),
	     "\n\n");
	$blockstart += $BLOCK_WIDTH;
	
    } while ($blockstart <= length $self->alignseq1());
     print STDERR "POINT 5 : EXIT\n";

    return $outstring;
}

sub site_url {

    my ($self, $site1, $site2, $site) = @_;
    	my $URL = WEBSERVER."/cgi-bin/jaspartf?ID=".$site1->Matrix->{ID}.
	    "&score1=".$site1->score().
		"&score2=".$site2->score().
		    "&pos1=".$site1->start().
			"&pos2=".$site2->start().
			    "&seq1=".$self->seq1name().
				"&seq2=".$self->seq2name();

    my $url = b({-style=>"background-color:#dfdfdf"},$site->siteseq()).":".
	a({-href=>WEBSERVER."/cgi-bin/jaspartf?ID=".$site->Matrix->{ID},
	   -onClick  =>"window.open('$URL', 'TFprofile', 'height=400,width=500,toolbar=no,menubar=no,scrollbars=yes,resizable=yes'); return false;", 

	   -onMouseOver=>
	       "status='TF: ".$site->Matrix->{name}. "; ".
	       "Structural class: ".$site->Matrix->{class}. "; ".
	       "Score ".$site->score()." (".
		   int(($site->score()-$site->Matrix->{min_score})/
		       ($site->Matrix->{max_score}-$site->Matrix->{min_score})
		       *100
		       +0.5).
	       "%)'; return true;",
		   -onMouseOut=>"status=''"},
	  $site->Matrix->{name}).
	"[".$site->strand()."]";
    return $url;

}	    


sub maparea_factor  {
    my ($self, $site, $site2, $x1, $y1, $x2, $y2) = @_;
    my $URL = WEBSERVER."/cgi-bin/jaspartf?ID=".$site->Matrix->{ID}.
	"&score1=".$site->score().
	    "&score2=".$site2->score().
		"&pos1=".$site->start().
		    "&pos2=".$site2->start().
			"&seq1=".$self->seq1name().
			    "&seq2=".$self->seq2name();
    my $onClick = "window.open('$URL', 'TFprofile', 'height=400,width=500,toolbar=no,menubar=no,scrollbars=yes,resizable=yes'); return false;";
    my $onMouseOver = "status='TF: ".$site->Matrix->{name}. "; ".
	       "Structural class: ".$site->Matrix->{class}. "; ".
	       "Score ".$site->score()." (".
		   int(($site->score()-$site->Matrix->{min_score})/
		       ($site->Matrix->{max_score}-$site->Matrix->{min_score})
		       *100
		       +0.5).
			   "%)'; return true;";
    return qq!<AREA COORDS="$x1,$y1,$x2,$y2"
	       SHAPE="RECT"
	       HREF="$URL"
	       TARGET="_blank"
	       ONMOUSEOVER="$onMouseOver"
	       ONMOUSEOUT="status=''"
	       ONCLICK="$onClick">!;

}
sub HTML_image  {
    my $self = shift;
    my ($image, $imagemap) = $self->_drawgene_png();
    open (OUT,">".TMPDIR.$self->{job}->jobID().".GENE.png");
    print OUT $image->png();
    close OUT;
    $image = $self->_drawpf_png();
    open (OUT,">".TMPDIR.$self->{job}->jobID().".PF.png");
    print OUT $image->png();
    close OUT;
    print STDERR join( "::", @{$self->conservation1()});
    return table(Tr(td($imagemap.
		       img({-src=>WEBSERVER."/TEMP/".$self->{job}->jobID().".GENE.png",
			    -usemap=>"#seq1map", -border=>0}))),
		 Tr(td(img({-src=>WEBSERVER."/TEMP/".$self->{job}->jobID().".PF.png"}))));
}


sub _drawgene_png  {
    my $self = shift;
    my $image_map_areas = "";
    my $factor = (IMAGE_WIDTH - 2*IMAGE_MARGIN)
	                        /
		 ($self->end_at()-$self->start_at() +1);
    my $FONT = gdSmallFont;
    my $SPACING = 1; #vertical spacing between labels
    $self->fsiteset1->sort_reverse();
    $self->fsiteset2->sort_reverse();
    my @labelrows;
    my ($leftmostpos, $leftmostrow);
    while (my $site = $self->fsiteset1->next_site())  {
	my $xpos = $site->start() * $factor + IMAGE_MARGIN;
	my $row=0;
	my $site2=$self->fsiteset2->next_site();
	while ($labelrows[$row] 
	       and 
	       (($xpos + (length($site->Matrix->{name}) +1) * $FONT->width)
	       >
	       $leftmostpos )#$labelrows[$row][0]->{pos}))
	       and
	       ($row<=$leftmostrow))
	{
	    $row++;
	}
	($leftmostpos, $leftmostrow) = ($xpos, $row);
	unshift @{$labelrows[$row]}, {pos => $xpos, site => $site, site2 =>$site2};
    }

    # create image

    my $OFFSET = IMAGE_MARGIN +(scalar(@labelrows)+2)*($FONT->height+$SPACING);
    my $im = GD::Image->new(IMAGE_WIDTH, 
			    $OFFSET + 6*3 + IMAGE_MARGIN);
    my $white = $im->colorAllocate(255,255,255);
    my $black = $im->colorAllocate(0,0,0);
    my $gray  = $im->colorAllocate(127,127,127);
    my $red   = $im->colorAllocate(255,0,0);
    my $blue  = $im->colorAllocate(0,0,255);
    my $lightyellow = $im->colorAllocate(255,255,127);
    
    # draw exons 
    
    my $pExonPdl = 
	$self->pdlindex->slice(':,1')->where
	    ($self->pdlindex->slice(':,3') == 1);
    $pExonPdl = $pExonPdl->where($pExonPdl>0);
    print STDERR "pExonPdl $pExonPdl\n";
    my $pExonStarts = 
	$pExonPdl->slice('1:-1')->where(
	    ($pExonPdl->slice('1:-1' )- $pExonPdl->slice('0:-2')>1));
    print STDERR "pExonStarts $pExonStarts\n";
    my $pExonEnds = 
	$pExonPdl->slice('0:-2')->where(
	    ($pExonPdl->slice('1:-1' )- $pExonPdl->slice('0:-2')>1));
    print STDERR "pExonEnds $pExonEnds\n";
    my @exon_starts = ((list($factor*$pExonPdl->slice(0) +IMAGE_MARGIN)),
		       (list ($pExonStarts*$factor +IMAGE_MARGIN)));
    print STDERR "ExonStarts @exon_starts\n";
    my @exon_ends = ((list ($pExonEnds*$factor +IMAGE_MARGIN)),
		      (list($factor*$pExonPdl->slice(-1) +IMAGE_MARGIN)));
    print STDERR "ExonEnds @exon_ends\n";

    while ( my $exonstart = splice (@exon_starts, 0, 1))  {
	my $exonend = splice (@exon_ends, 0,1);
	$im->filledRectangle($exonstart, $OFFSET-1, $exonend, $OFFSET+16, 
			     $lightyellow);
    }

    # backbone
    $im->filledRectangle(IMAGE_MARGIN, $OFFSET+2, 
			 IMAGE_WIDTH - IMAGE_MARGIN, $OFFSET+5, $black);

    # put the labels
    foreach my $row(0..$#labelrows)  {
	# line corner y value:
	my $lcy = $OFFSET - ($row+1.5)*($FONT->height()+$SPACING);
	# label y coordinate:
	my $laby = $OFFSET - ($row+2)*($FONT->height()+$SPACING);
	foreach my $label (@{$labelrows[$row]})  {
	    
	    $im->line($label->{pos}, $OFFSET+1, 
		      $label->{pos}, $lcy,
		      $gray);
	    $im->line($label->{pos}, $lcy,
		      $label->{pos} + $FONT->width()/2-1, $lcy,
		      $gray);
	    my $name = $label->{site}->Matrix->{name};
	    $im->string($FONT, $label->{pos} + $FONT->width() -1, $laby, 
			$name, $blue);
	    $image_map_areas .= $self->maparea_factor( $label->{site}, $label->{site2},
				     $label->{pos}+$FONT->width() -1,
				     $laby,
				     $label->{pos}+$FONT->width()*(length($name)+1)-1,
				     $laby+$FONT->height());
	}
    }
    

    # draw seq2 
    my $pSeq2Pdl = 
	$self->pdlindex->slice(':,1')->where
	    ($self->pdlindex->slice(':,2') > 0);
    $pSeq2Pdl = $pSeq2Pdl->where($pSeq2Pdl>0);
    print STDERR "pSeq2Pdl $pSeq2Pdl\n";
    my $pSeq2Starts = 
	$pSeq2Pdl->slice('1:-1')->where(
	    ($pSeq2Pdl->slice('1:-1' )- $pSeq2Pdl->slice('0:-2')>1));
    print STDERR "pSeq2Starts $pSeq2Starts\n";
    my $pSeq2Ends = 
	$pSeq2Pdl->slice('0:-2')->where(
	    ($pSeq2Pdl->slice('1:-1' )- $pSeq2Pdl->slice('0:-2')>1));
    print STDERR "pSeq2Ends $pSeq2Ends\n";
    my @seq2_starts = (list($pSeq2Pdl->slice(0)*$factor +IMAGE_MARGIN),
		       (list ($pSeq2Starts*$factor +IMAGE_MARGIN)));
    print STDERR "Seq2starts @seq2_starts\n";
    my @seq2_ends = ((list ($pSeq2Ends*$factor +IMAGE_MARGIN),
		      list($pSeq2Pdl->slice(-1)*$factor +IMAGE_MARGIN)));
    print STDERR "Seq2starts @seq2_starts\n";
    print STDERR "Seq2ends @seq2_ends\n";

    while ( my $seq2start = splice (@seq2_starts, 0, 1))  {
	my $seq2end = splice (@seq2_ends, 0,1);
	$im->filledRectangle($seq2start, $OFFSET+10, $seq2end, $OFFSET+13, $gray);
    }
	     
    return ($im, "<MAP name='seq1map'> $image_map_areas </MAP>");

}

sub _drawpf_png  {
    my $self = shift;
    my ($XDIM, $YDIM, $XOFFSET ,$YOFFSET) = 
	(IMAGE_WIDTH, 350, IMAGE_MARGIN, 300);
    my $FONT = gdSmallFont;
    my $WINDOW = $self->window();
    my $CUTOFF = $self->cutoff();
    my $BEGIN_AT = $self->start_at();
    my $END_AT   = $self->end_at();
    my $STRAND   = 1;
    my $factor = (IMAGE_WIDTH - 2*IMAGE_MARGIN)
	                        /
		 ($self->end_at()-$self->start_at() +1);

    my $image = new GD::Image(IMAGE_WIDTH, $YDIM);

    my $white       = $image->colorAllocate(255,255,255);
    my $black       = $image->colorAllocate(0,0,0);
    my $red         = $image->colorAllocate(255,0,0);
    my $yellow      = $image->colorAllocate(200,200,0);
    my $lightyellow = $image->colorAllocate(255,255,127);
    my $green       = $image->colorAllocate(0,200,0);
    my $ltgray      = $image->colorAllocate(200,200,200);

    my $i=0;

    my @graphslice = @{$self->conservation1()}[$BEGIN_AT..$END_AT];
    my @isexon     = (list $self->pdlindex->slice(':,3')->where($self->pdlindex->slice(':,1')>0))[$BEGIN_AT..$END_AT];

    $image->string($FONT, 
		   ($XOFFSET-length($BEGIN_AT)*$FONT->width/2), 
		   $YOFFSET+$FONT->height(), 
		   $BEGIN_AT, $black);
    $image->string($FONT,
		   ($XDIM-$XOFFSET-length($END_AT)*$FONT->width/2), 
		   $YOFFSET+$FONT->height(), 
		   $END_AT, $black);
    my ($prev_x, $prev_y);
    foreach $i (0..$#graphslice){
	$image->line ($XOFFSET+$factor*$i, $YOFFSET,
		      $XOFFSET+$factor*$i, 2, $lightyellow) 
	    if $isexon[$i] and  int($XOFFSET+$factor*$i)!=int($prev_x) ;
	if ($prev_x) {
	    $image->line ($prev_x, $prev_y, ($XOFFSET+$factor*$i),
			  ($YOFFSET-(($YOFFSET-2)/100)*$graphslice[$i]),$red); 
	}
	($prev_x, $prev_y) = ($XOFFSET+$factor*$i,
			      $YOFFSET-(($YOFFSET-2)/100)*$graphslice[$i]);
    }

    ##############
    # Draw the grid
    
    my $step = (int($STRAND*($END_AT-$BEGIN_AT)/1000)*100 or 100);
    my $firstoffset=($STRAND==1 
		     ? ($step - ($BEGIN_AT % $step))
		     :($BEGIN_AT % $step));

    # print "STEP:  $step\n";
    
    for ($i=$firstoffset; $i< $STRAND*($END_AT-$BEGIN_AT); $i+=$step)  {
	my $xpos=$XOFFSET + $i*$factor;
	$image->dashedLine($xpos, $YOFFSET, $xpos, 0, $ltgray);
	$image->stringUp($FONT, $xpos, $YOFFSET-10, 
			 int($BEGIN_AT+$i*$STRAND), $black);
    # print join ("\t",($xpos, $YOFFSET, $xpos, 0, "LINE")), "\n";
    } 

#############
# Draw the 70% cutoff
    $image->dashedLine($XOFFSET, (1-$CUTOFF/100)*$YOFFSET, $XDIM-$XOFFSET, (1-$CUTOFF/100)*$YOFFSET, $ltgray);

    
###########
# Draw rectangle around graph

    $image->filledRectangle($XOFFSET-1, $YOFFSET, 
			    $XOFFSET, $YOFFSET+3, $black);
    $image->filledRectangle($XDIM-$XOFFSET, $YOFFSET, 
			    $XDIM - $XOFFSET + 1, $YOFFSET+3, $black);
    $image->rectangle($XOFFSET, $YOFFSET, $XDIM-$XOFFSET, 1, $black);
    $image->rectangle($XOFFSET-1, $YOFFSET+1, $XDIM-$XOFFSET+1, 0, $black);

    return $image;
  
}    

sub _parse_alignment {
    my ($self, $alignment) = @_;
    
    my ($seq1, $seq2, $start);
    
    my @match;
    my @alnlines = split("\n", $alignment);
    shift @alnlines;shift @alnlines;shift @alnlines; # drop header
    
    while ($_=shift @alnlines) {
	$start=1;
	# print $_; 
	my ($label1, $string1) = split; # (/\s+/($_, 21,60);
	$self->seq1name($label1) unless $self->seq1name();
  	my ($label2, $string2) = split /\s+/, shift(@alnlines); 
 	$self->seq2name($label2) unless $self->seq2name();
                     #substr($string2,21,60);
	shift @alnlines; shift @alnlines; # skip asterisk line and blank line
	$seq1 .= $string1;
	$seq2 .= $string2;
    }

    $self->alignseq1($seq1);
    $self->alignseq2($seq2);
    $self->alignseq1array( my @seq1 = ("-", split('', $seq1) ));
    $self->alignseq2array( my @seq2 = ("-", split('', $seq2) ));
    
    my (@seq1index, @seq2index);
    my ($i1, $i2) = (0, 0);
    for (0..$#seq1) {
	my ($s1, $s2) = (0, 0);
	$seq1[$_] ne "-" and  $s1 = ++$i1;
	$seq2[$_] ne "-" and  $s2 = ++$i2;
	push @seq1index, $s1;
	push @seq2index, $s2;
    }

    $self->pdlindex( pdl [ [list sequence($#seq1+1)], [@seq1index], [@seq2index], [list zeroes ($#seq1+1)] ]) ;

    # mark exons in row 3
    my $seqobj =$self->{job}->seq1();
    my @features = $seqobj->top_SeqFeatures;
    for my $exon (@features)  {
	next unless ref($exon) =~ /::Exon/;
	my $exonstart = $self->pdlindex($exon->start(), 1=>0);
	my $exonend   = $self->pdlindex($exon->end(),   1=>0);

	my $sl = $self->pdlindex->slice("$exonstart:$exonend,3");
	$sl++;
    }
    return 1;

}

sub pdlindex {
    my ($self, $input, $p1, $p2) = @_ ;
    # print ("PARAMS ", join(":", @_), "\n");
    if (ref($input) eq "PDL")  {
	$self->{pdlindex} = $input;
    }

    unless (defined $p2)  {
	return $self->{pdlindex};
    }
    else {
	my @results = list 
	    $self->{pdlindex}->xchg(0,1)->slice($p2)->where(
	                             $self->{pdlindex}->xchg(0,1)->slice($p1)==$input
						 );
	wantarray ? return @results : return $results[0];
    }
}
	
			   
sub _calculate_conservation  {
    my ($self, $WINDOW) = @_;
    my @seq1 = $self->alignseq1array();
    my @seq2 = $self->alignseq2array();

    my @match;

    while ($seq1[0] eq "-")  {
	shift @seq1;
	shift @seq2;
    }

    for my $i (0..$#seq1) {
  	push (@match,( uc($seq1[$i]) eq uc($seq2[$i]) ? 1:0)) 
  	    unless ($seq1[$i] eq "-" or $seq1[$i] eq ".");
    }
    my @graph=($match[0]);
    for my $i (1..($#match+$WINDOW/2))  {
  	$graph[$i] = $graph[$i-1] 
  	           + ($i>$#match ? 0: $match[$i]) 
  		   - ($i<$WINDOW ? 0: $match[$i-$WINDOW]);
  	print STDERR ($i, "\t$match[$i]\t$graph[$i]\n");
    }

    # at this point, the graph values are shifted $WINDOW/2 to the right
    # i.e. the score at a certain position is the score of the window 
    # UPSTREAM of it: To fix it, we shoud discard the first $WINDOW/2 scores:
    $self->conservation1 ([]);
    foreach (@graph[int($WINDOW/2)..$#graph])  {
	push @{$self->conservation1()}, 100*$_/$WINDOW;
	print STDERR ">>>". 100*$_/$WINDOW."\n";
    }
    
}


sub _strip_gaps {
    # not OO
    my $seq = shift;
    $seq =~ s/\-|\.//g;
    return $seq;
}


sub _do_sitesearch  {
     my ($self, $MATRIXSET, $THRESHOLD, $CUTOFF) = @_;

    print STDERR "SEARCHING 1...\n";
    my $seqobj1 = Bio::Seq->new(-seq=>_strip_gaps($self->alignseq1()),
				-id => "Seq1");
    my $siteset1 = 
	$MATRIXSET->csearch(-seqobj => $seqobj1,
			    -threshold => $THRESHOLD);
    $siteset1->sort();

    print STDERR "SEARCHING 2...\n";
    my $seqobj2 = Bio::Seq->new(-seq=>_strip_gaps($self->alignseq2()),
				-id => "Seq2");
    my $siteset2 =
	$MATRIXSET->csearch(-seqobj => $seqobj2,
			    -threshold => $THRESHOLD);
    $siteset2->sort();
    
    # instead of $comp_of_main($i) use $self->pdlindex($i, 1 => 2)
    # instead of @main_sites[0,1,2], use $siteset1 object
    # instead of @comp_sites, use $siteset2 object

    my $site1 = $siteset1->next_site();
    my $site2 = $siteset2->next_site();

    # initialize filtered sitesets:

    $self->fsiteset1 (TFBS::SiteSet->new());
    $self->fsiteset2 (TFBS::SiteSet->new());

    while (defined $site1 and defined $site2) {
	#my $pos1    = $site1->{start};
	#my $pos1of2 = $self->pdlindex($site2->{start}, 2=>1) or 0;
	my $pos1_in_aln = $self->pdlindex($site1->{start}, 1=>0);
	my $pos2_in_aln = $self->pdlindex($site2->{start}, 2=>0);
	# print "pos1 $pos1  pos2 $pos1of2\n";
	#next unless $pos1of2;
	my $cmp = (($pos1_in_aln <=> $pos2_in_aln) 
		   or 
		   ($site1->Matrix->{name} cmp $site2->Matrix->{name}));
	if ($cmp==0) { ### match
	    if ($self->conservation1->[$site1->{start}]
		>=
		$self->cutoff()
		)
	    {
		$self->fsiteset1->add_site($site1);
		$self->fsiteset2->add_site($site2);
	    }
	    $site1 = $siteset1->next_site();
	    $site2 = $siteset2->next_site();
	}
	elsif ($cmp<0)  { ### $siteset1 is behind
	    $site1 = $siteset1->next_site();
	}
	elsif ($cmp>0)  { ### $siteset2 is behind
	    $site2 = $siteset2->next_site();
	}	    
    }    
}

1;

