package CONSITE::SingleSeq;


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
use CONSITE::Analysis;
use Class::MethodMaker
    get_set => [qw(cutoff window
		   alignseq1
		   fsiteset1
		   seq1name
		   seq1length
		   start_at end_at
		   REL_CGI_BIN_DIR
		   ABS_TMP_DIR
		   REL_TMP_DIR
		   )],
    list    => [qw(alignseq1array)],
    hash    => [qw(seq1siteset)],
    new_with_init => ['new'];


use strict;

use constant IMAGE_WIDTH =>600;
use constant IMAGE_MARGIN =>50;

my @seq1RGB   = (48, 176, 200);
my @seq2RGB   = (31, 225,   0);
my $seq1hex   = "#30B0C8";
my $seq2hex   = "#1FD100";
my $seq1bghex = "#CCFFFF";
my $seq2bghex = "#CCFFCC";

@ISA =('CONSITE::Analysis');

sub init  {

    my ($self, %args) = @_;
    $self->{job} = $args{-job};
    $self->window($args{-window} or 50);
    $self->REL_TMP_DIR($args{-rel_tmp_dir});
    $self->ABS_TMP_DIR($args{-abs_tmp_dir});
    $self->REL_CGI_BIN_DIR($args{-rel_cgi_bin_dir});
    #$self->alignseq1( $self->{job}->seq1->seq() );
    #$self->alignseq1array_clear();
    #$self->alignseq1array_push(split("", $self->alignseq1()));
    $self->_parse_seq();
    $self->seq1length($self->{job}->seq1->length());
    $self->_do_sitesearch(($args{-MatrixSet}), $self->{job}->threshold(),
			  $self->cutoff($args{-cutoff} or 60));

    $self->_set_start_end(%args);
    return $self;

}

sub DESTROY {
    my $self = shift;
    delete $self->{job};
}

sub HTML_table {
    my ($self) = @_;
    my $table = Tr (th({-rowspan=>2},"Transcription factor").
		    th({-colspan=>6, -bgcolor=>$seq1bghex}, $self->seq1name),
		    th({-colspan=>6, -bgcolor=>$seq2bghex}, $self->seq2name));
    $table .= Tr(td({-colspan=>2,-bgcolor=>$seq1bghex}, "Sequence").
		 td({-bgcolor=>$seq1bghex},[qw(From To Score Strand)]))."\n";

    my $iterator1 = $self->fsiteset1->Iterator(-sort_by=>'start');
    # my $site1 = $self->fsiteset1->next_site();
    my ($s, $e) = ($self->start_at(),$self->end_at());
    if ($self->fsiteset1->size == 0) {
		$table = Tr(td(p({-style=>"font-size:16px", -align=>"center"}, "No binding sites detected.")));
    }

    while (my $site1 = $iterator1->next())  {
		next if $site1->start < $s;
		last if $site1->start >= $e;
		my $color = "white";

		if ($self->pdlindex($site1->start(), 1=>3) == 1) { #exon
			$color = "#ffff7f";
		}
		elsif ($self->pdlindex($site1->start(), 1=>3) == 2) { #orf
			$color = "#ff7f7f";
		}

		my $Name = $site1->Matrix->{name};
		my $Url = $self->REL_CGI_BIN_DIR."/jaspartf?ID=".$site1->Matrix->{ID}.
			"&score1=".$site1->score().
			"&pos1=".$self->_abs_pos($site1->start(),1).
				"&seq1=".$self->seq1name().
				"&name=".$site1->Matrix->{name}.
				"&jobID=".$self->{job}->jobID()."";
		$table .= Tr(td({-bgcolor=>$color},
				a({ -onClick=>
					"window.open('$Url', 'TFprofile', ".
						"'height=400,width=500,toolbar=no,menubar=no,".
					"scrollbars=yes,resizable=yes'); return false;",
					 -href=>$Url, -target=>'_blank'
					 },$site1->Matrix->{name})).
				 td({-bgcolor=>$seq1bghex}, "&nbsp;").
				 td({-bgcolor=>$color},
				[
				 _plussiteseq($site1),
				 $self->_abs_pos($site1->start(),1),
				 $self->_abs_pos($site1->end(),1),
				 $site1->score(),
				 ($site1->strand()=="-1"?"-":"+")
				])

				 )
			."\n";
    }
    return table({-border=>1, -valign=>"top",-style=>"font-size:12px"}, $table);

}


sub HTML_seqview {
    my ($self) = @_;

    my $TITLE_AREA  = (length($self->seq1name) >length($self->seq2name)) ?
	length($self->seq1name) : length($self->seq2name);
    my $TICKS_EVERY = 10;
    my $BLOCK_WIDTH = 60;
    my $SITE_SHIFT  = 3;
    my ($lastindex1, $lastindex2, );
    my $outstring ="";
    my $iterator1 = $self->fsiteset1->Iterator(-sort_by => 'start');

    my ($s, $e) = ($self->start_at(),$self->end_at());
    my $site1;
	if ($self->fsiteset1->size > 0) {
		do {$site1 = $iterator1->next()} until $site1->start >=$s;
    }
    my $blockstart= ($s or 1);
	my @alnarray1 = ($self->alignseq1array)[0..$e];


    do {
		my ($seq1slice); #, $seq2slice, $identity_bars, $ruler_bars, $ruler);
		my $seq1_ruler = " "x($TITLE_AREA +4);
		my $seq1_ruler_bars;
		my $span_open=0;
			print STDERR "ALIGNSEQ1ARRAY".join("##",($self->alignseq1array))."\n";
		for my $i ($blockstart..($blockstart+$BLOCK_WIDTH-1)) {
			last unless defined ($alnarray1[$i]);
			print STDERR "POINT 1 I=$i\n";
			# record the actual sequences and conservation bar

			$seq1slice .= $alnarray1[$i];

			# seq1 ruler
			my $pos1 = $self->pdlindex($i, 0=>1);
			if ($pos1 !=0
			 and
			(my $abs1 = $self->_abs_pos($pos1, 1))  % $TICKS_EVERY == 0)
			{
				$seq1_ruler_bars .= '|';
				substr ($seq1_ruler, 1-length($abs1), length($abs1)) = $abs1;

			} else  {
				$seq1_ruler_bars .= " ";
				$seq1_ruler .= " ";
			}

		}

		my (@site1strings, %url1);
		while ($site1
			   and
			   ((my $sitepos = $self->pdlindex($site1->start(), 1=>0))
			<
			( $blockstart + $BLOCK_WIDTH - $SITE_SHIFT)))  {

			my $strlen; my $is_free = 0; my $offset=0;
			while (($strlen =(length($site1strings[$is_free]) or 0))
			   >
			   ( $offset = ($sitepos+$SITE_SHIFT -($s%$BLOCK_WIDTH)) % $BLOCK_WIDTH +1 )
			   )
			{
				$is_free++; last if $offset <0;
			}

			my $site_id;

			# pick an unused three-digit label for the site

			do
			{$site_id = sprintf("%03d", rand(1000))}
			while defined $url1{$site_id};

			$site1strings[$is_free] .= (" "x($offset-$strlen-1)).
			$site_id.($site1->strand()==-1?"-":"+")._plussiteseq($site1).
				":".$site1->Matrix->{name};

			$url1{$site_id} = $self->site_url($site1);

			($site1 = $iterator1->next()) || last;

		}

		foreach (0..$#site1strings)  {
			$site1strings[$_] =~ s/(\d{3})(\+|\-)\S+/$url1{$1} /g;
			print STDERR "S1S $_: $site1strings[$_]\n"
		}


		$outstring .= join
			("", hr,br,"&nbsp;\n",
			 (map {(" "x($TITLE_AREA + 4-$SITE_SHIFT). $_."\n")}
			  reverse(@site1strings)),
			 sprintf("\n<I><B><FONT color=$seq1hex>%".$TITLE_AREA.
				 "s</FONT></B></I>    %-".$BLOCK_WIDTH."s  %s\n",
				 $self->seq1name(), $seq1slice, " "),
			 sprintf("%".$TITLE_AREA."s    %-".$BLOCK_WIDTH."s  %s\n",
				 " ", $seq1_ruler_bars, " "),
			 "$seq1_ruler\n",

			 "\n\n");
		$blockstart += $BLOCK_WIDTH;

    } while ($blockstart <= $e);


    return $outstring;
}

sub _plussiteseq  {
    print STDERR "In _plussitestring \n";
    # a utility function
    my $site = shift;
    if ($site->strand == -1)  {
	return Bio::Seq->new(-seq=>$site->siteseq,
			     -moltype=>'dna')->revcom->seq;
    }
    else {
	return $site->siteseq;
    }
}


sub site_url {

    my ($self, $site) = @_;
    my $URL = $self->REL_CGI_BIN_DIR."/jaspartf?ID=".$site->Matrix->{ID}.
	"&score1=".$site->score().
	    "&pos1=".$self->_abs_pos($site->start(), 1).
		"&seq1=".$self->seq1name().
		    "&name=".$site->Matrix->{name}.
		    "&jobID=".$self->{job}->jobID()."";

    my $url = b({-style=>"background-color:#dfdfdf"},_plussiteseq($site)).":".
	a({-href=>$self->REL_CGI_BIN_DIR."/jaspartf?ID=".$site->Matrix->{ID}.
	       "&jobID=".$self->{job}->jobID()."",
	   -onClick  =>
	       "window.open('$URL', 'TFprofile', ".
	       "'height=400,width=500,toolbar=no,menubar=no,".
		   "scrollbars=yes,resizable=yes'); return false;",

	   -onMouseOver=>
	       "status='TF: ".$site->Matrix->{name}. "; ".
	       "Structural class: ".$site->Matrix->{class}. "; ".
	       "Score ".$site->score()." (".
		   int(($site->score()-$site->Matrix->{min_score})/
		       ($site->Matrix->{max_score}-$site->Matrix->{min_score})
		       *100
		       +0.5).
	       "%)';",
		   -onMouseOut=>"status=''"},
	  $site->Matrix->{name}).
	"[".($site->strand()==-1 ?"-":"+")."]";
    return $url;

}


sub maparea_factor  {
    my ($self, $which, $site, $site2, $x1, $y1, $x2, $y2) = @_;
    my $URL = $self->REL_CGI_BIN_DIR."/jaspartf?ID=".$site->Matrix->{ID}.
	"&score1=".$site->score().
		"&pos1=".$self->_abs_pos($site->start(),$which).
			"&seq1=".$self->seq1name().
			    "&name=".$site->Matrix->{name}.
			    "&jobID=".$self->{job}->jobID()."";
    my $onClick = "window.open('$URL', 'TFprofile',
    'height=400,width=500,toolbar=no,menubar=no,scrollbars=yes,resizable=yes');
     return false;";
    my $onMouseOver = "TF: ".$site->Matrix->{name}. "; ".
	       "Structural class: ".$site->Matrix->{class}. "; ".
	       "Score ".$site->score()." (".
		   int(($site->score()-$site->Matrix->{min_score})/
		       ($site->Matrix->{max_score}-$site->Matrix->{min_score})
		       *100
		       +0.5).
			   "%)";
    return qq!<AREA COORDS="$x1,$y1,$x2,$y2"
	       SHAPE="RECT"
	       HREF="$URL"
	       TARGET="_blank"
	       ALT="$onMouseOver"
	       onMouseOver="status='$onMouseOver';"
	       onMouseOut="status=''"
	       ONCLICK="$onClick">!;

}

sub HTML_image  {
    my $self = shift;
    my $page_id = sprintf("%04d", int rand(10000));
    my ($image, $image1map) = $self->_drawgene_png();
    print STDERR "IMAGE PATH: ".$self->ABS_TMP_DIR."/".$self->{job}->jobID().$page_id.".GENE1.png";
    open (OUT,">".$self->ABS_TMP_DIR."/".$self->{job}->jobID().$page_id.".GENE1.png");
    print OUT $image->png();
    close OUT;
    return table
	(Tr(td({-align=>"center"},
	       hr,font({-size=>"+1"},
		       b("Putative transcription factor binding sites ".
			 "found along ".
			 font({-color=>$seq1hex},i($self->seq1name)))))),
	 Tr(td($image1map.
	       img({-src=>$self->REL_TMP_DIR."/".$self->{job}->jobID().
			$page_id.".GENE1.png",
			-usemap=>"#seq1map", -border=>0}))));
}

sub _drawaln_png {
    my $self = shift;
    my $image_map_areas = "";

    # this is factor for the whole alignment: the denominator is the
    # entire length of the alignment

    my $factor = (IMAGE_WIDTH - 2*IMAGE_MARGIN)
	                        /
		 ($self->pdlindex->getdim(0));
    my $FONT = gdSmallFont;

    # create image

    my $OFFSET = IMAGE_MARGIN;# +(scalar(@labelrows)+2)*($FONT->height+$SPACING);
    my $im = GD::Image->new(IMAGE_WIDTH,
			    $OFFSET + 6*3 + IMAGE_MARGIN);
    my $lightgray = $im->colorAllocate(217,217,217);
    my $white = $im->colorAllocate(255,255,255);
    my $black = $im->colorAllocate(0,0,0);
    my $gray  = $im->colorAllocate(127,127,127);
    my $red   = $im->colorAllocate(255,0,0);
    my $blue  = $im->colorAllocate(0,0,255);
    my $lightred = $im->colorAllocate(255,127,127);
    my $lightyellow = $im->colorAllocate(255,255,127);
    my $yellow = $im->colorAllocate(255,255,0);
    my $seq1color = $im->colorAllocate(@seq1RGB);
    my $seq2color = $im->colorAllocate(@seq2RGB);


    # draw position box

    my ($s, $e) = ($self->start_at(),$self->end_at());
    $im->filledRectangle($s * $factor+IMAGE_MARGIN, 0, #$OFFSET-10,
			 $e * $factor+IMAGE_MARGIN, $OFFSET + 6*3 + IMAGE_MARGIN,#$OFFSET+25,
			 $white);

    # draw exons
    my @color = ($white, $yellow, $lightred);
    foreach my $exontype(1,2)  {

	my $pExonPdl =
	    $self->pdlindex->slice(':,0')->where
		($self->pdlindex->slice(':,3') >= $exontype);
	$pExonPdl = $pExonPdl->where($pExonPdl>0);
	if ($pExonPdl->getdim(0)>3) {
	    print STDERR "pExonPdl $pExonPdl\n";
	    my $pExonStarts =
		$pExonPdl->slice('1:-1')->where(
	     	($pExonPdl->slice('1:-1' )- $pExonPdl->slice('0:-2')>1));
	    print STDERR "pExonStarts $pExonStarts\n";
	    my $pExonEnds =
		$pExonPdl->slice('0:-2')->where
		    (($pExonPdl->slice('1:-1' )- $pExonPdl->slice('0:-2')>1));
	    print STDERR "pExonEnds $pExonEnds\n";
	    my @exon_starts = ((list($factor*$pExonPdl->slice(0) +IMAGE_MARGIN)),
			       (list ($pExonStarts*$factor +IMAGE_MARGIN)));
	    print STDERR "ExonStarts @exon_starts\n";
	    my @exon_ends = ((list ($pExonEnds*$factor +IMAGE_MARGIN)),
			     (list($factor*$pExonPdl->slice(-1) +IMAGE_MARGIN)));
	    print STDERR "ExonEnds @exon_ends\n";

	    while ( my $exonstart = splice (@exon_starts, 0, 1))  {
		my $exonend = splice (@exon_ends, 0,1);
		$im->filledRectangle($exonstart, $OFFSET-1,
				     $exonend, $OFFSET+16,
				     $color[$exontype]);
	    }
	}
    }
    # draw seq1

    my $pSeq1Pdl =
	$self->pdlindex->slice(':,0')->where
	    ($self->pdlindex->slice(':,1') > 0);
    $pSeq1Pdl = $pSeq1Pdl->where($pSeq1Pdl>0);
    print STDERR "pSeq1Pdl $pSeq1Pdl\n";
    my $pSeq1Starts =
	$pSeq1Pdl->slice('1:-1')->where(
	    ($pSeq1Pdl->slice('1:-1' )- $pSeq1Pdl->slice('0:-2')>1));
    print STDERR "pSeq1Starts $pSeq1Starts\n";
    my $pSeq1Ends =
	$pSeq1Pdl->slice('0:-2')->where(
	    ($pSeq1Pdl->slice('1:-1' )- $pSeq1Pdl->slice('0:-2')>1));
    print STDERR "pSeq1Ends $pSeq1Ends\n";
    my @seq1_starts = (list($pSeq1Pdl->slice(0)*$factor +IMAGE_MARGIN),
		       (list ($pSeq1Starts*$factor +IMAGE_MARGIN)));
    print STDERR "Seq1starts @seq1_starts\n";
    my @seq1_ends = ((list ($pSeq1Ends*$factor +IMAGE_MARGIN),
		      list($pSeq1Pdl->slice(-1)*$factor +IMAGE_MARGIN)));
    print STDERR "Seq1starts @seq1_starts\n";
    print STDERR "Seq1ends @seq1_ends\n";

    # numbers
    my ($abs1_start, $abs1_end) =
	($self->_abs_pos(1,1), $self->_abs_pos($self->seq1length, 1));
    $im->string(gdTinyFont, $seq1_starts[0] - gdTinyFont->width*length($abs1_start)/2,
		$OFFSET - (gdTinyFont->height),
		$abs1_start, $black);
    $im->string(gdTinyFont,
		$seq1_ends[-1] - gdTinyFont->width*length($abs1_end)/2,
		$OFFSET - (gdTinyFont ->height),
		$abs1_end, $black);
    while ( my $seq1start = splice (@seq1_starts, 0, 1))  {
	my $seq1end = splice (@seq1_ends, 0,1);
	$im->filledRectangle($seq1start, $OFFSET+2, $seq1end, $OFFSET+5, $seq1color);
    }


    # draw seq2
    my $pSeq2Pdl =
	$self->pdlindex->slice(':,0')->where
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


    my ($abs2_start, $abs2_end) =
	($self->_abs_pos(1,2), $self->_abs_pos($self->seq2length, 2));
    $im->string(gdTinyFont, $seq2_starts[0] - gdTinyFont->width*length($abs2_start)/2,
		$OFFSET+15,
		$abs2_start, $black);
    $im->string(gdTinyFont,
		$seq2_ends[-1]- gdTinyFont->width*length($abs2_end)/2,
		$OFFSET+15,
		$abs2_end, $black);
    while ( my $seq2start = splice (@seq2_starts, 0, 1))  {
	my $seq2end = splice (@seq2_ends, 0,1);
	$im->filledRectangle($seq2start, $OFFSET+10, $seq2end,
			     $OFFSET+13, $seq2color);
    }

    # draw labels for the sequences
    $im->string(gdTinyFont,
		IMAGE_MARGIN -
		     (length($self->seq1name())+1)* gdTinyFont->width,
		$OFFSET-1,
		$self->seq1name(), $black);

    $im->string(gdTinyFont,
		IMAGE_MARGIN -
		     (length($self->seq2name())+1) * gdTinyFont->width,
		$OFFSET+8,
		$self->seq2name(), $black);



    return $im;


}


sub _drawgene_png  {
    my ($self, $which)  = @_;
    my $image_map_areas = "";
    my $FONT = gdSmallFont;
    my ($s, $e, $BEGIN_AT, $END_AT, $ref_site_iterator, $other_site_iterator,
	$ref_seq_color, $other_seq_color, $ref_seq_yoffset, $other_seq_yoffset);
    my $SPACING = 1; #vertical spacing between labels

    # set parameters and draw

    #$self->fsiteset1->sort_reverse();
    #my $iterator 1 = $self->fsiteset1 (-sort_by =>'start', -reverse =>1);
    $BEGIN_AT = $self->higher_pdlindex($self->start_at(), 0=>1);
    $END_AT   = $self->lower_pdlindex($self->end_at(),   0=>1);
    ($s, $e) = ($self->start_at(), $self->end_at());
    $ref_site_iterator = $self->fsiteset1->Iterator(-sort_by => 'start',
						    -reverse =>1);
    #$other_site_iterator = $self->fsiteset2->Iterator(-sort_by => 'start',
#						      -reverse =>1);
    $which = 1;

    my $factor = (IMAGE_WIDTH - 2*IMAGE_MARGIN)
	                        /
		 ($END_AT-$BEGIN_AT +1);
    my @labelrows;
    my ($leftmostpos, $leftmostrow);
    my @leftmostpos_in_row; # for labels
    while (my $site = $ref_site_iterator->next())  {

	next if $site->end > $END_AT;
	last if $site->start< $BEGIN_AT;
	my $xpos = ($site->start() - $BEGIN_AT) * $factor + IMAGE_MARGIN;
	my $row=0;
	while ($labelrows[$row]
	       and
	       (($xpos + (length($site->Matrix->{name}) +1) * $FONT->width)
	       >
	       $leftmostpos )
	       and
	       (($xpos + (length($site->Matrix->{name}) +1) * $FONT->width)
	       >
	       $leftmostpos_in_row[$row]))

	{
	    $row++;
	}
	for my $i (0..$row) {$leftmostpos_in_row[$i] = $xpos;}
	($leftmostpos, $leftmostrow) = ($xpos, $row);
	unshift @{$labelrows[$row]}, {pos => $xpos, site => $site} ;
    }

    # create image

    my $OFFSET =
	IMAGE_MARGIN/4
	+(scalar(@labelrows)+2)*($FONT->height+$SPACING);
    my $im = GD::Image->new(IMAGE_WIDTH,
			    $OFFSET + 6*3 + IMAGE_MARGIN/2);
    my $white = $im->colorAllocate(255,255,255);
    my $black = $im->colorAllocate(0,0,0);
    my $gray  = $im->colorAllocate(127,127,127);
    my $red   = $im->colorAllocate(255,0,0);
    my $blue  = $im->colorAllocate(0,0,255);
    my $yellow      = $im->colorAllocate(200,200,0);
    my $lightyellow = $im->colorAllocate(255,255,127);
    my $lightred = $im->colorAllocate(255,127,127);
    my $seq1color = $im->colorAllocate(@seq1RGB);

    $ref_seq_yoffset = $OFFSET+5;
    $ref_seq_color = $seq1color;

    # alignment positions of $self->start and $self->end

    ($s, $e) = ($self->pdlindex($BEGIN_AT, $which=>0),
		   $self->pdlindex($END_AT,   $which=>0));

    # reference sequence (continuous)
    $im->filledRectangle(IMAGE_MARGIN, $ref_seq_yoffset,
			 IMAGE_WIDTH - IMAGE_MARGIN -1, $ref_seq_yoffset +3, $ref_seq_color);
    # put the labels

    foreach my $row(0..$#labelrows)  {
	# line corner y value:
	my $lcy = $OFFSET - ($row+1.5)*($FONT->height()+$SPACING);
	# label y coordinate:
	my $laby = $OFFSET - ($row+2)*($FONT->height()+$SPACING);
	foreach my $label (@{$labelrows[$row]})  {

	    $im->line($label->{pos}, $ref_seq_yoffset-1,
		      $label->{pos}, $lcy,
		      $gray);
	    $im->line($label->{pos}, $lcy,
		      $label->{pos} + $FONT->width()/2-1, $lcy,
		      $gray);
	    my $name = $label->{site}->Matrix->{name};
	    $im->string($FONT, $label->{pos} + $FONT->width() -1, $laby,
			$name, $blue);
	    $image_map_areas .= $self->maparea_factor
		( $which, $label->{site}, $label->{site2},
		  $label->{pos}+$FONT->width() -1,
		  $laby,
		  $label->{pos}+$FONT->width()*(length($name)+1)-1,
		  $laby+$FONT->height());
	}
    }
    # draw exons

    print STDERR $self->pdlindex->slice("$s:$e,$which");
    foreach my $exontype (1,2)  {
	my $pExonPdl =
	    $self->pdlindex->slice("$s:$e,$which")->where
		($self->pdlindex->slice("$s:$e,3") >= $exontype);
	$pExonPdl = $pExonPdl->where($pExonPdl>0);
	if ($pExonPdl->getdim(0) > 3) {
	    my $pExonStarts =
		$pExonPdl->slice('1:-1')->where
		    (($pExonPdl->slice('1:-1' )- $pExonPdl->slice('0:-2')>1));
	    my $pExonEnds =
		$pExonPdl->slice('0:-2')->where
		    (($pExonPdl->slice('1:-1' )- $pExonPdl->slice('0:-2')>1));
	    my @exon_starts =
		((list($factor*($pExonPdl->slice(0)-$BEGIN_AT) +IMAGE_MARGIN)),
		 (list (($pExonStarts-$BEGIN_AT)*$factor +IMAGE_MARGIN)));
	    my @exon_ends = ((list (($pExonEnds-$BEGIN_AT)*$factor +IMAGE_MARGIN)),
			     (list($factor*($pExonPdl->slice(-1)-$BEGIN_AT) +IMAGE_MARGIN)));

	    while ( my $exonstart = splice (@exon_starts, 0, 1))  {
		my $exonend = splice (@exon_ends, 0,1);
		$im->filledRectangle($exonstart, $OFFSET, $exonend, $OFFSET+11,
				 ($lightyellow, $lightred)[$exontype-1]);
		$im->rectangle($exonstart-1, $OFFSET-1, $exonend+1, $OFFSET+12,
				 $seq1color) if $exontype == 1;
	    }
	}

    }


    return ($im, "<MAP name='seq".$which."map'> $image_map_areas </MAP>");

}


sub _parse_seq {
    my ($self) = @_;

    my ($seq1, $start);
    $seq1 = $self->{job}->seq1->seq();

    $self->alignseq1($seq1);

    $self->alignseq1array( my @seq1 = ("-", split('', $seq1) ));

    my (@seq1index, @seq2index);

    $self->pdlindex( pdl [ [list sequence($#seq1+1)],
			    [list sequence($#seq1+1)], #[@seq1index],
			   [list zeroes ($#seq1+1)],
			   [list zeroes ($#seq1+1)]
			 ]) ;

    # mark exons in row 3
    my $seqobj =$self->{job}->seq1();
    $self->seq1name($self->{job}->seq1->id());
    #$self->seq3name($self->{job}->seq3->id()) if $self->{job}->seq3();

    my @features = $seqobj->top_SeqFeatures;
    for my $exon (@features)  {
	next unless ref($exon) =~ /::Exon/;
	my $exonstart = $self->pdlindex($exon->start(), 1=>0);
	my $exonend   = $self->pdlindex($exon->end(),   1=>0);

	my $sl = $self->pdlindex->slice("$exonstart:$exonend,3");
	$sl++;
    }

    # mark ORF
    my @cdnafeatures = ();
    if ($self->{job}->seq3())  {
	@cdnafeatures = $self->{job}->seq3()->top_SeqFeatures;
    }
    my ($cdnaexonstart, $cdnaexonend) = (0,0);
    my $orf;
    for my $feat (@cdnafeatures)  {
	$orf = $feat
	    if (ref($feat) eq "Bio::SeqFeature::Generic"
		and $feat->primary_tag eq 'orf');
	next unless ref($feat) =~ /::Exon/;
	$cdnaexonstart = $feat->start() if $cdnaexonstart > $feat->start();
	$cdnaexonend   = $feat->end()   if $cdnaexonend   < $feat->end();
    }
    if ($orf)  {
	my ($orf_start, $orf_end, $orf_strand) =
	    ($orf->start(), $orf->end(), $orf->strand());
	my $cdnaPdl =
	    $self->pdlindex->slice(':,3')->where
		($self->pdlindex->slice(':,3')==1);
	print STDERR "\ncdnaPdl:".$cdnaPdl."\n";
	my $orfPdl;

	$orf_end = $cdnaexonend if $orf_end>$cdnaexonend;
	$orf_start = $cdnaexonstart if $orf_start < $cdnaexonstart;

	if ($orf_strand eq "1") {
	    $orf_end = $cdnaexonend if $orf_end>$cdnaexonend;
	    $orfPdl = $cdnaPdl->slice
		(($orf_start-$cdnaexonstart).":".($orf_end-$cdnaexonstart-1));
	}
	else  {
	    $orfPdl = $cdnaPdl->slice
		(($self->{job}->seq3->length() - $orf_end + $cdnaexonstart).
		 ":".
		 ($self->{job}->seq3->length()-$orf_start+$cdnaexonstart-1));
	}
	print STDERR "\norfPdl:".sprintf($orfPdl)."\n";
	$orfPdl++; # set it to 2
    }

    print STDERR "\nROW#3:".$self->pdlindex->slice(':,3');
    return 1;

}





sub _do_sitesearch  {
     my ($self, $MATRIXSET, $THRESHOLD, $CUTOFF) = @_;

     print STDERR "SEARCHING 1...\n";
     my $seqobj1 = Bio::Seq->new(-seq=>$self->alignseq1(),
				-id => "Seq1");
     my $siteset1 =
	$MATRIXSET->search_seq(-seqobj => $seqobj1,
			    -threshold => $THRESHOLD);
     #$siteset1->sort();
     my $iterator1 = $siteset1->Iterator (-sort_by => "start");

    # initialize filtered sitesets:

    $self->fsiteset1 (TFBS::SiteSet->new());

    while (my $site1 = $iterator1->next()) {

	if (	    # exclude ORF test
	    (!($self->{job}->exclude_orf()
	       and $self->pdlindex($site1->start(), 1=>3) == 2
	       )
	     )
	    )
	{
	    $self->fsiteset1->add_site($site1);
	}
    }
}



1;


