# AT::GFX::MulChr module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD


package AT::GFX::MulSynChr;

use strict;
use vars '@ISA';
use Carp;
use AT::GFX::MulChr;

@ISA = qw(AT::GFX::MulChr);


sub new
{
    my ($caller, %args) = @_;

    my $self = $caller->SUPER::new(%args);

    $self->{syn_chr_width} = $args{syn_chr_width} || $self->chr_width;
    $self->{proximity_group_threshold} = $args{proximity_group_threshold}
	|| 50000;

    # increase bottom margin so we can fit a legend there
    $self->{bottom_margin} = $self->{bottom_margin} + 20;

    return $self;
}


sub draw {
    my ($self, $gendb, $ref_chr_list, $ref_tracks,
	$syn_file, $syn_chr_list, $syn_tracks) = @_;

    # Create canvas
    $self->_new_img(scalar @$ref_chr_list);

    # Read synteny data
    my $syn_data = $self->_read_syn_file($syn_file);

    # Allocate colors for secondary ('syntenic') chromosomes
    $self->_add_syn_chr_cols($syn_chr_list);

    # Create features representing secondary chromosomes along
    # reference chromosomes
    my @syn_bands;
    foreach my $syn_part (@$syn_data) {
	push @syn_bands,
	    { chr => $syn_part->{chrA},
	      start => $syn_part->{startA},
	      end => $syn_part->{endA},
	      color => 'synchr'.$syn_part->{chrB},
	      id => $syn_part->{chrA}.'_'.$syn_part->{startA} };
    }
    my $syn_chr_track = {
	side => 1,
	feat_width => $self->syn_chr_width,
	feat_bump_threshold => -1,
	feat_min_span => 1,
	features => \@syn_bands
    };

    # Recalculate positions for features on the secondary chromosomes
    # so the features can be drawn along the reference chromosomes
    my @new_syn_tracks;
    foreach my $track (@$syn_tracks) {
	my @new_feats;
	my $i = 0;
	foreach my $feat (sort { $a->{chr} cmp $b->{chr} or
			         $a->{start} <=> $b->{start} }
			  @{$track->{features}})
	{
	    my $chr = $feat->{chr};
	    my $start = $feat->{start};
	    my $end = $feat->{end};
	    #print "Feat: $chr:$start-$end\n";
	    #print "Synt: ",$syn_data->[$i]->{chrB},":",
		#        $syn_data->[$i]->{startB},"-",
		#	$syn_data->[$i]->{endB},"\n"
		#if($i < @$syn_data);
	    while($i < @$syn_data and
		  ($syn_data->[$i]->{chrB} cmp $chr or
		   $syn_data->[$i]->{endB} <=> $start) < 0) {
		 $i++;
		 #print "Synt: ".$syn_data->[$i]->{chrB}.":",
			#        $syn_data->[$i]->{startB}."-",
			#	$syn_data->[$i]->{endB},"\n"
			#if($i < @$syn_data);				
	    }
	    my @my_new_feats;
	    for(my $j = $i; $j < @$syn_data and
			    $syn_data->[$j]->{chrB} eq $chr and
			    $syn_data->[$j]->{startB} <= $end; $j++)
	    {
		my $syn = $syn_data->[$j];
		if($syn->{endB} >= $start) {
		    my $k = ($syn->{endA} - $syn->{startA}) /
			    ($syn->{endB} - $syn->{startB});
		    my ($new_start, $new_end);
		    if($syn->{strand} eq '+') {
			$new_start =
			 int(($start - $syn->{startB}) * $k + $syn->{startA} + 0.5);
			$new_end =
			int(($end - $syn->{startB}) * $k + $syn->{startA} + 0.5);
		    }
		    else {
			$new_start =
			 int(-($end - $syn->{startB}) * $k + $syn->{endA} + 0.5);
			$new_end =
			 int(-($start - $syn->{startB}) * $k + $syn->{endA} + 0.5);
		    }
		    $new_start = $syn->{startA} if ($new_start < $syn->{startA});
		    $new_end = $syn->{endA} if ($new_end > $syn->{endA});
		    #  ^^ make sure rounding errors do not make feature
		    #     stick outside the synteny block
		    push @my_new_feats, {%$feat,
				      chr => $syn->{chrA},
				      start => $new_start,
				      end => $new_end,
				      original_feature => $feat };
#print "FEAT\t",$syn->{chrA},"\t$new_start\t$new_end\t",$new_end-$new_start,
#	"\t$chr\t$start\t$end\t",$end-$start,"\n"
#	if($syn->{chrA} == 11 and $chr == 17 and $new_end <= 77497466);
		}
	    }
	    if(@my_new_feats) {
		$self->_merge_close_features(\@my_new_feats);
		push @new_feats, @my_new_feats;
	    }
	    else {
		warn "Could not recalculate feature postion $chr:$start-$end (no synteny data)\n";
	    }
	}
	push @new_syn_tracks, { %$track,
				features => \@new_feats }
    }
    

    # Draw image
    my $tracks = [ @$ref_tracks, $syn_chr_track, @new_syn_tracks ];
    $self->_draw($gendb, $ref_chr_list, $tracks);

    # Add a legend describing secondary chromosome colors
    $self->_draw_synchr_legend($syn_chr_list);
}


sub _merge_close_features
{
    my ($self, $feats) = @_;

    # Group features; features that are sufficiently close
    # are put in the same group
    my @feat_groups;
    my $i = 0;
    my $prev_f;
    foreach my $f (sort { $a->{chr} cmp $b->{chr} or
			  $a->{start} <=> $b->{start} } @$feats) {
	if($prev_f) {
	    my $d = $f->{start} - $prev_f->{end};
	    $i++ if($f->{chr} ne $prev_f->{chr} or
		    $d > $self->proximity_group_threshold);
	}
	push @{$feat_groups[$i]}, $f;
	$prev_f = $f;
    }

    # Convert multi-feature groups to multi-part features
    @$feats = ();
    foreach my $fg (@feat_groups) {
	if(@$fg == 1) {
	    push @$feats, $fg->[0];
	}
	else {
	    push @$feats, { %{$fg->[0]},
			    end => $fg->[-1]->{end},
			    subfeatures => $fg };
	}
    }
}


sub _draw_synchr_legend
{
# we could make this method call MulChr::_draw_legend()
    my ($self, $chr_list, $chr_cols) = @_;
    
    my $img = $self->img;
    my $x = $self->left_margin;
    my $y = $self->height - $self->bottom_margin;
    my $i = 0;
    foreach my $chr (@$chr_list) {
	my $col = $self->colors->{"synchr$chr"};
	$img->filledRectangle($x, $y, $x+15, $y+4, $col);
	$img->string(gdMediumBoldFont, $x, $y+9, $chr, $self->colors->{black});
	$x += 25;
	$i++;
    }
}


sub _read_syn_file {
    my ($self, $syn_file) = @_;

    # Read file into array
    open SYN_FILE, $syn_file or croak("Could not open file $syn_file");
    my @syn_data;
    while(my $line = <SYN_FILE>) {
	chomp $line;
	my @fields = split /\t/, $line;
	push @syn_data,
	    {
		chrA => ($fields[1] =~ /chr(.*)/)[0],
		startA => $fields[2],
		endA => $fields[3],
		chrB => ($fields[4] =~ /chr(.*)/)[0],
		startB => $fields[7],
		endB => $fields[8],
		strand => $fields[6]
	    };
    }
    close SYN_FILE;

    # Make sure array is sorted by chrB, startB
    #  (since will do lookups from secondary chromosome positions)
    @syn_data = sort { $a->{chrB} cmp $b->{chrB} or
		       $a->{startB} <=> $b->{startB} } @syn_data;

    # Return ref to array
    return \@syn_data;
}


sub _add_syn_chr_cols
{
    my ($self, $chr_list) = @_;

    # read color codes
    my @cols;
    while (<DATA>) {
	chomp;
	last if /^__END__/;
	my ($r,$g,$b) = split /\s+/;
	push @cols, [hex $r,hex $g,hex $b];
    }

    # allocate colors
    my $i = 0;
    foreach my $chr (@$chr_list) {
	my $col;
	if($i < @cols) {
	    $col = $cols[$i];
	}
	else {
	    warn "Out of colors for syntenic chromosomes.";
	    $col = [0x00, 0x00, 0x00];
	}
	$self->_add_color("synchr$chr", @$col);
	$i++;
    }
   
}

1;


__DATA__
C0	00	00
FF	00	00
FF	B0	B0
60	90	60
00	F0	00
B0	FF	B0
60	60	A0
00	00	FF
B0	B0	FF
90	90	60
C0	C0	00
FF	FF	00
A0	60	A0
FF	00	FF
FF	B0	FF
00	B0	B0
00	FF	FF
B0	FF	FF
B2	40	00
B8      86      0B
BD      B7      6B
00	00	00
90	90	90
D0	D0	D0
__END__
