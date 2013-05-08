# AT::GFX::Locus module
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



package AT::GFX::Locus;

use strict;
use vars '@ISA';
use DBI;
use Carp;
use AT::Root;
use AT::DB::Binner;
use Bio::Graphics::Panel;
use Bio::SeqFeature::Generic;
use Benchmark;


@ISA = qw(AT::Root);


sub new
{
    my ($caller, %args) = @_;

    my $mapping_db = $args{track_db} || die "No track db";
    my $organism = $args{organism} || die "No organism";

    my $self = bless {
	organism => $organism,
	track_db => $args{track_db},
	tracks => ($args{tracks} || die "No tracks"),
	width => ($args{width} || 1250),
	pad_left => ($args{pad_left} || 30),
	pad_right => ($args{pad_right} || 220),
	pad_top => ($args{pad_top} || 10),
	pad_bottom => ($args{pad_bottom} || 10),
	flip => ($args{flip} || 0),
	min_exon_size => ($args{min_exon_size} || 0),
	svg => $args{svg},
	_binner => AT::DB::Binner->new(),
	_panel => undef,
	_errstr => ''
    }, ref $caller || $caller;

    return $self;   
}


sub errstr {
    return shift->{_errstr};
}


sub panel
{
    return shift->{_panel};
}


sub gd_image
{
    my $panel = shift->{_panel};
    return $panel ? $panel->gd : undef;
}


#sub svg
#{
#    my $panel = shift->{_panel};
#    return $panel ? $panel->svg : undef;
#}


sub draw_locus
{
    my ($self, $chr, $start, $end) = @_;

    my $organism = $self->organism || "";

    my @g;

    my $panel = Bio::Graphics::Panel->new
         ( -length => $end-$start + 1,
	   -offset => $start,
	   -width => $self->{width},
	   #-grid => 1,
	   -pad_left => $self->{pad_left},
	   -pad_right => $self->{pad_right},
	   -pad_top => $self->{pad_top},
	   -pad_bottom => $self->{pad_bottom},
	   -flip => $self->{flip},
	   -key_style => 'between',
	   -image_class => $self->{svg} ? 'GD::SVG' : 'GD'
    );

    my $full_length = Bio::SeqFeature::Generic->new(-start=> $panel->start,
						-end=> $panel->end,
						-primary=>'ruler');
    my $title = "$organism $chr:$start-$end:".($self->{flip}?'-':'+');
    $panel->add_track('arrow',
		      $full_length,
		      -tick => 2,
		      -fgcolor => 'black',
		      -double => 1,
		      -key => $title);

    foreach my $track (@{$self->tracks}) {
	$self->_add_track($panel, $chr, $start, $end, $track);
    }

    $self->{_panel} = $panel;
    return 1;
}


sub _add_track
{
    my ($self, $panel, $chr, $start, $end, $track) = @_;

    my $type = $track->{type};
    my $name = $track->{name};  # display name
    my $glyph = $track->{glyph};     # Bio::Graphics glyph
    my $bgcolor = $track->{bgcolor};     # color
    my $fgcolor = $track->{fgcolor};     # color

    my $features;
    if($type eq 'table') { 
	my $tables = $track->{tables};
	my ($table1,$table2) = @$tables;
	my $format = $track->{format};
	if($format =~ s/Split$//) {
	    $table1 = "${chr}_${table1}";
	}
	if($format eq 'simple') {
	    $features = $self->_get_simple_features($chr, $start, $end, $table1, 0, 0);
	}
	elsif($format eq 'simpleBinned') {
	    $features = $self->_get_simple_features($chr, $start, $end, $table1, 1, 0);
	}
	elsif($format eq 'gene') {
	    $features = $self->_get_gene_features($chr, $start, $end, $table1);
	}
	elsif($format eq 'refGene') {
	    $features = $self->_get_refGene_features($chr, $start, $end, $table1);
	}
	elsif($format eq 'knownGene') {
	    $features = $self->_get_knownGene_features($chr, $start, $end, $table1, $table2);
	}
	elsif($format eq 'flyBaseGene') {
	    $features = $self->_get_flyBaseGene_features($chr, $start, $end, $table1, $table2);
	}
	elsif($format eq 'psl') {
	    $features = $self->_get_psl_features($chr, $start, $end, $table1, 0);
	}
	elsif($format eq 'pslProtToDna') {
	    $features = $self->_get_psl_features($chr, $start, $end, $table1, 1);
	}
	elsif($format eq 'rmsk') {
	    $features = $self->_get_rmsk_features($chr, $start, $end, $table1);
	}
	else {
	    croak "Unknown format $format";
	}
    }
    elsif($type eq 'memory') {
	my $data = $track->{data}; # data should be in format [name, chr, strand, [[start,end]]]
	$features = $self->_get_features_from_array($chr, $start, $end, $data);
    }
    elsif($type eq 'file') {
	my $fn = $track->{filename};
	my $format = $track->{format};
	if($format eq 'bed') {
	    $features = $self->_get_bed_file_features($chr, $start, $end, $fn);
	}
	else {
	    croak "Only bed format supported for annotation files";
	}
    }
    else {
	croak "Unknown track type $type";
    }

    # draw track with features
    $panel->add_track
	    ( $features || [],
	      -label => 1,
	     -glyph => $glyph,
	     -bgcolor => $bgcolor,
	     -fgcolor => $fgcolor,
	     -key => $name,
	    );
}


sub timediff_str
{
    my ($b1, $b2) = @_;
    my $t = timestr(timediff($b2, $b1));
    return "[$t sec]";
}


sub _get_simple_features {
    my ($self, $chr, $start, $end, $table, $is_binned) = @_;

    my $db = $self->track_db or croak 'No mapping db';
    my @locations = $db->get_simple_feature_locations_in_region(chr => $chr,
								start => $start,
								end => $end,
								table => $table,
								is_binned => $is_binned);
    my @feats = map
	{ Bio::SeqFeature::Generic->new(-start => $_->[0], -end => $_->[1]) }
	@locations;

    return \@feats;
}
    

sub _get_features_from_array
{
    my ($self, $view_chr, $view_start, $view_end, $all_elements) = @_;

    my @feats;
    foreach my $element (@$all_elements) {
	my ($id, $chr, $strand, $pos) = @$element;
	next unless ($chr eq $view_chr and $pos->[0][0] <= $view_end and $pos->[-1][-1] >= $view_start);
	my $main = Bio::SeqFeature::Generic->new(-start => $pos->[0][0],
						 -end => $pos->[-1][1],
						 -strand => $strand);
	foreach my $iv (@$pos) {
	    my $subfeature = Bio::SeqFeature::Generic->new(-start => $iv->[0],
							   -end => $iv->[1],
							   -strand => $strand);
	    $main->add_SeqFeature($subfeature);
	}
	push @feats, $main;
    }

    return \@feats;
}


sub _get_gene_features {
    my ($self, $chr, $start, $end, $table) = @_;

    my $db = $self->track_db or croak 'No mapping db';
    my @mappings = $db->get_features_in_region(chr => $chr,
					       start => $start,
					       end => $end,
					       types => [$table]);

    my @feats;
    foreach my $m (@mappings) {
	my $strand = $m->strand_numeric;
	my $main = Bio::SeqFeature::Generic->new(-display_name => $m->qName,
						 -start => $m->tStart,
						 -end => $m->tEnd,
						 -strand => $strand);
	foreach my $iv ($m->abs_HSP_tPos) {
	    my $subfeature = Bio::SeqFeature::Generic->new(-start => $iv->[0],
							   -end => $iv->[1],
							   -strand => $strand);
	    $main->add_SeqFeature($subfeature);
	}
	push @feats, $main;
    }

    return \@feats;
}


sub _get_refGene_features {
    my ($self, $chr, $start, $end, $table) = @_;

    my @feats;
    my $db = $self->track_db or croak 'No mapping db';
    my $bin_string = $self->{_binner}->bin_restriction_string($start,$end);
    my $query = "SELECT distinct strand, exonStarts, exonEnds, name2, name FROM $table
                 WHERE $bin_string and chrom = ? AND txEnd >= ? AND txStart < ?";
    my $sth = $db->dbh->prepare($query) or die "could not prepare query [$query]";
    $sth->execute($chr,$start,$end);
    return $self->_get_gene_features_from_sth($sth);
}


sub _get_knownGene_features {
    my ($self, $chr, $start, $end, $table1, $table2) = @_;

    my $db = $self->track_db or croak 'No mapping db';
    croak "Need two tables to make a knownGene track" unless($table1 and $table2); 
    # Note: known genes are not binned
    my $query = "SELECT distinct a.strand, a.exonStarts, a.exonEnds, b.geneSymbol, name FROM $table1 a, $table2 b
                 WHERE a.name = b.kgID and a.chrom = ? AND a.txEnd >= ? AND a.txStart < ?";
    my $sth = $db->dbh->prepare($query) or die "could not prepare query [$query]";
    $sth->execute($chr,$start,$end);
    return $self->_get_gene_features_from_sth($sth);
}


sub _get_flyBaseGene_features {
    my ($self, $chr, $start, $end, $table1, $table2) = @_;

    my $db = $self->track_db or croak 'No mapping db';
    croak "Need two tables to make a flyBaseGene track" unless($table1 and $table2); 
    my $query = "SELECT distinct a.strand, a.exonStarts, a.exonEnds, b.symbol, a.name FROM $table1 a, $table2 b
                 WHERE a.name = b.name and a.chrom = ? AND a.txEnd >= ? AND a.txStart < ?";
    my $sth = $db->dbh->prepare($query) or die "could not prepare query [$query]";
    $sth->execute($chr,$start,$end);
    return $self->_get_gene_features_from_sth($sth);
}


sub _get_gene_features_from_sth 
{
    my ($self, $sth) = @_;

    my @feats;
    while(my ($strand,$exonStart_str,$exonEnd_str,$symbol,$accession) = $sth->fetchrow_array) {
	my @exonStarts = map { $_ + 1 } split /,/, $exonStart_str;
	my @exonEnds = split /,/, $exonEnd_str;
	$strand = $strand eq '+' ? 1 : -1;
	my $main = Bio::SeqFeature::Generic->new(-display_name => $symbol || $accession,
						 -start => $exonStarts[0],
						 -end => $exonEnds[-1],
						 -strand => $strand);
	while(@exonStarts) {
	    my $subfeature = Bio::SeqFeature::Generic->new(-start => shift @exonStarts,
							   -end => shift @exonEnds,
							   -strand => $strand);
	    $main->add_SeqFeature($subfeature);
	}
	push @feats, $main;
    }

    return \@feats;
}


sub _get_psl_features {
    my ($self, $chr, $start, $end, $table, $query_is_prot) = @_;

    my $db = $self->track_db or croak 'No mapping db';
    my $bin_string = $self->{_binner}->bin_restriction_string($start,$end);
    my $query = "SELECT qName, strand, tSize, tStarts, blockSizes FROM $table
                 WHERE $bin_string and tName = ? AND tEnd >= ? AND tStart < ?";
    my $sth = $db->dbh->prepare($query) or die "could not prepare query [$query]";
    $sth->execute($chr,$start,$end);

    my @feats;
    while(my ($name, $strand, $chr_size, $blockStart_str, $blockSize_str) = $sth->fetchrow_array) {
	my @blockStarts = split /,/, $blockStart_str;
	my @blockSizes = split /,/, $blockSize_str;
	@blockSizes = map { $_ * 3 } @blockSizes if($query_is_prot);
	my ($qStrand, $tStrand) = split(//, $strand);
	$qStrand = $qStrand eq '-' ? -1 : 1;
	$tStrand = $tStrand eq '-' ? -1 : 1;
        # Note: in the above calculations, it is important that we default to 1, since
        # tStrand may be an empty string and should then be interpreted as '+'.
	$strand = $qStrand * $tStrand;
	if($tStrand == -1) {
	    @blockStarts = reverse @blockStarts;
	    @blockSizes = reverse @blockSizes;
	    for my $i (0..@blockStarts-1) {
		$blockStarts[$i] = $chr_size - $blockStarts[$i] - $blockSizes[$i]; 
	    }
	}
	my $ftStart = $blockStarts[0]+1;
	my $ftEnd = $blockStarts[-1]+$blockSizes[-1];
	my $main = Bio::SeqFeature::Generic->new(-display_name => $name,
						 -start => $ftStart,
						 -end => $ftEnd,
						 -strand => $strand);
	while(@blockStarts) {
	    my $blockStart = shift @blockStarts;
	    my $blockSize = shift @blockSizes;
	    my $subfeature = Bio::SeqFeature::Generic->new(-start => $blockStart+1,
							   -end => $blockStart+$blockSize,
							   -strand => $strand);
	    #print STDERR join(" ", $name, $ftStart, $ftEnd, $blockStart+1, $blockStart+$blockSize), "\n";
	    $main->add_SeqFeature($subfeature);
	}
	push @feats, $main;
    }

    return \@feats;


}


sub _get_rmsk_features {
    my ($self, $chr, $start, $end, $table) = @_;

    my @feats;
    my $db = $self->track_db or croak 'No track db';
    my $bin_string = $self->{_binner}->bin_restriction_string($start,$end);
    my $query = "SELECT repName, strand, genoStart, genoEnd FROM $table
                 WHERE $bin_string and genoName = ? AND genoEnd >= ? AND genoStart < ?";
    my $sth = $db->dbh->prepare($query) or die "could not prepare query [$query]";
    $sth->execute($chr, $start, $end);
    while(my ($name, $strand, $start, $end) = $sth->fetchrow_array) {
	push @feats, Bio::SeqFeature::Generic->new(-display_name => "", #$name,
						   -start => $start,
						   -end => $end,
						   -strand => $strand);
    }

    return \@feats;
}


sub _get_bed_file_features 
{
    my ($self, $query_chr, $query_start, $query_end, $fn) = @_;

    my @feats;
    open IN, $fn or die "could not open $fn";
    while(my $line = <IN>) {
	chomp $line;
	next if($line =~/^browser/ or $line =~ /^track/ or $line =~ /^\s*$/);
	my ($chr, $start, $end, $name, undef, $strand, undef, undef, undef, undef, $blockStarts_str, $blockSizes_str) = split /\t/, $line;
	unless(defined($chr) and defined($start) and defined($end)) {
	    die "Error parsing file $fn, line [$line].\n";
	}
	$start++;
	next unless($chr eq $query_chr and $start <= $query_end and $end >= $query_start);
	if(defined $strand) {
	    if($strand eq '+') {
		$strand = 1;
	    }
	    elsif($strand eq '-') {
		$strand = -1;
	    }
	    else {
		die "Invalid strand value $strand in file $fn";
	    }
	}
	else {
	    $strand = 1;
	}
	my $main = Bio::SeqFeature::Generic->new(-display_name => $name || '',
						 -start => $start,
						 -end => $end,
						 -strand => $strand);
	if(defined($blockStarts_str) and defined($blockSizes_str)) {
	    my @blockStarts = split /,/, $blockStarts_str;
	    my @blockSizes = split /,/, $blockSizes_str;
	    while(@blockStarts) {
		my $blockStart = shift @blockStarts;
		my $blockSize = shift @blockSizes;
		die "Invalid blockSizes $blockSizes_str in file $fn" unless($blockSize);
		my $subfeature = Bio::SeqFeature::Generic->new(-start => $blockStart+1,
							       -end => $blockStart+$blockSize,
							       -strand => $strand);
		$main->add_SeqFeature($subfeature);
	    }
	}

	print STDERR "added bed file feature ".$main->start."-".$main->end."\n";

	push @feats, $main;
    }
    close IN;

    return \@feats;
}


sub _get_abs_exon_coords
{
    my ($self, $start, $exon_str) = @_;
    my @e;
    $start--;
    foreach my $se (split /;/, $exon_str) {
	my ($s,$e) = split /,/, $se;
	push @e, [$s+$start, $e+$start];
	#print STDERR "* $s\t$e\n";
    }
    return \@e;
}


1;

