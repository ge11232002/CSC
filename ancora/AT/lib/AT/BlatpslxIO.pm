# AT (Alternative Transcripts) module for AT::BlatpslxIO
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::BlatpslxIO - read and write psl/pslx (blat output format)

=head1 SYNOPSIS

use AT::BlatpslxIO;

# Filter out mappings with less than 100 matching bases

 my $map_in = AT::BlatpslxIO->new(-fh => \*STDIN);
 my $map_out = AT::BlatpslxIO->new(-fh => \*STDOUT);

 while(my $m = $map_in->next_mapping) {
     $map_out->write_mapping($m)
 	if(($m->matches + $m->repMatches) >= 100);
 }


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut



package AT::BlatpslxIO;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use AT::Mapping;
use AT::HSP;

@ISA = qw(AT::Root);

=head2 new

 Title     : new
 Usage     : my $psl_in = AT::BlatpslxIO->new(-file => 'in.psl');
 Function  : Constructor.
 Returns   : AT::BlatpslxIO
 Args      : -fh	  	filehandle
	     -file		filename
	     -has_bin_field	optional; use when reading some
				UCSC database dumps
	     Give either -fh or -file. If -file is given,
	     the method will try to open the file with that name.
	     Prefix the filename with '>' if the file is to be
	     opened for output.

=cut

sub new {
    my ($caller, %args) = @_;

    my $fh = $args{'-fh'};
    unless ($fh) {
	unless (exists $args{'-file'})  {
	    $fh = \*STDIN;
	}
	else  {
	    open FILE, $args{'-file'} 
	    or croak "Could not open ".$args{'-file'};
            $fh = \*FILE;
	}
    }

    # The following option is for compability with UCSC database dumps.
    # Generally, the UCSC blat result tables follow psl format. However,
    # some contain an extra initial field (called 'bin').
    my $has_bin_field = $args{'-has_bin_field'} || 0;

    my $self = bless { _fh => $fh,
		       has_bin_field => $has_bin_field
		       }, ref $caller || $caller;
    return $self;
}

=head2 next_mapping

 Title     : next_mapping
 Usage     : my $m = $psl_in->next_mapping();
 Function  : Gets the next mapping in the stream
 Returns   : AT::Mapping
 Args      : -

=cut

sub next_mapping {
    my ($self) = @_;
    my $fh = $self->{'_fh'};
    my $blatline = <$fh> || return undef;
    return $self->next_mapping 
	if ($blatline =~ /^\s*match/ 
            or $blatline =~ /^\-\-\-/
            or $blatline =~ /^psLayout/
            or $blatline eq "\n")
        ;
    return $self->_blat_line_to_mapping($blatline);
    
}


=head2 write_mapping

 Title     : write_mapping
 Usage     : $psl_out->write_mapping($mapping);
 Function  : Writes a mapping to the stream
 Returns   : -
 Args      : AT::Mapping

=cut


sub write_mapping  {
    my ($self, $m) = @_;
    croak "write_mapping needs a mapping" unless(defined $m);
    my $fh = $self->{'_fh'};
    my ($blockSizes, $qStarts, $tStarts) = ('','','');
    foreach my $hsp ($m->all_HSPs) {
	$blockSizes .= $hsp->blockSize.',';
	$qStarts .= ($hsp->qStart-1).',';
	$tStarts .= ($hsp->tStart-1).',';
    }
    my $line = join "\t",
	($self->has_bin_field() ? $m->bin : (),
	 $m->matches,
	 $m->misMatches,
	 $m->repMatches,
	 $m->nCount,
	 $m->qNumInsert,
	 $m->qBaseInsert,
	 $m->tNumInsert,
	 $m->tBaseInsert,
	 $m->strand,
	 $m->qName,
	 $m->qSize,
	 $m->qStart-1,
	 $m->qEnd,
	 $m->tName,
	 $m->tSize,
	 $m->tStart-1,
	 $m->tEnd,
	 $m->blockCount,
	 $blockSizes,
	 $qStarts,
	 $tStarts);
    print $fh $line, "\n";
}


sub _blat_line_to_mapping  {
    my ($self, $blatline) = @_;
    chomp $blatline; $blatline =~ s/\s+$//;
    my @blatline_fields = split /\t/, $blatline;
    my $bin = ($self->{has_bin_field}) ? shift @blatline_fields : undef;
    my ($matches,
	$misMatches,
	$repMatches,
	$nCount,
	$qNumInsert,
	$qBaseInsert,
	$tNumInsert,
	$tBaseInsert,
	$strand,
	$qName,
	$qSize,
	$qStart,
	$qEnd,
	$tName,
	$tSize,
	$tStart,
	$tEnd,
	$blockCount,
	$blockSizes,
	$qStarts,
	$tStarts,
	$qSeqs,
	$tSeqs) = @blatline_fields;
    
    $tStart++;
    $qStart++;

    my @HSPs;
    $blockSizes =~ s/,\s*$//;
    $qStarts    =~ s/,\s*$//;
    $tStarts    =~ s/,\s*$//;
    $qSeqs    =~ s/,\s*$// if $qSeqs;
    $tSeqs    =~ s/,\s*$// if $tSeqs;    
    my @blockSizeList = split /,\s*/, $blockSizes;
    my @qStartList    = split /,\s*/, $qStarts;
    my @tStartList    = split /,\s*/, $tStarts;
    my @qSeqList      = split /,\s*/, $qSeqs if $qSeqs;
    my @tSeqList      = split /,\s*/, $tSeqs if $tSeqs;
    for my $i (0..$blockCount-1)  {
	push @HSPs, AT::HSP->new
	    (qStart    => $qStartList[$i] + 1,
	     qEnd      => $qStartList[$i] + $blockSizeList[$i],
	     tStart    => $tStartList[$i] + 1,
	     tEnd      => $tStartList[$i] + $blockSizeList[$i],
	     blockSize => $blockSizeList[$i],
	     qSeq      => $qSeqList[$i],
	     tSeq      => $tSeqList[$i]
	     );
    }
    return AT::Mapping->new(bin		 => $bin,
			    matches      => $matches,
			    misMatches   => $misMatches,
			    repMatches   => $repMatches,
			    nCount       => $nCount,
			    qNumInsert   => $qNumInsert,
			    qBaseInsert  => $qBaseInsert,
			    tNumInsert   => $tNumInsert,
			    tBaseInsert  => $tBaseInsert,
			    strand       => $strand,
			    qName        => $qName,
			    qSize        => $qSize,
			    qStart       => $qStart,
			    qEnd         => $qEnd,
			    tName        => $tName,
			    tSize        => $tSize,
			    tStart       => $tStart,
			    tEnd         => $tEnd,
			    blockCount   => $blockCount,
			    HSPs         => [@HSPs]);
}

1;
