#
# Copyright Par Engstrom and Boris Lenhard
#
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::Mapping - mapping of a transcript sequence to genomic sequence

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::Mapping;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use Data::Dumper;


@ISA = qw(AT::Root);

=head2 new

 Title     : new
 Usage     :
 Function  : Constructor
	     Do not call this method directly. Create an object using
	     an AT::BlatpslxIO or AT::DB::GenomeMapping object.
 Returns   : AT::Mapping
 Args      : lots!

=cut

sub new  {
    my ($caller, %args) = @_;
    my $self = bless { db => undef,
		       bin =>  undef,
	               mapping_id => undef,
		       target_db => '',
		       matches => undef,
		       misMatches => undef,
		       repMatches => undef,
		       nCount => undef,
		       qNumInsert => undef,
		       qBaseInsert => undef,
		       tNumInsert => undef,
		       tBaseInsert => undef,
		       strand => undef,
		       qName => undef,
		       qVersion => undef,
		       qSize => undef,
		       qStart => undef,
		       qEnd => undef,
       		       qType => undef,
		       tName => undef,
		       tSize => undef,
		       tStart => undef,
		       tEnd => undef,
		       blockCount => undef,
		       qStarts => undef,
		       tStarts => undef,
		       qSeqs => undef,
		       tSeqs => undef,
		       #library => undef,
		       #mrnaClone => undef,
		       query_seq => undef,
		       #refSeqStatus => undef,
		       HSPs => undef,
		       percentId => undef,
		       score => undef,
		       qReadDir => 0,
		       qInitTs => 0,
		       qTermAs => 0,
		       status => '',
		       status_note => '',
		       mRNAInfo => undef,
		       _HSPs_are_sorted => 0,
		       %args
    }, ref $caller || $caller;
    $self->{_HSPs} = $self->{HSPs} || [];
    delete $self->{HSPs};
    my ($qStrand, $tStrand) = split(//,$self->{strand});
    delete $self->{strand};
    $self->{qStrand} = $qStrand || croak "no strand";
    $self->{tStrand} = $tStrand || ''; # must default to empty string
    return $self;
}


=head2 strand, qStrand, tStrand

 Title     : strand, qStrand, tStrand
 Usage     : $m->strand('-');
 Function  : Get/set strand of mapping as sign
 Returns   : string (e.g. '+','-','++','+-')
 Args      : string

=cut

sub strand {
    my ($self, $value) = @_;
    if($value) {
	my ($qStrand, $tStrand) = split(//, $value);
	$self->qStrand($qStrand);
	$self->tStrand($tStrand);
    }
    return ($self->qStrand).($self->tStrand);
}


=head2 strand_numeric, qStrand_numeric, tStrand_numeric

 Title     : strand_numeric, qStrand_numeric, tStrand_numeric
 Usage     : $m->strand_numeric();
 Function  : Get strand of mapping as number
 Returns   : -1 or 1
 Args      : -

=cut

sub strand_numeric { return $_[0]->qStrand_numeric * $_[0]->tStrand_numeric; }
sub qStrand_numeric { return $_[0]->qStrand eq '-' ? -1 : 1; }
sub tStrand_numeric { return $_[0]->tStrand eq '-' ? -1 : 1; }
# Note: in the above methods, it is important that we default to 1, since
# tStrand may be an empty string and should then be interpreted as '+'.


sub qVersion {
    warn "AT::Mapping->qVersion is deprecated. Use AT::Mapping->mRNAInfo->version\n";
    return $_[0]->mRNAInfo->version;
}

=head2 all_HSPs

 Title     : all_HSPs
 Usage     : my @HSPs = $mapping->all_HSPs;
 Function  : Returns all HSPs of the mapping ordered after tStart *
	     tStrand_numeric (i.e. left->right on the genomic plus strand)
 Returns   : Array of AT::HSP objects
 Args      : -

=cut

sub all_HSPs {
    my ($self) = @_;
    $self->_sort_HSPs_if_required();
    return @{$self->{_HSPs}};
}


=head2 HSP

 Title     : HSP
 Usage     : my $hsp = $mapping->hsp(2);
 Returns   : AT::HSP object
 Args      : Number from 1 to $mapping->nr_hsps

=cut

sub HSP {
    my ($self, $nr) = @_;
    $self->_sort_HSPs_if_required();
    return $self->{_HSPs}->[$nr-1];
}


=head2 nr_HSPs

 Title     : nr_HSPs
 Usage     : print "mapping has ", $mapping->nr_hsps, " HSPs";
 Returns   : Integer
 Args      : -

=cut

sub nr_HSPs { return scalar(@{$_[0]->{_HSPs}}); }


=head2 set_HSPs

 Title     : set_HSPs
 Usage     : $mapping->set_HSPs($hsp_list_ref);
 Function  : Set HSPs
 Returns   : -
 Args      : Reference to array of AT::HSP objects

=cut

sub set_HSPs {
    my ($self, $HSPs);
    $self->{_HSPs} = $HSPs;
    $self->{_HSPs_are_sorted} = 0;
}


sub HSPs {
    croak "Use Mapping methods all_HSPs and set_HSPs to access HSPs; do not "
	."attempt to call a method named HSPs";
}


sub _sort_HSPs_if_required
{
    my ($self) = @_;
    unless($self->{_HSPs_are_sorted}) {
	@{$self->{_HSPs}} = sort {$a->tStart <=> $b->tStart} @{$self->{_HSPs}};
	@{$self->{_HSPs}} = reverse(@{$self->{_HSPs}}) if ($self->tStrand eq '-');
	$self->{_HSPs_are_sorted} = 1;
    }
}


=head2 loc_str

 Title     : loc_str
 Usage     : print  "Mapping at ", $m->loc_str, "\n";
 Function  : Return the the chromosomal location of the mapping
	     as a string.
 Returns   : A string on the form 'chr14:1112000-1114000:+'
 Args      : -

=cut

sub loc_str
{
    my ($self) = @_;

    my $str = $self->tName.":".$self->tStart."-".$self->tEnd.":".
	$self->strand;
    return $str;
}


=head2 abs_HSP_tPos

 Title     : abs_HSP_tPos
 Usage     : my @hsp_coords = $mapping->abs_HSP_tPos;
 Function  : Return absolute (plus-strand) genomic coords of HSPs in
	     increasing order.
	     This is useful since, for mappings from translated searches,
	     the $hsp->tStart and $hsp->tEnd methods may give coords that
	     refer to the minus strand.
 Returns   : 2-dimensional array.
 Args      :
 
=cut

sub abs_HSP_tPos {
    my ($self) = @_;
    if($self->tStrand eq '-') {
	return (map {[$self->tSize - $_->tEnd + 1, $self->tSize - $_->tStart + 1]}
		$self->all_HSPs);
    }
    else {
	return (map {[$_->tStart, $_->tEnd]} $self->all_HSPs);
    }
}


=head2 tInsert_coords

 Title     : tInsert_coords
 Usage     : my @insert_coords = $mapping->tInsert_coords
		    (min_tInsert_length => 6);
 Function  : Return genomic coordinates of target inserts
	     (putative introns).
	     An array of coordinate pairs is returned. Each pair is
	     an array of two integers, the first being start of the
	     insert and the second the end.
	     The order of coordinate pairs follows the direction of
	     transcription (as indicated by the strand of the mapping).
	     However, for each coord pair start is always < end.
 Returns   : 2-dimensional array.
 Args      : min_tInsert_length    Min length an insert must have to
		be included. Optional. Defaults to 1.

=cut

sub tInsert_coords
{
    my ($self, %args) = @_;
    my $min_len = ($args{min_tInsert_length} || 1);
    my @coords;
    my @hsp_coords = $self->abs_HSP_tPos;
    for(my $i = 1; $i < @hsp_coords; $i++) {
	my $start = $hsp_coords[$i-1]->[1] + 1;
	my $end = $hsp_coords[$i]->[0] - 1;
	my $len = $end - $start + 1;
	next if ($len < $min_len);
	push @coords, [ $start, $end ];
    }
    @coords = reverse @coords if ($self->strand_numeric == -1);
    return @coords;
}


=head2 attach_query_seq

 Title     : attach_query_seq
 Usage     : $mapping->attach_query_seq($seq);
 Function  : Attach a Bio::Seq object corresponding to the
	     cDNA/EST sequence mapped.
 Returns   : 1
 Args      : Bio::Seq

=cut

sub attach_query_seq {
    my ($self, $query_seq) = @_;
    $self->{query_seq} = $query_seq;
    return 1;
}


=head2 query_seq

 Title     : query_seq
 Usage     : my $seq = $mapping->query_seq;
 Function  : Return attached query sequence if any
 Returns   : Bio::Seq
 Args      : -

=cut

sub query_seq {
    my ($self) = @_;
    unless ($self->{query_seq}) {
	#warn "automatic attachment of query seq to mapping currently unsupported";
	#return undef;
	return undef unless ($self->{db});
	my $version = $self->mRNAInfo ? $self->mRNAInfo->version : undef;
	my $query_seq = $self->db->get_query_seq($self->qName, $version);
	$self->attach_query_seq($query_seq);
    }
    return $self->{query_seq};
}


=head2 set_hsp_flanks

 Title     : set_hsp_flanks
 Usage     : $mapping->set_hsp_flanks($assembly_db);
 Function  : Retrieves flanking dinucleotide sequences (putative
	     splice junctions) for each HSP from the genome assembly
	     db given as argument.
	     To get the sequences, use the lFLank and rFlank
	     methods in the HSP class.
	     The flanking sequence retrieval should be made
	     more automatic, e.g. rFlank and lFlank could be
	     made fields in the HSP table of mapping databases.
 Returns   : -
 Args      : AT::DB::GenomeAsssembly, AT::DB::GenomeAssemblyNibs
	     or compatible.

=cut

sub set_hsp_flanks {
    my ($self, $assembly_db) = @_;
    my @hsps = $self->all_HSPs;
    my @hsp_pos = $self->abs_HSP_tPos;
    foreach my $i (0..@hsps-1) {
	my $left = $assembly_db->get_genome_seq_str(chr => $self->tName,
						    start => $hsp_pos[$i]->[0] - 2,
						    end => $hsp_pos[$i]->[0] - 1);
	my $right = $assembly_db->get_genome_seq_str(chr => $self->tName,
		 				     start => $hsp_pos[$i]->[1] + 1,
						     end => $hsp_pos[$i]->[1] + 2);
	$hsps[$i]->lFlank($left);
	$hsps[$i]->rFlank($right);
    }
}


sub set_hsp_flanks_from_seq {
    my ($self, $seq) = @_;
    my @hsps = $self->all_HSPs;
    my @hsp_pos = $self->abs_HSP_tPos;
    my $offset = $seq->start - 1;
    $hsps[0]->lFlank('nn');	# always setting this to nn avoids problems with mappings having tStart=1
    foreach my $i (0..@hsps-2) {
	my $gap_start = $seq->subseq($hsp_pos[$i]->[1] + 1 - $offset,
				 $hsp_pos[$i]->[1] + 2 - $offset);
	my $gap_end = $seq->subseq($hsp_pos[$i+1]->[0] - 2 - $offset,
				$hsp_pos[$i+1]->[0] - 1 - $offset);
	$hsps[$i]->rFlank($gap_start);
	$hsps[$i+1]->lFlank($gap_end);
    }
    $hsps[-1]->rFlank('nn');	# always setting this to nn avoids problems with mappings having tEnd=tSize
}


=head2 mRNAInfo

 Title     : mRNAInfo
 Usage     : $mapping->mRNAInfo($info);
 Function  : Get/set mRNAInfo.
	     mRNAInfo is an instance of AT::mRNAInfo containing
	     various data for the mRNA sequence mapped, e.g. version,
	     clone, library, cds.
 Returns   : AT::mRNAInfo
 Args      : AT::mRNAInfo (optional)

=cut

sub mRNAInfo {
    my ($self, $value) = @_;
    if($value) { $self->{mRNAInfo} = $value; }
    elsif(!defined($self->{mRNAInfo})) {
	$self->{mRNAInfo} = $self->db->get_mRNA_info($self->qName);
    }
    return $self->{mRNAInfo} || undef;
}


sub _get_library_and_clone_ids {
    my ($self, $value) = @_;
    my ($lib, $clone) = $self->db->get_library_and_clone_ids($self->qName);
    $self->{library} = ($lib || 0) unless (defined $self->{library});
    $self->{mrnaClone} = ($clone || 0) unless (defined $self->{mrnaClone});
}


# 'exact' means that this method will return undef if the given tPos is in a gap; 
# this is all we need now anyway
# should implement this as a binary search
sub qPos_from_tPos_exact {
    my ($self, $tPos, %args) = @_;
    $tPos = $self->tSize - $tPos + 1 if($self->tStrand eq '-');
    foreach my $hsp ($self->all_HSPs) {
	return undef if($hsp->tStart > $tPos); ## <- ok?
	if($hsp->tEnd >= $tPos) {
	    return $hsp->qPos_from_tPos_exact($tPos);
        }
    }
    return undef;
}


sub tRange_from_qRange {
    my ($m, $qRangeStart, $qRangeEnd) = @_;
    my ($tRangeStart, $tRangeEnd);
    my @hsps = $m->all_HSPs;
    if($m->strand eq '-') {
	my $qSize = $m->qSize;
	($qRangeStart, $qRangeEnd) = ($qSize - $qRangeEnd + 1, $qSize - $qRangeStart + 1);
    }
    foreach my $hsp (@hsps) {
        if($hsp->qStart >= $qRangeStart) {
	    $tRangeStart = $hsp->tStart;
	    last;
	}
	elsif($hsp->qEnd >= $qRangeStart) {
	    $tRangeStart = $hsp->tPos_from_qPos_exact($qRangeStart);
	    last;
	}
    }
    return undef unless $tRangeStart;
    foreach my $hsp (reverse @hsps) {
	if($hsp->qEnd <= $qRangeEnd) {
	    $tRangeEnd = $hsp->tEnd;
	    last;
	}
	elsif($hsp->qStart <= $qRangeEnd) {
	    $tRangeEnd = $hsp->tPos_from_qPos_exact($qRangeEnd);
	    last;
	}
    }
    return undef unless $tRangeEnd;
    return $tRangeStart <= $tRangeEnd ? ($tRangeStart, $tRangeEnd) : undef;
}


sub qRange_from_tRange {
    my ($m, $tRangeStart, $tRangeEnd) = @_;
    my ($qRangeStart, $qRangeEnd);
    my @hsps = $m->all_HSPs;
    foreach my $hsp (@hsps) {
        if($hsp->tStart >= $tRangeStart) {
	    $qRangeStart = $hsp->qStart;
	    last;
	}
	elsif($hsp->tEnd >= $tRangeStart) {
	    $qRangeStart = $hsp->qPos_from_tPos_exact($tRangeStart);
	    last;
	}
    }
    return undef unless $qRangeStart;
    foreach my $hsp (reverse @hsps) {
	if($hsp->tEnd <= $tRangeEnd) {
	    $qRangeEnd = $hsp->qEnd;
	    last;
	}
	elsif($hsp->tStart <= $tRangeEnd) {
	    $qRangeEnd = $hsp->qPos_from_tPos_exact($tRangeEnd);
	    last;
	}
    }
    if($m->strand eq '-') {
	my $qSize = $m->qSize;
	($qRangeStart, $qRangeEnd) = ($qSize - $qRangeEnd + 1, $qSize - $qRangeStart + 1);
    }
    return undef unless $qRangeEnd;
    return $qRangeStart <= $qRangeEnd ? ($qRangeStart, $qRangeEnd) : undef;
}



# This method gives percent-id as proposed by Kent.
# This way of calculation may not be optimal when working with cDNA
# mappings, as target inserts are scored negatively.

sub percent_id {
    my ($self) = @_;
    my $a = $self->matches + $self->repMatches;
    my $b = $a + $self->misMatches + $self->qNumInsert+ $self->tNumInsert;
    return 100 * $a / $b;
}


=head2 score_A1

 Title     : score_A1
 Usage     : my $score = $mapping->score_A1;
 Function  : Calculate a score:
		matches + repMatches - misMatches
		- a gap penalty
		+ an intron reward
	     This is an experimental score that seems to work
	     slightly better than UCSCs various scores (below).
	     The gap penalty is currently log2(gap length)+1.
	     If a gap looks like an intron (i.e. gap is only in query,
	     at least 13 bp, and has GT-AG, GC-AG or AT-AC splice-
	     junctions), a reward of 5 points is given instead of a
	     penalty.
	     To make it easier to handle, the resulting score is
	     multiplied by 10 and rounded to the nearest integer.
	     NOTE: this method requires that set_hsp_flanks has
	           been called.
 Returns   : Integer

=cut

sub score_A1 {
    my ($self, %args) = @_;
    my $seq = $args{'seq'};
    my ($intron_reward, $insert_penalty) = (0,0);
    my @hsps = $self->all_HSPs;
    @hsps = reverse(@hsps) if($self->tStrand eq '-');
    for my $i (0 .. @hsps - 2) {
	my $qInsert_len = $hsps[$i+1]->qStart - $hsps[$i]->qEnd - 1;
	my $tInsert_len = $hsps[$i+1]->tStart - $hsps[$i]->tEnd - 1;
	if($qInsert_len == 0 and $tInsert_len >= 13 and
	   $self->_canonical_splice_junctions($hsps[$i], $hsps[$i+1])) {
	    $intron_reward += 5;
	}
	else {
	    my $max_insert_len = $qInsert_len > $tInsert_len ? $qInsert_len : $tInsert_len;
	    $insert_penalty += log($max_insert_len)/log(2) + 1
		if($max_insert_len); # need this if() since blat sometimes gives empty gaps
	}
    }
    #print STDERR "[$intron_reward $qInsert_penalty $tInsert_penalty]\n";
    my $score = $self->matches + $self->repMatches - $self->misMatches
	    + $intron_reward - $insert_penalty;
    return int(10*$score+.5);
}


sub _canonical_splice_junctions {
    my ($self, $hsp1, $hsp2) = @_;
    my $pair = lc($hsp1->rFlank.$hsp2->lFlank);
    if($pair eq 'gtag' or	# Major
       $pair eq 'ctac' or
       $pair eq 'gcag' or	# Minor
       $pair eq 'ctgc' or
       $pair eq 'atac' or	# Very minor
       $pair eq 'gtat')
    {
	return 1;
    }
    return 0;
}


=head2 psl_score

 Title     : psl_score
 Usage     : my $psl_score = $mapping->psl_score;
 Function  : Calculate a simple score:
		matches + repMatches
		- misMatches - qNumInsert - tNuminsert
	     This should equal the score given by UCSC's
	     online blat service.
	     (Ported from Kent's psl.c)
 Returns   : Integer

=cut

sub psl_score {
    my ($self) = @_;
    return $self->matches + $self->repMatches - $self->misMatches
	- $self->qNumInsert - $self->tNumInsert;
}


=head2 psl_milli_score

 Title     : psl_milli_score
 Usage     : my $psl_score = $mapping->psl_milli_score;
 Function  : Calculate a somewhat more sophisticated score.
	     This should equal the 'alignment ratio' calculated
	     by Kent_s pslReps.
 Returns   : Integer between 0 and 1000.
 Args	   : If the query is a genomic sequence rather than an
	     mRNA give 'unspliced' => 1.

=cut

sub psl_milli_score { 1000 - psl_milli_bad(@_); }


=head2 psl_sized_score

 Title     : psl_sized_score
 Usage     : my $psl_score = $mapping->psl_sized_score;
 Function  : An extension of the psl_milli_score that
	     takes alignment length and intron content
	     into account.
	     This should equal the score used to rank mappings
	     in Kent_s pslReps.
 Returns   : If the query is a genomic sequence rather than an
	     mRNA give 'unspliced' => 1.

=cut

sub psl_sized_score { psl_milli_score(@_) + _psl_size_factor(@_); }


=head2 psl_milli_bad

 Title     : psl_milli_bad
 Usage     : my $psl_milli_bad = $mapping->psl_milli_bad;
 Function  : Equal to 1000 - psl_milli_score
 Returns   : Integer between 0 and 1000.

=cut

sub psl_milli_bad {
    my ($self, %args) = @_;
    my $is_mrna = !($args{'unspliced'} || 0);
    my $q_ali_size = $self->qEnd - $self->qStart + 1;
    my $t_ali_size = $self->tEnd - $self->tStart + 1;
    my $ali_size = ($q_ali_size > $t_ali_size) ? $q_ali_size : $t_ali_size;
    die "Invalid qStart/qEnd or tStart/tEnd for mapping_id ".$self->mapping_id
	if($ali_size <= 0);
    my $size_dif = $q_ali_size - $t_ali_size;
    if($size_dif < 0) {
	$size_dif = $is_mrna ? 0 : -$size_dif;
    }
    my $insert_factor = $self->qNumInsert;
    $insert_factor += $self->tNumInsert if ($is_mrna);
    my $milli_bad =
	(1000 *
	    ($self->misMatches + $insert_factor + int(3*log(1+$size_dif)+.5)))
	/
	($self->matches + $self->repMatches + $self->misMatches);
    return int($milli_bad + .5);
}


# Figure out bonus for having introns (from pslReps.c)
sub _psl_intron_factor
{
    my ($self) = @_;
    my $bonus = 0;
    my @hsps = $self->all_HSPs;
    my @hsp_pos = $self->abs_HSP_tPos;
    for my $i (0 .. @hsps - 2) {
	$bonus += 3 if($hsps[$i+1]->qStart - $hsps[$i]->qEnd == 1 and
		       $hsp_pos[$i+1]->[0] - $hsp_pos[$i]->[1] >= 30);
    }
    return $bonus;
}


# Return a factor that will favor longer alignments (from pslReps.c)
sub _psl_size_factor
{
    my ($self, %args) = @_;
    my $score = 4 * int(sqrt($self->matches + $self->repMatches / 4)+.5);
    $score += $self->_psl_intron_factor unless ($args{'unspliced'});
    return $score;
}

=head2 transcript_feature_object

 Title     : transcript_feature_object
 Usage     : my $transcript = $mapping->transcript_feature_object
 Arguments : none
 Function  : converts the mapping to bioperl
               with HSPs as Bio::SeqFeature::Gene::Exon subfeatures
 Returns   : a Bio::SeqFeature::Gene::Transcript object

=cut

sub transcript_feature_object {
    my ($self, $merge_limit) = (@_, 4);
    #my @exons;
    my @hsps = sort {$a->tStart<=>$b->tStart} ($self->all_HSPs);
    my $transcript = Bio::SeqFeature::Gene::Transcript->new
                        (-start => $hsps[0]->tStart, -end => $hsps[-1]->tEnd,
                         -source_tag=>$self->qName,
			 -strand => ($self->strand eq "-" ? -1 : 1),
			 -tag  => {template_tag => $self->tName,
				    query_tag    => $self->qName });
    while (my $hsp = shift @hsps) {
        my $start = $hsp->tStart+1;
        while ($hsps[0] and $hsps[0]->tStart- $hsp->tEnd <$merge_limit)  {
            $hsp = shift @hsps;
        }
        my $end = $hsp->tEnd;


        $transcript->add_exon(
        Bio::SeqFeature::Gene::Exon->new
            (-start => $start,
             -end   => $end,
             -source => $self->qName,
             -strand => $self->strand eq "+" ? 1 : -1));
                #it would be much cleaner if we could do something like
                #   $transcript->add_exon($hsp->exon_feature_object)
                # but we need to merge the hsps with short gaps between them

    }
     #return \@exons;
     return $transcript;
}




=head2 more methods!

In addition to the methods mentioned above,
the class also has the following get/set methods corrseponding to
columns in the PSL format and MAPPING SQL table
    mapping_id
    bin
    target_db
    matches
    misMatches
    repMatches
    nCount
    qNumInsert
    qBaseInsert
    tNumInsert
    tBaseInsert
    qName
    qSize
    qStart
    qEnd
    qType
    tName
    tSize
    tStart
    tEnd
    blockCount
    
=cut

1;

