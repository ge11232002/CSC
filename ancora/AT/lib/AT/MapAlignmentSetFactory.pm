# AT::MapAlignmentSetFactory module
#
# Copyright Boris Lenhard
# You may distrubute this code under the same terms as perl itself
#
# POD documentation:

=head1 NAME

AT::MapAlignmenSetFactory - factory object that takes a set of MapAlignment
objects and aligns their genomic regions to produce a MapAlignmentSet object

=head2 SYNOPSIS

my $alnsset_factory = AT::MapAlignmentSetFactory->new();
my $aln_set = $alnset_factory->create_set
    (alignments => [$human_aln, $mouse_aln]);

=head2 DESCRIPTION

Given a set of cDNA-genomic alignments for two or more homologous genes, one
may want to find corresponding exons and introns. This module is intended to
be the first step towards that goal, by producing an object
(AT::MapAligmentSet) in which the cDNA-genomic alignments are related
by an alignment of the genomic regions. Genomic regions that overlap
according to their coordinates are not aligned, but simply merged.

Currently, the module does not handle more than two disjunct geomic
regions, i.e. it only does pairwise alignment. However, it should be ok
to pass it many overlapping genomic regions. The only alignment program
currently supported is Blastz (via the AT::Tools::Run::Blastz interface).
Support for bl2seq is halfway implemented, but was interrupted as it
was realized that Blastz performs better.

=head2 TO DO

This module could probably be fused with AT::CMP::CrossMatchFactory, which
is written for a similar purpose.

Currently the genomic alignment stored in AT::MapAlignmentSet is an 
AT::Alignment object. It has not been decided whether we should use
this type of alignment object or a Bio::SimpleAlign object.

The coordinates in the genomic alignment start from 1, i.e. they are not
absolute coordinates with respect to the chromosomes. This should probably be
changed, but perhaps that change should be implemented in the Blastz parser.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# The code begins HERE

package AT::MapAlignmentSetFactory;

use strict;
use vars '@ISA';
use Carp;
use File::Temp qw (tempfile);
use AT::Root;
use AT::MapAlignmentSet;
use AT::Alignment::Composite;
use AT::Alignment::Subalignment;
use AT::Tools::Run::Blastz;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Bio::Tools::Run::StandAloneBlast;
#use Bio::SeqFeature::Generic;
#use Bio::SeqFeature::Gene::Exon;

use constant BLAST_PARAMS => ('program' => 'blastn');

@ISA = qw(AT::Root);

=head2 new

 Title     : new
 Usage     : my $alnset_factory = AT::MapAlignmentSetFactory->new();
 Function  : Constructor
 Returns   : An AT::MapAlignmentSetFactory object
 Args      : Blastz is currently the default and only supported
             alignmentfactory for producing the genomic alignment.
             There are two options intented for specifying a
             different alignmentfactory:
             aln_factory - alignmentfactory object to use.
              (note that an AT::Tools::Run::Blastz object can not
               be passed, as its constructor requires the sequences
               that are to be aligned)
             aln_program - name of the alignment program to
              use, e.g. "blastz".

=cut

sub new  {
    my ($caller, %args) = @_;
    my $self = bless { aln_factory => ($args{aln_factory} or undef) },
               ref $caller || $caller;
}


=head2 create_set

 Title     : create_set
 Usage     : my $alnset = $alnset_factory->create_set();
 Function  : Does all the work
 Returns   : An AT::MapAlignmentSet object
 Args      : Same as for the constructor.

=cut

sub create_set {
    my ($self, %args) = @_;

    # Get alignments from arg
    my $alignments_ref = $args{alignments};
    if(!defined($alignments_ref) || ref($alignments_ref) ne "ARRAY") {
	croak "Missing or invalid 'alignments' argument";
    }
    my @alignments = @$alignments_ref;

    # Get genomic sequence for each alignment.
    # Overlapping or adjacent genomic sequences are replaced with one.
    my @genomic_seqs = $self->_get_genomic_seqs($alignments_ref);

    # Get alignment factory for aligning genomic sequences
    my $aln_factory;
    if($args{aln_factory}) {
	$aln_factory = $args{aln_factory};
    }
    elsif($args{aln_program} || !defined($self->{aln_factory})) {
	$aln_factory = 
	    $self->_create_aln_factory($args{aln_program},@genomic_seqs);
    }
    else {
	$aln_factory = $self->{aln_factory};
    }

    # Align genomic sequences
    my $genomic_aln = $self->_align_genomic_seqs($aln_factory, @genomic_seqs)
	|| return undef;

    # Make set
    my $set = AT::MapAlignmentSet->new( alignments => $alignments_ref,
					genomic_seqs => \@genomic_seqs,
					genomic_aln => $genomic_aln );
    return $set;
}


sub _get_genomic_seqs
{
    my ($self, $alignments) = @_;

    # Get the target seq of each alignment
    my @seqs = map { ($_->each_seq)[0] } @$alignments;

    # Sort seqs by 1.id (database & chromosome), 2.strand, 3.start
    @seqs = sort {$a->id cmp $b->id ||
		      $a->strand <=> $b->strand ||
		      $a->start <=> $b->start} @seqs;

    #print "Input seqs:\n"; _debug_pr_seqs(@seqs);
    
    # The loop below does the following until the list @seqs is empty.
    # 1. Take the first seq off the list
    #    Set the cluster to span that seq.
    # 2. If there is a next seq in the list and it overlaps with the cluster,
    #    extend the cluster with that seq and take it off the list. Repeat
    #    until the next seq in the list does not overlap with the cluster.
    #    (Same strand is a requirement for overlap. Adjacency is counted as
    #    overlap.)
    # 3. Create a sequence that spans the cluster and add it to
    #    @resulting_seqs.

    # May want to add the following functionality: an option which specifies
    # that sequences on the same chromosome should be clustered even if
    # they do not overlap, but are less than x bases apart. To implement
    # that, the method needs the database object for each target sequence,
    # so that the gaps can be filled.

    my @resulting_seqs;
    while (@seqs) {
	my $seq = shift @seqs;
	my $cluster_end = $seq->end;
	my $cluster_seqstr = $self->_ungap($seq->seq);
	my $next_seq;
	while(($next_seq = $seqs[0]) &&
	      $next_seq->id eq $seq->id &&
	      $next_seq->strand == $seq->strand &&
	      $next_seq->start <= $cluster_end+1) {
	    if ($cluster_end < $next_seq->end) {
		$cluster_seqstr .= substr($self->_ungap($next_seq->seq),
					  $cluster_end - $next_seq->start + 1);
		$cluster_end = $next_seq->end;
	    }
	    shift @seqs;
	}
	my $composed_seq = $seq->new(-id => $seq->id,
				     -seq => $cluster_seqstr,
				     -start => $seq->start,
				     -end => $cluster_end,
				     -strand => $seq->strand);
	push @resulting_seqs, $composed_seq;
    }

    #print "Output seqs:\n"; _debug_pr_seqs(@resulting_seqs);

    return @resulting_seqs;
}


sub _ungap
{
    my ($self, $str) = @_;
    $str =~ tr/-//d;
    return $str;
}


sub _debug_pr_seqs
{
    foreach my $s (@_) {
	print 
	    $s->id, "\t",
	    $s->strand, "\t",
	    $s->start, "\t",
	    $s->end, "\n"; 
	#if($s->length != ($s->end - $s->start + 1)) {
	#    warn ("length is ".$s->length."!");
	#}
    }
}


sub _create_aln_factory
{
    my ($self, $program, @seqs) = @_;
    my $factory;

    $program = lc ($program) if defined ($program);

    if(scalar @seqs == 2) {
	if(!defined($program) || $program eq "blastz") {
	    # Default to blastz for 2 sequences
	    $factory = AT::Tools::Run::Blastz->new(-seq1 => $seqs[0],
						   -seq2 => $seqs[1]);
	}
	elsif($program eq "blast") {
	    $factory = Bio::Tools::Run::StandAloneBlast->new(BLAST_PARAMS);
	}
	else {
	    croak "Alignment program $program not supported.";
	}
    }
    elsif(scalar @seqs > 2) {
	croak "More than two target regions not yet supported.";
    }

    return $factory;
}


sub _align_genomic_seqs
{
    my ($self, $factory, @seqs) = @_;
    my $aln;

    if (scalar(@seqs) == 1) {
	$aln = Bio::SimpleAlign->new();
	$aln->add_seq(@seqs);
    }
    elsif($factory->isa("AT::Tools::Run::Blastz")) {
	if (scalar(@seqs) != 2) {
	    croak "Can only use blast alignment factory for two sequences";
	}
	$factory->run;
	$factory->parse;
	my $blastz_alnset = $factory->get_alignmentset;
	if(defined $blastz_alnset->get_alignments) {
	    $blastz_alnset->calculate_score;
	    $blastz_alnset->sort_by_score;
	    $aln = $blastz_alnset->get_best_alignment;
	    if($aln->get_strand eq "-") {
		warn ("Best blastz alignment is +/- for sequences ".
		      $seqs[0]->get_nse." and ".$seqs[1]->get_nse);
		undef $aln;
	    }
	}
    }
    elsif($factory->isa("Bio::Tools::Run::StandAloneBlast")) {
	# Note: this part is unfinished, but should work to some extent.
	# I stopped developing it since blastz obviously performs better.
	if (scalar(@seqs) != 2) {
	    croak "Can only use blast alignment factory for two sequences";
	}
	#print ">seq1\n", $seqs[0]->seq, "\n";
	#print ">seq2\n", $seqs[1]->seq, "\n";
	my ($temp_fh, $temp_fn) = tempfile();
	close($temp_fh);
	print "tempfile: $temp_fn\n";
	$factory->outfile($temp_fn);
	$factory->Strands(1);   # only need to align plus against plus
	# also specify which query strand to use
	my $report = $factory->bl2seq(@seqs);
	$report->queryName($seqs[0]->id);
	my @alns = $self->_bl2seq_report_to_alignments($report);
	$aln = $alns[0];
	unlink($temp_fn);
    }
    else {
	croak (ref($factory)." not supported");
    }

    return $aln;
}       


# There is supposed to be a way of doing this using Bio::AlignIO, but it does
# not work since the Bio::AlignIO::bl2seq->next_aln method causes the program
# to terminate when there are no more alignments in the report.
# Note: this method is under development
sub _bl2seq_report_to_alignments
{
    my ($self, $report, $name1) = @_;
    my @alns;

    while(my $hsp = $report->next_feature) {

	my $seqchar1 = $hsp->querySeq;
	my $start1 = $hsp->query->start;
	my $end1 = $hsp->query->end;
	my $seqchar2 = $hsp->sbjctSeq;
	my $start2 = $hsp->hit->start;
	my $end2 = $hsp->hit->end;
	unless ($seqchar1 && $start1 && $end1 &&
		$seqchar2 && $start2 && $end2) {
	    croak "Error parsing bl2seq report for query ".$report->queryName;
	}

	my $seq1 = new Bio::LocatableSeq('-seq'=>$seqchar1,
				     '-id'=>$report->queryName,
				     '-start'=>$start1,
				     '-end'=>$end1,
				     );
	my $seq2 = new Bio::LocatableSeq('-seq'=>$seqchar2,
				     '-id'=>$report->sbjctName,
				     '-start'=>$start2,
				     '-end'=>$end2,
				     );
	# Note: strands are not set.

	my $aln = Bio::SimpleAlign->new(-source => 'bl2seq');
	$aln->add_seq($seq1);
	$aln->add_seq($seq2);
	push @alns, $aln;

#	print "A ", $seq1->start, "-", $seq1->end,
#	"\t(", $seq1->end-$seq1->start+1,")\t",
#	$seq2->start, "-", $seq2->end, "\t(", $seq2->end-$seq2->start+1,")\t",
#	$hsp->score,"\t",
#	$hsp->bits,"\t",
#	$hsp->percent,"\t",
#	$hsp->P,"\t",
#	$hsp->match,"\t",
#	$hsp->positive,"\t",
#	$hsp->length, "\n";

    }

    return @alns;
}


# For an AT::Alignment object
sub _debug_print_aln
{
    my ($a) = @_;
    print "--- ALIGNMENT ".$a->get_strand." ".$a->get_score." ---\n";
    my $sa = $a->get_subalignments;
    foreach my $aln (@$sa) {
	my ($hs, $he) = @{$aln->get_pos_seq1};
	my ($ms, $me) = @{$aln->get_pos_seq2};
	print "A ", $hs, "-", $he,
	    "\t(", $he-$hs+1,")\t",
	    $ms, "-", $me, "\t(", $me-$ms+1,")\n";
    }
}


# For an array of SimpleAlign objects
sub _debug_print_alns
{
    my @a = @_;
    @a = sort { $a->get_seq_by_pos(1)->start <=> $b->get_seq_by_pos(1)->start}
    @a;
    foreach my $aln (@a) {
	my ($hs, $ms) = $aln->each_seq;
	if ($hs->end - $hs->start > 16) {
	    print "A ", $hs->start, "-", $hs->end,
	    "\t(", $hs->end-$hs->start+1,")\t",
	    $ms->start, "-", $ms->end, "\t(", $ms->end-$ms->start+1,")\n";
	}
    }
}
