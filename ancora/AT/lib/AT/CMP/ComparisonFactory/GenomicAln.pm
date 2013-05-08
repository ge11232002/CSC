# AT::CMP::ComparisonFactory::GenomicAln module
#
# Copyright Boris Lenhard
# You may distrubute this code under the same terms as perl itself
#
# POD documentation:

=head1 NAME

AT::CMP::ComparisonFactory::GenomicAln

=head2 SYNOPSIS

 # Create factory and compare two gene predictions
 my $cmp_fac = AT::CMP::ComparisonFactory::GenomicAln->new;
 my $cmp = $cmp_fac->run(objects => [$gene1, $gene2]);
 unless ($cmp) { die "No comparison possible"; }

 # Find the region of gene2 corresponding to the region 86-121 of gene1
 my $cc = $cmp->get_corresp_coords(ref_pos => 1,
                                  rel_start => 86,
                                  rel_end => 121)
 unless ($cc) { print "No corresponding region\n"; }

 # Output absolute and relative coords for corresponding regions
 print join '-', $cc->abs_coords(pos => 1);
 print "\t";
 print join '-', $cc->rel_coords(pos => 1);
 print "\n";
 print join '-', $cc->abs_coords(pos => 2);
 print "\t";
 print join '-', $cc->rel_coords(pos => 2);
 print "\n";


=head2 DESCRIPTION

Creates a comparison of two or more of AT::Prediction::Genomic-compliant objects
by aligning their genomic sequences. Currently, only pairwise local alignments
(using blastz) are supported. More than two objects can be included in a
comparison, but the number of disjunct genomic regions objects are located
on cannot exceed two. If the locations of all objects overlap, the resulting
comparison will have an alignment containing a single sequence. It
should be possible to use a comparison with such an alignment in
the same way as a comparison with a pairwise alignment, but this has
not been thoroughly tested.

=head2 TO DO

=head1 APPENDIX

=cut

# The code begins HERE

package AT::CMP::ComparisonFactory::GenomicAln;

use strict;
use vars '@ISA';
use Carp;
use File::Temp qw (tempfile);
use AT::Root;
use AT::Alignment::Composite;
use AT::Alignment::Subalignment;
use AT::CompositeAlignmentProxy;
use AT::Tools::Run::Blastz;
use AT::Tools::SeqHandler;
use AT::Tools::Colinearizer;
use AT::CMP::Comparison::Generic;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Data::Dumper;
#use Bio::Tools::Run::StandAloneBlast;

use constant BLAST_PARAMS => ('program' => 'blastn');

@ISA = qw(AT::Root);

=head2 new

 Title     : new
 Usage     : my $cmp_factory = AT::CMP::ComparisonFactory::GenomicAln->new();
 Function  : Constructor
 Returns   : An AT::CMP::ComparisonFactory::GenomicAln object
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


=head2 run

 Title     : run
 Usage     : my $comparison = $cmp_factory->run(objects => [ $gene1, $gene2 ]);
 Function  : Does all the work
 Returns   : An AT::CMP::Comparison::Generic-compliant object
 Args      : objects - Ref to array of objects to compare
		       The objects should be AT::Prediction::Genomic-compliant.

=cut


sub run {
    my ($self, %args) = @_;

    # Get objects to compare from arg
    my $objects_ref = $args{objects};
    if(!defined($objects_ref) || ref($objects_ref) ne "ARRAY") {
	croak "Missing or invalid 'objects' argument";
    }
    my @objects = @$objects_ref;

    # Get object handler
    #my $object_handler = _get_object_handler($objects_ref);

    # Get genomic sequence for each object.
    # Overlapping or adjacent genomic sequences are replaced with one.
    my ($genomic_seqs, $object_gseq_map) =
	$self->_nonred_genomic_seqs($objects_ref);
    #my @genomic_seqs = (map {$_->oriented_genomic_seq} (@$objects_ref));
    $self->{_comparing} = $genomic_seqs; # just for debug

    # Get alignment factory for aligning genomic sequences
    my $aln_factory;
    if($args{aln_factory}) {
	$aln_factory = $args{aln_factory};
    }
    elsif($args{aln_program} || !defined($self->{aln_factory})) {
	$aln_factory = 
	    $self->_create_aln_factory($args{aln_program},$args{aln_program_file},
                                       @$genomic_seqs);
    }
    else {
	$aln_factory = $self->{aln_factory};
    }

    # Align genomic sequences
    my $genomic_aln =
	$self->_align_genomic_seqs($aln_factory, @$genomic_seqs)
	or return undef;
   
    if (0) {
	foreach my $obj (@$objects_ref) {
	    print STDERR "$obj\t",$obj->seq->id,"\t",$obj->length,"\n",$obj->seq->seq,"\n";
	}
	foreach my $obj (@$genomic_seqs) {
	    print STDERR "$obj\t",$obj->id,"\t",$obj->length,"\n",$obj->seq,"\n";
	}
	print STDERR Dumper($object_gseq_map);
    }

    # Create comparison object
    my $comparison = AT::CMP::Comparison::Generic->new(
	compared_obj_list => $objects_ref,
	genomic_seq_list => $genomic_seqs,
	genomic_aln => $genomic_aln,
	#genomic_aln_strand => $genomic_aln_strand,
	obj_gseq_map => $object_gseq_map);
    return $comparison;
}


sub _nonred_genomic_seqs
{
    my ($self, $objects_ref) = @_;

    # Sort seqs by 1.id (database & chromosome), 2.strand, 3.start
    my @objects = sort
	{$a->genomic_seq->id cmp $b->genomic_seq->id ||
	 $a->genomic_seq->strand <=> $b->genomic_seq->strand ||
	 $a->genomic_seq->start <=> $b->genomic_seq->start} @$objects_ref;

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
    my %obj2seqidx_map;
    for(my $i = 0; $i < @objects;) {
	my $seq = $objects[$i]->genomic_seq;
	my $cluster_end = $seq->end;
	my $cluster_seqstr = $seq->seq;
	my $j;
	for($j = $i+1; $j < @objects; $j++) {
	    my $next_seq = $objects[$j]->genomic_seq;
	    last if ($next_seq->id ne $seq->id or
		     $next_seq->strand != $seq->strand or
		     $next_seq->start > $cluster_end+1);
	    if ($cluster_end < $next_seq->end) {
		$cluster_seqstr .=
		    $next_seq->subseq($cluster_end - $next_seq->start + 1,
				      $next_seq->length);
		$cluster_end = $next_seq->end;
	    }
	}
	my $composed_seq = $seq->new(-id => $seq->id,
				     -seq => $cluster_seqstr,
				     -start => $seq->start,
				     -end => $cluster_end,
				     -strand => 1);
	    # note: strand defaults to 1; it will be changed later if
	    # the sequence is aligned in the reverse orientation
	push @resulting_seqs, $composed_seq;
	for(my $k = $i; $k < $j; $k++) {
	    $obj2seqidx_map{$objects[$k]} = $#resulting_seqs;
	}
	$i = $j;
    }

    #print "Output seqs:\n"; _debug_pr_seqs(@resulting_seqs);
    return (\@resulting_seqs, \%obj2seqidx_map);
}


#sub _make_object_gseq_map
#{
#    my ($self, $objects, $gseqs) = @_;
#    my %map;

#    foreach my $object (@$objects) {
#	my $obj_gseq = $object->genomic_seq;
#	for(my $i = 0; $i < @$gseqs; $i++) {
#	    if($obj_gseq->id eq $gseq->id and
#	       $obj_gseq->start >= $gseq->start and
#	       $obj_gseq->end <= $gseq->end) {
#		$map{object} = $gseq;
#	    }
#	}
#    }
#    
#}


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
    my ($self, $program, $executable, @seqs) = @_;
    my $factory;

    $program = lc ($program) if defined ($program);

    if(scalar @seqs == 2) {
	if(!defined($program) || $program eq 'blastz') {
	    # Default to blastz for 2 sequences
	    $factory = AT::Tools::Run::Blastz->new(-executable => $executable || undef,
                                                   -seq1 => $seqs[0],
						   -seq2 => $seqs[1],
						   -C => 2, -K => 2500);
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
	    croak "Can only use blastz alignment factory for two sequences";
	}
	$factory->run;
	$factory->parse;
	my $blastz_alnset = $factory->get_alignmentset;
	if(defined $blastz_alnset->get_alignments) {
	    foreach my $raw_aln (@{$blastz_alnset->get_alignments}) {
		my $colin_aln = AT::Tools::Colinearizer
		    ->new(-CA => $raw_aln, -max_overlap => 0)->get_colinear();
		if (!defined($aln) or $colin_aln->get_score > $aln->get_score) {
		    $aln = $colin_aln;
		}
	    }
	    $seqs[1]->strand(-1) if($aln->get_strand eq "-");
	    #$aln = AT::CompositeAlignmentProxy->new(aln => $aln,
	    #					    seqs => \@seqs);
   	    #$self->_debug_print_aln($aln);
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
    my ($self, $a) = @_;
    print STDERR "--- ALIGNMENT ".$a->get_strand." ".$a->get_score." ---\n";
    my $sa = $a->get_subalignments;
    foreach my $aln (@$sa) {
	my $alnobj = $aln->get_alignment_obj();
	my ($hs, $he) = @{$aln->get_pos_seq1};
	my ($ms, $me) = @{$aln->get_pos_seq2};
	print STDERR "A ", $hs, "-", $he,
	    "\t(", $he-$hs+1,")\t",
	    $ms, "-", $me, "\t(", $me-$ms+1,")\t";
	if ($a->get_strand eq '-') {
	    ($ms, $me) = ($self->{_comparing}->[1]->length+1-$me,
			  $self->{_comparing}->[1]->length+1-$ms);
	}
	($hs, $he) = AT::Tools::SeqHandler->rel2abs
		    ($self->{_comparing}->[0], $hs, $he);
	($ms, $me) = AT::Tools::SeqHandler->rel2abs
		    ($self->{_comparing}->[1], $ms, $me);
	print STDERR "A ", $hs, "-", $he,
	    "\t(", $he-$hs+1,")\t",
	    $ms, "-", $me, "\t(", $me-$ms+1,")\n";
	#($hs, $he) = map {$_->start, $_->end} $alnobj->get_seq_by_pos(1);
	#($ms, $me) = map {$_->start, $_->end} $alnobj->get_seq_by_pos(2);
	#print "A ", $hs, "-", $he,
	#    $ms, "-", $me, "\t(", $me-$ms+1,")\n";
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
