# AT::MapAlignmentFactory module
#
# Copyright Boris Lenhard
# You may distrubute this code under the same terms as perl itself
#
# POD documentation:

=head1 NAME

AT::MapAlignmentFactory - factory object that produces 
Bio::SeqAlign objects from mappings and sequence objects

=head2 SYNOPSIS

=head2 DESCRIPTION

The MapAlignmentFactory converts target sequence and mapping information into
AT::MapAlignment objects (subclass of Bio::SimpleAlign objects with added
methods
for intron and exon retrieval. For convenience, it is equipped with methods
to produce alignments from accession numbers and database objects, too.


=cut

# The code begins HERE

package AT::MapAlignmentFactory_old;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use AT::MapAlignment;
use Bio::SimpleAlign;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Gene::Exon;

@ISA = qw(AT::Root);




sub new  {
    my ($caller, %args) = @_;
    my $self = bless {compact => ($args{compact} or 0),
		      target_db   => ($args{target_db} or undef),
		      query_db    => ($args{query_db} or undef),
		      upstr_context => ($args{upstr_context} or 0),
		      downstr_context => ($args{downstr_context} or 0),
		      _resultlist => []},
               ref $caller || $caller;
}


sub run {
    my ($self, %args) = @_;
    my $acc; 
    my $query_db = $args{query_db} 
                    || $self->{query_db}
                    || croak ("No query database set.");
    my $target_db = $args{target_db} 
                    || $self->{target_db}
                    || croak ("No target database set.");
    my @mappings;
    if (my $mapping = $args{mapping}) {
        @mappings = ($mapping);
    }
    elsif (my $mappingsref = $args{mappings})  {
        if (ref($mappingsref) eq "ARRAY")  {
		@mappings = @$mappingsref;
        }
        
        else {
            croak ("mappings argument not an array reference");
        }
    }
    elsif ($acc= $args{acc})  {
        @mappings = $query_db->get_mappings_for_acc($acc);
    }
    else  {
        croak ("no acc, mappings or mapping argument");
    } 
    my $upstr_context = $args{upstr_context} || $self->upstr_context;
    my $downstr_context = $args{downstr_context} || $self->downstr_context;
    $self->_resultlist([]);
    foreach my $mapping (@mappings)  {
        my $qSeqObj = $query_db->get_query_seq($acc or $mapping->qName);
	my ($genome_seq_start, $genome_seq_end);
	if ($mapping->strand eq "+") {
	    $genome_seq_start = $mapping->tStart - $upstr_context;
	    $genome_seq_end = $mapping->tEnd + $downstr_context;
	}
	else {
	    $genome_seq_start = $mapping->tStart - $downstr_context;
	    $genome_seq_end = $mapping->tEnd + $upstr_context;
	}
        my $tSeqObj = $target_db->get_genome_seq(%args, 
						 chr    => $mapping->tName,
                                                 start  => $genome_seq_start,
                                                 end    => $genome_seq_end,
                                                 strand => $mapping->strand);
	push @{$self->_resultlist}, 
	 [ $self->_align_sequences_create_exons_introns
	   (target_seq => $tSeqObj,
	    query_seq  => $qSeqObj,
	    mapping    => $mapping,
	    compact    => ($args{'compact'} or $self->{compact}))
	   ];
    }
    return @{$self->_resultlist};
   
}



sub get_alignments_in_region {
    
}


sub get_alignments {
    my ($self, %args) = @_;
    unless (@{$self->_resultlist}) { $self->run(%args); }
    my @alignments;
    foreach my $resultref ( @{$self->_resultlist} ) {
	my ($aln, $query_exon_ref, $target_exon_ref, $intron_ref) = @$resultref;
	push @alignments, $aln;
    }
    return @alignments;
}



# sub get_query_exons {
#     my ($self, %args) = @_;
#     unless (@{$self->_resultlist}) { $self->run(%args); }
#     my @exons;
#     foreach my $resultref ( @{$self->_resultlist} ) {
# 	my ($aln, $query_exon_ref, $target_exon_ref, $intron_ref) = @$resultref;
# 	push @exons, $query_exon_ref;
#     }
#     return @exons;
# }

# sub get_target_exons {
#     my ($self, %args) = @_;
#     unless (@{$self->_resultlist}) { $self->run(%args); }
#     my @exons;
#     foreach my $resultref ( @{$self->_resultlist} ) {
# 	my ($aln, $query_exon_ref, $target_exon_ref, $intron_ref) = @$resultref;
# 	push @exons, $target_exon_ref;
#     }
#     return @exons;
# }

# sub get_introns {
#     my ($self, %args) = @_;
#     unless (@{$self->_resultlist}) { $self->run(%args); }
#     my @introns;
#     foreach my $resultref ( @{$self->_resultlist} ) {
# 	my ($aln, $query_exon_ref, $target_exon_ref, $intron_ref) = @$resultref;
# 	push @introns, $intron_ref;
#     }
#     return @introns;
# }








#sub construct_alignment {
#    my ($self, %args) = @_;
#    my ($aln, @query_exons, @target_exons, @introns) = 
#	$self->_align_sequences_create_exons_introns(%args); 
#    return $aln;
# }


sub _align_sequences_create_exons_introns  {
    my ($self, %args) = @_;
    my $tSeqobj = $args{'target_seq'} || croak "No genome_seq";
    my $qSeqobj  = $args{'query_seq'} || croak "No query_seq";
    my $mapping      = $args{'mapping'} || croak "No mapping argument";
    my ($min_intron_length) = ($args{'min_intron_length'} or 12);

    my ($prev_qStart, $prev_qEnd, $prev_tEnd);
    my @hsps = ($mapping->all_HSPs);
    # if ($mapping->strand eq "-") { @hsps = reverse @hsps; }
    
    my (@query_exons, @target_exons, @introns, @tSeqbits, @qSeqbits);

    foreach my $hsp (@hsps)  {

        my ($qs, $qe, $ts, $te) = ($hsp->qStart,
                                   $hsp->qEnd,
                                   $hsp->tStart,
                                   $hsp->tEnd);
	if ($tSeqobj->strand == -1) {
		my $rqs = $qSeqobj->length - $qe +1;
		$qe = $qSeqobj->length - $qs +1;
		$qs = $rqs;
	}
        my ($qGap, $tGap) = (0, 0);

        if (defined($prev_qEnd)) {
	    if($tSeqobj->strand == 1 and $qs - $prev_qEnd > 1) {
		# gap in target seq
		push @qSeqbits, $qSeqobj->subseq($prev_qEnd+1, $qs-1);
		$tGap += $qs - $prev_qEnd -1;
		# correct previuos exon length - not pretty
		$query_exons[-1]->end( $query_exons[-1]->end + $tGap );
		$target_exons[-1]->end( $target_exons[-1]->end + $tGap );
	    }
	    elsif($tSeqobj->strand == -1 and $prev_qStart - $qe > 1) {
		# gap in target seq
		push @qSeqbits, $qSeqobj->subseq($qe+1, $prev_qStart-1);
		$tGap += $prev_qStart - $qe -1;
		# correct previuos exon length - not pretty
		$query_exons[-1]->start( $query_exons[-1]->start - $tGap );
		$target_exons[-1]->start( $target_exons[-1]->start - $tGap );
	    }
	}
        if (defined($prev_tEnd) and $ts - $prev_tEnd >1) {
	    # gap in query seq (maybe intron)
            if ($args{compact} and $ts - $prev_tEnd > 15)  {
                push @tSeqbits,
		correct_strand ( $tSeqobj->subseq
				 ($prev_tEnd-$tSeqobj->start+2,
				  $prev_tEnd-$tSeqobj->start+6),
				 $tSeqobj->strand);
		push @tSeqbits, ".....";
		push @tSeqbits,
		correct_strand ( $tSeqobj->subseq($ts-$tSeqobj->start-4,
						  $ts-$tSeqobj->start),
				 $tSeqobj->strand);
		$qGap += 15;                     
            }
            else  {
                push (@tSeqbits,
		      correct_strand( $tSeqobj->subseq
				      ($prev_tEnd-$tSeqobj->start+2,
				       $ts-$tSeqobj->start),
				      $tSeqobj->strand));
                $qGap += $ts - $prev_tEnd -1;
            }
        }  
	# correct gaps
	my $mingap = $qGap<$tGap ? $qGap : $tGap;
	$qGap -= $mingap;
	$tGap -= $mingap;

	# preliminary intron generation 
        if (defined($prev_tEnd) and $qGap>0) {
	    push @introns , Bio::SeqFeature::Generic->new
		(-primary_tag => "intron",
		 -start => $prev_tEnd-$tSeqobj->start+2,
		 #          + abs($qs - $prev_qEnd) -1
		 # ^^^ this row is quick fix for unmatching gaps
		 -end   => $ts-$tSeqobj->start,
		 -strand => $tSeqobj->strand,
		 -tag =>   {abs_start => $prev_tEnd+1,
			    #                     + abs($qs - $prev_qEnd) -1
			    # ^^^ this row is quick fix for unmatching gaps
			    abs_end  => $ts-1 }
		 );
	    $introns[-1]->attach_seq($tSeqobj);
	}
	
        if ($args{compact} and $ts - ($prev_tEnd or 0) > 15 and $qGap == 15) {
            push @qSeqbits, "-----.....-----";
        }
        else {
            push @qSeqbits, "-" x $qGap;
        }
        push @tSeqbits, "-" x $tGap;
         
	push @tSeqbits, correct_strand
	    ($tSeqobj->subseq($ts - $tSeqobj->start + 1,
			      $te - $tSeqobj->start + 1),
	     $tSeqobj->strand );
        push @qSeqbits, $qSeqobj->subseq($qs,$qe);

	# create exons
	push @target_exons , Bio::SeqFeature::Gene::Exon->new
	    (-start => $ts - $tSeqobj->start + 1, 
	     -end => $te - $tSeqobj->start + 1,
	     -strand => $tSeqobj->strand,
	     -tag =>{abs_start => $ts,
		     abs_end   => $te
		     }
	     );
	$target_exons[-1]->attach_seq($tSeqobj);
       
	push @query_exons , Bio::SeqFeature::Gene::Exon->new(-start => $qs, 
							      -end   =>$qe);
	$query_exons[-1]->attach_seq($qSeqobj);
	
	# bookkeeping for the following intron

	$prev_qStart = $qs;
        $prev_qEnd = $qe;
        $prev_tEnd = $te;
    }
    my $context1_length = ($mapping->tStart - $tSeqobj->start);
    my $context1 = ($context1_length) ?
	correct_strand($tSeqobj->subseq(1, $context1_length),
		       $tSeqobj->strand):"";
    my $context2_length = $tSeqobj->end - $mapping->tEnd;
    my $context2 = ($context2_length) ?
	correct_strand($tSeqobj->subseq($tSeqobj->length-$context2_length+1,
					$tSeqobj->length),$tSeqobj->strand):"";
    @tSeqbits = ($context1, @tSeqbits, $context2);
    @qSeqbits = ('-' x $context1_length, @qSeqbits, '-' x $context2_length);
    if ($mapping->strand eq "-") { 
	@tSeqbits = reverse @tSeqbits; 
	@qSeqbits = reverse @qSeqbits; 
    }
    my $tSeqstring = join "", @tSeqbits;
    my $qSeqstring = join "", @qSeqbits;
    
    my $qAlnSeq = Bio::LocatableSeq->new(-id    =>$mapping->qName,
					 -seq   => $qSeqstring,
					 -start => $mapping->qStart,
					 -end   => $mapping->qEnd);
    my $tAlnSeq = Bio::LocatableSeq->new
	(-id    => $tSeqobj->id.($tSeqobj->strand==-1?"(complement)":""),
	 -seq   => $tSeqstring,
	 -start => $tSeqobj->start,
	 -end   => $tSeqobj->end,
	 -strand => $tSeqobj->strand);
    my ($qeL, $teL, $iL) = $self->_coalesce($min_intron_length, \@query_exons,
					    \@target_exons, \@introns);

    if ($mapping->strand eq "-") {
	($qeL, $teL, $iL) = ([reverse @$qeL], [reverse @$teL], [reverse @$iL]);
    }

    my $aln = AT::MapAlignment->new(query_exonL  => $qeL,
				    target_exonL => $teL,
				    intronL      => $iL );
    $aln->add_seq($tAlnSeq);
    $aln->add_seq($qAlnSeq);

    return ($aln);

}



sub correct_strand {
	my ($seq, $strand) = @_;
	if ($strand == -1 or $strand eq "-")  {
	    $seq = join('', reverse split('', $seq));
	    $seq =~ tr/ACGTacgt/TGCAtgca/;
	}
	return $seq;
}


sub _coalesce  {
    # this method eliminates "introns" that are actually small gaps in query sequence
    my ($self, $min_intron_length, $query_exon_ref, $target_exon_ref, $intron_ref) = @_;
    my (@query_exons, @target_exons, @introns);
    my (@query_exon_buffer, @target_exon_buffer);
    my $intron_counter = 0;
#    }
    while (my $curint = $intron_ref->[$intron_counter])  {
	if ($target_exon_ref->[0]->start < $curint->start ) {
	    push @target_exon_buffer, shift @$target_exon_ref;
	    push @query_exon_buffer, shift @$query_exon_ref;
	}
	elsif ($curint->length < $min_intron_length)  {
	    $intron_counter++;
	}
	else {
	    push @introns, $curint;
	    push @target_exons, _merge_exons(@target_exon_buffer) if @target_exon_buffer;
	    push @query_exons,  _merge_exons(@query_exon_buffer)  if @query_exon_buffer ;
	    @target_exon_buffer = ();
	    @query_exon_buffer  = ();
	    $intron_counter++;
	}
    }
    my @remaining_target_exons = (@target_exon_buffer, @$target_exon_ref);
    my @remaining_query_exons  = (@query_exon_buffer, @$query_exon_ref);
	    
    push @target_exons, _merge_exons(@remaining_target_exons) if @remaining_target_exons;
    push @query_exons,  _merge_exons(@remaining_query_exons)  if @remaining_query_exons ;

    return (\@query_exons, \@target_exons, \@introns);

}

sub _merge_exons  {
    my (@exons) = @_;
    return undef unless @exons;
    # print join(" ", ">>>", (map {$_->start."-". $_->end} @exons), "\n") if @exons>1;
    @exons = sort {$a->start <=> $b->start} @exons;
    my $merged_exon = Bio::SeqFeature::Gene::Exon->new(-start=>$exons[0]->start,
						       -end  =>$exons[-1]->end,
						       -strand => $exons[0]->strand);
    if ($exons[0]->has_tag("abs_start") and $exons[-1]->has_tag("abs_end"))  {
	$merged_exon->add_tag_value("abs_start", $exons[0]->each_tag_value("abs_start"));
	$merged_exon->add_tag_value("abs_end", $exons[-1]->each_tag_value("abs_end"));
    }
    $merged_exon->attach_seq($exons[0]->entire_seq);
    return $merged_exon;
						       
}
    

1;
