# AT::MapAlignmentFactory module
#
# Copyright Pär Engström and Boris Lenhard
#
# You may distrubute this code under the same terms as perl itself
#
# POD documentation:

=head1 NAME

AT::MapAlignmentFactory - factory object that produces Bio::SimpleAlign-compliant objects from mappings and sequence objects

=head2 SYNOPSIS

 # Connect to databases

 my $target_db = AT::DB::GenomeAssembly->connect
    ( -dbname => "HS_JUN02",
      -dbhost => "nautilus.cgb.ki.se",
      -dbuser => "at_read",
      -dbpass => "***");

 my $query_db = AT::DB::GenomeMapping->connect
    ( -dbname => "AT_HS_JUN02",
      -dbhost => "nautilus.cgb.ki.se",
      -dbuser => "at_read",
      -dbpass => "***");

 # Create an aligment factory

 my $af = AT::MapAlignmentFactory_mul->new
    ( compact => 1,
      target_db => $target_db,
      query_db => $query_db,
      upstr_context => 50,
      downstr_context => 50 );

 # Create and retrieve alignments for all mappings of cDNA AB046859

 $af->run( acc => 'AB046859' );

 my @alns = $af->get_alignments();

=head2 DESCRIPTION

The MapAlignmentFactory converts target sequence and mapping information into
AT::MapAlignment objects (AT::MapAlignment is a subclass of Bio::SimpleAlign
with added methods for intron and exon retrieval). For convenience, it is
equipped with methods to produce alignments from accession numbers and database objects, too.

The feature annotation is directly derived from the mappings. HSPs are
annotated as exons. Query gaps (where two HSPs are adjacent in the genome,
but not in the cDNA) are annotated as introns. In addition, target
gaps (where two HSPs are adjacent in the cDNA, but not in the genome) and
mismatches (where two HSPs are not adjacent in the genome nor the cDNA)
are annotated.
If query gaps are smaller than a minimum intron length (12 bp by default), they
are fused with the surrounding HSPs into one exon. Similarly, there is
a minimum length for target gaps and mismatches (also 12 bp by default).
If the genomic length of a mismatch is longer or equal to this minimum length,
but the cDNA length of the mismatch is not, the mismatch is annotated as an intron in the genomic sequence and as a (short) mismatch in the cDNA sequence.

For several mappings that all map to the same target sequence, a multiple
sequence alignment can be made. Note that this is not a true multiple
sequence alignment, simply a combination of pairwise cDNA-genome alignments.
When a msa is made, exons, introns and other features are not annotated.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# The code begins HERE

package AT::MapAlignmentFactory;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use AT::MapAlignment;
use Bio::SimpleAlign;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Gene::Exon;

use constant DEFAULT_MIN_INTRON_LENGTH => scalar 12;
use constant DEFAULT_MIN_MISMATCH_LENGTH => scalar 12;

use constant Q_UNDEF    => scalar 0;
use constant Q_GAP      => scalar 1;
use constant Q_MATCH    => scalar 2;
use constant Q_MISMATCH => scalar 4;
use constant Q_TGAP     => scalar 8;

@ISA = qw(AT::Root);



=head2 new

 Title     : new
 Usage     : my $af = AT::MapAlignmentFactory->new(target_db => $db_object);
 Function  : Constructor
 Returns   : An AT::MapAlignmentFactory object
 Args      : target_db - target database to use
               (an AT::DB::GenomeAssembly object)
             query_db - query database to use
               (an AT::DB::GenomeMapping object)
	     query_seq_db - query sequence database to use
	       (an AT::DB::TranscriptSeq object)
             compact - set to 1 for compact mode; default 0
               In compact mode, query gaps are shortened.
             upstr_context - number of bases of upstream genomic
               context to include in alignment; default 0
             downstr_context - number of bases of downstream genomic
               context to include in alignment; default 0
	     revcom_string - string to append to identifiers of
	       revcom'd sequences; default: '(complement)'

=cut

sub new  {
    my ($caller, %args) = @_;

    my $self = bless {compact     => ($args{compact} or 0),
		      target_db   => ($args{target_db} or undef),
		      query_db    => ($args{query_db} or undef),
		      query_seq_db => ($args{query_seq_db} or undef),
		      upstr_context   => ($args{upstr_context} or 0),
		      downstr_context => ($args{downstr_context} or 0),
		      revcom_string => (defined $args{revcom_string} ? $args{revcom_string} : '(complement)'),
		      _resultlist => []},
               ref $caller || $caller;
}


=head2 run

 Title     : run
 Usage     : $af->run(mapping => $mapping);
 Function  : Creates alignments
 Returns   : A list of results that propably should not be used.
             Instead, use the method get_alignments to retrieve the
             alignments.
 Args      : Accepts all the arguments that the constructor accepts.
             Arguments given to run overrides arguments given
             to the constructor.

             In addition, one (and only one) of the following three
             arguments must be supplied:
             acc - query sequence accession number
               All mappings for the accession nr will be retrieved
               and used to construct alignments. A query database
               object must be supplied for this to work.
             acca - reference to array of query sequence accession
	       numbers. Alignments are made for all mappings of the
	       accessions numbers.
             mapping - mapping object to make an alignment for
             mappings - reference to an array of mapping objects
               to make aligments for

             Additional optional arguments:
             combine - for multiple seq alignment, set to
	        'all' to combine all mappings into one alignment
		'overlapping' to combine overlapping mappings
		Leave false/undef to make one alignment per mapping
             min_intron_length - see class description above
             min_mismatch_length - see class description above

=cut

sub run {
    my ($self, %args) = @_;

    # Get some values from args
    my $upstr_context = $args{upstr_context} || $self->upstr_context;
    my $downstr_context = $args{downstr_context} || $self->downstr_context;
    my $qseq_db = $args{query_seq_db} || $self->query_seq_db;

    # Get sets of mappings from args
    # One alignment will be made for each set of mappings
    my $mapping_sets = $self->_get_mapping_sets(%args);

    # Empty the list of results
    $self->_resultlist([]);

    # Generate an alignment for each mapping set
    foreach my $mapping_set (@$mapping_sets) {
	my @mappings = @$mapping_set;

	# Make sure query sequence is attached to each mapping
	foreach my $m (@mappings) {
	    unless($m->query_seq) {
		unless($qseq_db) { croak "Need query_seq_db"; }
		my $qseq = $qseq_db->get_seq($m->qName,
				  	     $m->mRNAInfo ? $m->mRNAInfo->version : undef);
		unless($qseq) {
		    warn("No sequence for ".($m->qName)." in db ".($qseq_db->dbname));
		    # should remove mapping from set - program will fail later
		}
		else {
		    $m->attach_query_seq($qseq);
		}
	    }
	}

	# Go through mappings to get genomic boundaries
	my $tStart = $mappings[0]->tStart;
	my $tEnd = $mappings[0]->tEnd;
	for(my $i = 1; $i < @mappings; $i++) {
	    $tStart = $mappings[$i]->tStart if ($mappings[$i]->tStart < $tStart);
    	    $tEnd = $mappings[$i]->tEnd if ($mappings[$i]->tEnd > $tEnd);
	}

	# Expand genomic boundaries using context values
	if ($mappings[0]->strand eq "+") {
	    $tStart -= $upstr_context;
	    $tEnd += $downstr_context;
	}
	else {
	    $tStart -= $downstr_context;
	    $tEnd += $upstr_context;
	}

	if($tEnd - $tStart > 10e6) {
	    croak "Genomic region is too large (>10 Mb)";
	}

	# Get genomic sequence
	my $tSeqObj;
	if($tSeqObj = $args{'target_seq'}) {
	    unless ($tSeqObj->start <= $tStart and $tSeqObj->end >= $tEnd) {
		# Should also check target_db name here
		croak "Target sequence does not cover mappings";
	    }
	}
	else {
	    my $target_db = $args{target_db} || $self->{target_db} ||
		croak("No target database set");
	    $tSeqObj = $target_db->get_genome_seq(%args,
						  chr => $mappings[0]->tName,
						  start => $tStart,
						  end => $tEnd,
						  strand =>
						  $mappings[0]->strand);
	}

	# Run!
	push @{$self->_resultlist},
	[ $self->_align_sequences_create_exons_introns
	  (target_seq => $tSeqObj,
	   mappings   => \@mappings,
	   compact    => ($args{'compact'} or $self->{compact}))
	  ];
    }

    # Return list of results
    return @{$self->_resultlist};
}


=head2 get_alignments

 Title     : get_alignments
 Usage     : my @alns = $af->get_alignments();
 Function  : Retrieves alignments made by a previous call to the
             run method. If no alignments have been made, calls
             run using the supplied arguments.
 Returns   : An array of AT::MapAlignment objects
 Args      : See the run method.

=cut


sub get_alignments {
    my ($self, %args) = @_;
    unless (@{$self->_resultlist}) { $self->run(%args); }
    my @alignments;
    foreach my $resultref ( @{$self->_resultlist} ) {
	my ($aln, $query_exon_ref, $target_exon_ref, $intron_ref) =
	    @$resultref;
	push @alignments, $aln;
    }
    return @alignments;
}


#
# PRIVATE METHODS FOLLOW
#


sub _get_mapping_sets
# Get sets of mappings (one alignment is made for each set)
{
    my ($self, %args) = @_;

    # Check for misuse of args
    unless($args{acc} xor $args{accs} xor $args{'mappings'} xor $args{'mapping'}) {
	croak("Need one (and only one) of arguments 'acc', 'accs', 'mappings' or 'mapping'");
    }

    # Get mappings to use
    my $mappings;  # $mappings will be reference to array of mappings
    if($args{'acc'}) {
	my $acc = $args{'acc'};
	my $query_db = $args{query_db} || $self->{query_db} ||
	    croak ("No query database set");
	$mappings = [ $query_db->get_mappings_for_acc($acc) ];

    }
    elsif($args{'accs'}) {
	my $accsref = $args{'accs'};
        unless (ref($accsref) eq "ARRAY")  {
            croak ("accs argument not an array reference");
        }
	my $query_db = $args{query_db} || $self->{query_db} ||
	    croak ("No query database set");
	$mappings = [ map { $query_db->get_mappings_for_acc($_) } @$accsref ];
    }
    elsif($args{'mapping'}) {
        $mappings = [ $args{'mapping'} ];
    }
    elsif($args{'mappings'})  {
	my $mappingsref = $args{'mappings'};
        unless (ref($mappingsref) eq "ARRAY")  {
            croak ("Mappings argument not an array reference");
        }
	$mappings = $mappingsref;
    }
    else {
	croak "Argument parsing error";
    }

    # Put mappings into sets
    my $mapping_sets;
    my $combine = $args{'combine'};
    if($combine) {
	if($combine eq 'overlapping') {
	    # Create one set for each group of overlapping mappings
            $mapping_sets = $self->_get_sets_of_overlapping_mappings($mappings);   
	}
	else {
	    # Put all mappings in one set
	    warn "Invalid value of 'combined' argument; assuming value 'all'"
		unless($combine eq 'all');
	    croak "Cannot combine mappings to different chromosomes into one alignment"
		unless($self->_mappings_have_same_target($mappings));
	    $mapping_sets = [ $mappings ];   
	}
    }
    else {
	# Create one set for each mapping
        $mapping_sets = [ map { [$_] } @$mappings ];
    }

    # Return mapping sets
    return $mapping_sets;
}


sub _get_sets_of_overlapping_mappings
{
    my ($self, $mappings_ref) = @_;

    my @mapping_sets;
    
    # Sort mappings by 1.assembly, 2.chromosome, 3.start
    my @mappings = sort {$a->target_db cmp $b->target_db ||
		         $a->tName cmp $b->tName ||
		         $a->tStart <=> $b->tStart} @$mappings_ref;

    # Initialize first set with first mapping
    my $set = [shift @mappings];
    my $set_end = $set->[0]->tEnd;
    push @mapping_sets, $set;

    # Make sets
    while(@mappings) {
	my $m = shift @mappings;
	# Add mapping to set if it overlaps with set
	if($m->target_db eq $set->[0]->target_db and
	   $m->tName eq $set->[0]->tName and
	   $m->tStart <= $set_end) {
	    push @$set, $m;
	    $set_end = $m->tEnd if($set_end < $m->tEnd);
	}
	# Else make new set from mapping
	else {
	    $set = [$m];
	    my $set_end = $m->tEnd;
	    push @mapping_sets, $set;
	}
    }

    # Return sets
    return \@mapping_sets;
}


sub _mappings_have_same_target
# Go through mappings to check that they are all to the same chromosome
{
    my ($self, $mappings) = @_;
    my $target_db = $mappings->[0]->target_db;
    my $chr = $mappings->[0]->tName;
    for(my $i = 1; $i < @$mappings; $i++) {
	if($target_db ne $mappings->[$i]->target_db or $chr ne $mappings->[$i]->tName) {
	    return 0;
	}
    }
    return 1;    
}


sub _align_sequences_create_exons_introns  {
    my ($self, %args) = @_;
    my $tSeqObj = $args{'target_seq'} || croak "No genome_seq";
    my $mappings = $args{'mappings'} || croak "No mappings";

    my $strand = $tSeqObj->strand;

    # Create description of alignment for internal use
    my $chunks = $self->_get_aligned_chunks($mappings, $strand);

    # Create gapped sequences
    my @seqs = $self->_create_gapped_sequences(chunks => $chunks, %args);

    # Create exons, introns etc
    my %feats;
    unless(scalar(@$mappings) > 1) {
	my $feature_data = $self->_get_feature_data_from_chunks
	    (chunks => $chunks, %args);
	%feats = $self->_create_features_from_feature_data
	    (feature_data => $feature_data, %args);
    }

    # Make alignment object
    my $aln = AT::MapAlignment->new(%feats);
    foreach my $seq (@seqs) {
	$aln->add_seq($seq);
    }

    return $aln;
}


# It could have been worthwile to make the chunks objects (e.g. instances
# of an AT::MapAlignmentFactory::_Chunk class), however I didn't... :)

sub _get_aligned_chunks {
    my ($self, $mappings, $strand) = @_;

    # Initialize
    my @chunks;
    my $num_mappings = scalar(@$mappings);
    my ($tStart, $tEnd);
    #my $strand = $mappings->[0]->strand;
    my @query_info;
    foreach my $m (@$mappings) {
	$tStart = $m->tStart if(!defined($tStart) or $m->tStart < $tStart);
	$tEnd = $m->tEnd if(!defined($tEnd) or $m->tEnd > $tEnd);
	my @hsps = $m->all_HSPs;  # sorted by tStart
	push @query_info, { pos => $hsps[0]->qStart,
			    strand => $m->strand_numeric,
			    HSPnr => 0,
			    HSPs => \@hsps };
    }
    my $tPos = $tStart;

    # This loop generates one chunk in each iteration
    while($tPos <= $tEnd) {
	my $tGap = 0;
	my $chunk_length;
	my @query_status = (Q_GAP) x $num_mappings;

	# Determine the size of this chunk
        #       and the status of each query in this chunk
	for(my $i = 0; $i < $num_mappings; $i++) {
	    my $q = $query_info[$i];

            # Get ref to current HSP
	    # (If we are not currently in an HSP, this will be the next HSP)
	    my $hsp = $q->{HSPs}[$q->{HSPnr}];

	    # Have we reached the end of this mapping?
	    next unless($hsp);

	    # Here comes the tricky stuff...
	    my $my_chunk_length;
	    my $in_hsp =       # Are we inside an HSP?
		($hsp->qStart <= $q->{pos} and $hsp->qEnd >= $q->{pos}) ? 1:0;
	    if($in_hsp) {
		# In HSP
		my $my_tPos = $q->{pos} - $hsp->qStart + $hsp->tStart;
		my $pos_rel_to_target = $my_tPos - $tPos;
		if ($pos_rel_to_target == 0) {
		    # In line with target (possible exon)
		    # Set chunk_length to remaining length of this HSP
		    # (unless we're ahead of another query, then do nothing)
		    unless($tGap) {
			$my_chunk_length = $hsp->qEnd - $q->{pos} + 1;
			$query_status[$i] = Q_MATCH;
		    }
		}
		elsif ($pos_rel_to_target > 0) {
		    # Ahead of target (possible intron)
		    # set chunk_length to distance to target
		    $my_chunk_length = $pos_rel_to_target;
		}
		else {
		    croak("Behind target when inside an HSP!");
		}
	    }
	    else {
		# Not in HSP
		if ($tPos == $hsp->tStart) {
		    # Behind target (there is a target gap)
		    # Set chunk_length to distance to next HSP
		    $my_chunk_length = $hsp->qStart - $q->{pos};
		    unless($tGap) {
			$tGap = 1;
			@query_status = (Q_GAP) x $num_mappings;
		    }
		    $query_status[$i] = Q_TGAP;
		}
		elsif ($tPos < $hsp->tStart) {
		    # In line with target (there is a mismatch)
		    # Set chunk_length to length of mismatch
		    # (unless we're ahead of another query, then do nothing)
		    unless($tGap) {
			my $length1 = $hsp->tStart - $tPos;
			my $length2 = $hsp->qStart - $q->{pos};
			$my_chunk_length =
			    ($length1 < $length2) ? $length1 : $length2;
			$query_status[$i] = Q_MISMATCH;
		    }
		}
		else {
		    croak("Ahead of target but not inside an HSP!");
		}
	    }
	    $chunk_length = $my_chunk_length
		if (defined ($my_chunk_length) and
		    (!defined($chunk_length) or
		     $chunk_length > $my_chunk_length));
	   }

	#print "chunk_length = $chunk_length\n";

	# Create chunk, increment counters
	my %chunk = ( length => $chunk_length,   # Length in bases
		      tGap => $tGap,             # Flag for target gap
		      tStart => undef,           # start in genome
		      qStatus => \@query_status, # status of each cDNA
		      qStarts => []);            # start in each cDNA
	unless($tGap) {
	    $chunk{tStart} = $tPos;
	    $tPos += $chunk_length;
	}
	for(my $i = 0; $i < $num_mappings; $i++) {
	    unless($query_status[$i] == Q_GAP) {
		my $q = $query_info[$i];
		if ($strand == 1) {
		    $chunk{qStarts}->[$i] = $q->{pos};
		}
		else {
		    $chunk{qStarts}->[$i] =
			$mappings->[$i]->qSize - $q->{pos} - $chunk_length + 2;
		}
		$q->{pos} += $chunk_length;
		my $hsp = $q->{HSPs}->[$q->{HSPnr}];
		$q->{HSPnr}++ if ($q->{pos} > $hsp->qEnd); # Move to next HSP?
	    }
	}

	# Add chunk to list
	push @chunks, \%chunk;
    }

    # Make sure chunks are in 5'->3' order with respect to the transcripts
    @chunks = reverse @chunks if ($strand == -1);

    # Return chunks
    return \@chunks;
}


sub _create_gapped_sequences  {
    my ($self, %args) = @_;
    my $tSeqObj = $args{'target_seq'} || croak "No genome_seq";
    my $chunks = $args{'chunks'} || croak "No chunks";
    my $mappings = $args{'mappings'} || croak "No mappings";
    my $num_queries = scalar(@$mappings);
    my $strand = $tSeqObj->strand;
    my $revcom_str = $self->revcom_string;

    # Get query sequence strings.
    # Revcom them if their strand is differenct from the target sequence strand
    my @q_seq_strs;
    foreach my $mapping (@$mappings) {
	my $str = $mapping->query_seq->seq;
	if($mapping->strand_numeric != $strand) {
	    $str = $self->_revcom($str);
	}
	push @q_seq_strs, $str;
    }

    # Variables to hold alignment strings
    my $t_aln_str = "";
    my @q_aln_strs = ("") x $num_queries;

    # Create sequence strings; $c->length bases is added in each iteration
    foreach my $c (@$chunks) {
	if ($args{compact} and $c->{length} > 15 and
	    $self->_all_queries_gap($c)) {
	    # Target seq string
	    my $subseq_start = $c->{tStart} - $tSeqObj->start + 1;
	    my $subseq_end = $subseq_start + $c->{length} - 1;
	    my $subseq =
		$tSeqObj->subseq ($subseq_start, $subseq_start + 4)
		.".....".
		$tSeqObj->subseq ($subseq_end - 4, $subseq_end);
	    $t_aln_str .= ($tSeqObj->strand == 1) ?
		$subseq : $self->_revcom($subseq);
	    # Query seq strings
	    for(my $i = 0; $i < $num_queries; $i++) {
		    $q_aln_strs[$i] = $q_aln_strs[$i]."-----.....-----";
		}
	}
	else {
	    my $gap_str = '-' x $c->{length};
	    # Target seq string
	    if($c->{tGap}) {
		$t_aln_str .= $gap_str;
	    }
	    else {
		my $subseq_start = $c->{tStart} - $tSeqObj->start + 1;
		my $subseq = $tSeqObj->subseq ($subseq_start,
					       $subseq_start + $c->{length}-1);
		$t_aln_str .= ($tSeqObj->strand == 1) ?
		    $subseq : $self->_revcom($subseq);
	    }
	    # Query seq strings
	    for(my $i = 0; $i < $num_queries; $i++) {
		if($c->{qStatus}->[$i] == Q_GAP) {
		    $q_aln_strs[$i] = $q_aln_strs[$i].$gap_str;
		}
		else {
		    my $subseq_start = $c->{qStarts}->[$i] - 1;
		    my $subseq = substr($q_seq_strs[$i], $subseq_start, $c->{length});
#		    my $subseq = $mappings->[$i]->query_seq->subseq
#			($subseq_start, $subseq_start + $c->{length} - 1);
		    $q_aln_strs[$i] = $q_aln_strs[$i].$subseq;
		}
	    }
	}
    }

    # Add genomic context
    my ($chk1, $chk2) = ($tSeqObj->strand == 1) ?
	($chunks->[0], $chunks->[-1]) : ($chunks->[-1], $chunks->[0]);
    my $pos1 = $chk1->{tStart} - $tSeqObj->start;
    my $pos2 = $chk2->{tStart} + $chk2->{length} - $tSeqObj->start + 1;
    my $context1 = ($pos1 >= 1) ? $tSeqObj->subseq(1, $pos1) : "";
    my $context2 = ($pos2 <= $tSeqObj->length) ?
	$tSeqObj->subseq($pos2, $tSeqObj->length) : "";
    my ($fp_context, $tp_context) = ($tSeqObj->strand == 1) ?
	($context1, $context2) :
	($self->_revcom($context2), $self->_revcom($context1));
    my $fp_gapStr = '-' x length($fp_context);
    my $tp_gapStr = '-' x length($tp_context);
    $t_aln_str = $fp_context . $t_aln_str . $tp_context;
    foreach my $q_aln_str (@q_aln_strs) {
        $q_aln_str = $fp_gapStr . $q_aln_str . $tp_gapStr;
    }

    # Create seqs
    my @seqs;
    push @seqs, Bio::LocatableSeq->new   # target seq
	(-id    => $tSeqObj->id .(($strand==1)?"":$revcom_str),
	 -seq   => $t_aln_str,
	 -start => $tSeqObj->start,
	 -end  => $tSeqObj->end,
	 -strand => $tSeqObj->strand );
    for(my $i = 0; $i < $num_queries; $i++) {   # query seq
	my $m = $mappings->[$i];
	my $id = $m->query_seq->accession_number || $m->query_seq->display_id || '';
	$id .= $revcom_str if($m->strand_numeric != $strand);
	push @seqs, Bio::LocatableSeq->new
	    (-id    => $id,
	     -seq   => $q_aln_strs[$i],
	     -start => $m->qStart,
	     -end   => $m->qEnd);
    }

    # Return sequences
    return @seqs;
}


# The feature creation part is a bit messy, but it seems to work fine...
# Instead of inventing a separate data structure for the feature data,
# it could have been better to make the chunk structure more flexible
# and, in the method below, "condense" away small chunks

sub _get_feature_data_from_chunks  {
    my ($self, %args) = @_;

    my $tSeqObj = $args{'target_seq'} || croak "No target_seq";
    my $chunks_ref = $args{'chunks'} || croak "No chunks";
    my $min_intron_length = $args{'min_intron_length'} ||
	DEFAULT_MIN_INTRON_LENGTH;
    my $min_mismatch_length = $args{'min_mismatch_length'} ||
	DEFAULT_MIN_MISMATCH_LENGTH;

    my $strand = $tSeqObj->strand;
    my @chunks = @$chunks_ref;

    if(@{$chunks[0]->{qStatus}} != 1)
    { croak "Features only supported for one query sequence."; }

    my (@feat_data, $f1);

    while(@chunks) {
	my $f2;
	my $c1 = shift @chunks;
	my $stat1 = $c1->{qStatus}->[0];
	if($stat1 == Q_MATCH) {
	    $f2 = { type => Q_MATCH,
		    tStart => $c1->{tStart},
		    qStart => $c1->{qStarts}->[0],
		    tLen => $c1->{length},
		    qLen => $c1->{length}
		};
	}
	else {
	    my $c2 = $chunks[0];
	    my $stat2 = $c2->{qStatus}->[0];
	    if($chunks[0]->{qStatus}->[0] == Q_MATCH) {
		$f2 = { type => $stat1,
			tStart => $c1->{tStart},
			qStart => $c1->{qStarts}->[0],
			tLen => ($stat1 == Q_TGAP) ? 0 : $c1->{length},
			qLen => ($stat1 == Q_GAP) ? 0 : $c1->{length}
		    };
	    }
	    else {
		shift @chunks;
		if($stat1 == Q_TGAP and $stat2 == Q_MISMATCH) {
		    $f2 = { type => Q_MISMATCH,
			    tStart => $c2->{tStart},
			    qStart => $c1->{qStarts}->[0],
			    tLen => $c2->{length},
			    qLen => $c1->{length} + $c2->{length}
			};
		}
		elsif($stat1 == Q_MISMATCH and $stat2 == Q_TGAP) {
		    $f2 = { type => Q_MISMATCH,
			    tStart => $c1->{tStart},
			    qStart => $c1->{qStarts}->[0],
			    tLen => $c1->{length},
			    qLen => $c1->{length} + $c2->{length}
			};
		}
		else {
		    my ($qgap_chunk, $mism_chunk);
		    if($stat1 == Q_GAP and $stat2 == Q_MISMATCH) {
			($qgap_chunk, $mism_chunk) = ($c1, $c2);
		    }
		    elsif($stat1 == Q_MISMATCH and $stat2 == Q_GAP) {
			($qgap_chunk, $mism_chunk) = ($c2, $c1);
		    }
		    else { croak "Illegal query status(es)"; }
		    $f2 = { tStart => (($strand == 1)?$c1:$c2)->{tStart},
			    qStart => $mism_chunk->{qStarts}->[0],
			    tLen => $c1->{length} + $c2->{length},
			    qLen => $mism_chunk->{length} };
		    $f2->{type} =
			($mism_chunk->{length} >= $min_mismatch_length or
			 $qgap_chunk->{length} < $min_intron_length) ?
			 Q_MISMATCH : Q_GAP;
		}
	    }
	    if($chunks[0]->{qStatus}->[0] == Q_MATCH and
	       (($f2->{type} == Q_GAP and $f2->{tLen} < $min_intron_length)
		or
		($f2->{type} & (Q_TGAP|Q_MISMATCH) and
		 $f2->{qLen} < $min_mismatch_length and
		 $f2->{tLen} < $min_mismatch_length)))
	    {
		my $c3 = shift @chunks;
		$f1->{tLen} += $f2->{tLen} + $c3->{length};
		$f1->{qLen} += $f2->{qLen} + $c3->{length};
		$f1->{tStart} = $c3->{tStart} unless($strand == 1);
		$f2 = $f1;
		undef $f1;
	    }
	}
	push @feat_data, $f1 if (defined($f1));
	$f1 = $f2;
    }
    push @feat_data, $f1 if (defined($f1));

    return \@feat_data;
}


sub _create_features_from_feature_data  {
    my ($self, %args) = @_;
    my $feat_data = $args{'feature_data'} || croak "No feature_data";
    my $tSeqObj = $args{'target_seq'} || croak "No target_seq";
    my $mappings = $args{'mappings'} || croak "No mappings";

    my $qSeqObj = $mappings->[0]->query_seq;

    my %feats = (query_exonL => [],
		 target_exonL => [],
		 intronL => [],
		 mismatchL => [],
		 target_gapL => []);

    for(my $i = 0; $i < @$feat_data; $i++) {
	my $f = $feat_data->[$i];

	# Get positions
	my $tAbsStart = $f->{tStart} || 0;
	my $tAbsEnd = $tAbsStart + $f->{tLen} - 1;
	my $tStart = $tAbsStart - $tSeqObj->start + 1;
	my $tEnd = $tAbsEnd - $tSeqObj->start + 1;
	my $qStart = $f->{qStart} || 0;
	my $qEnd = $qStart + $f->{qLen} - 1;

#	my $fuzz;
#	if($f->{type} == Q_GAP and $f->{qLen}) {
#	    $fuzz = $f->{qLen};
#	    $tAbsStart .= '.'.($tAbsStart + $fuzz);
#	    $tAbsEnd = ($tAbsEnd - $fuzz).'.'.$tAbsEnd;
#	    $tStart .= '.'.($tStart + $fuzz);
#	    $tEnd = ($tEnd - $fuzz).'.'.$tEnd;
#	}
#	elsif($i and $feat_data->[$i-1]->{type} == Q_GAP and
#	      $feat_data->[$i-1]->{qLen}) {
#	    $fuzz = $feat_data->[$i-1]->{qLen};
#	    $tAbsStart .= '.'.($tAbsStart - $fuzz);
#	    $tStart .= '.'.($tStart - $fuzz);
#	    $qStart .= '.'.($qStart - $fuzz);
#	}
#	elsif($i+1 < @$feat_data and $feat_data->[$i+1]->{type} == Q_GAP and
#	      $feat_data->[$i+1]->{qLen}) {
#	    $fuzz = $feat_data->[$i+1]->{qLen};
#	    $tAbsEnd .= '.'.($tAbsEnd - $fuzz);
#	    $tEnd .= '.'.($tEnd - $fuzz);
#	    $qEnd .= '.'.($qEnd - $fuzz);
#	}

	# Create feature objects
	my $stat = $f->{type};
	my ($qFeat, $tFeat);
	if ($stat == Q_MATCH) {
	    push @{$feats{target_exonL}},
	    $tFeat = Bio::SeqFeature::Gene::Exon->new();
	    push @{$feats{query_exonL}},
	    $qFeat = Bio::SeqFeature::Gene::Exon->new();
	}
	elsif ($stat == Q_GAP) {
	    push @{$feats{intronL}},
	    $tFeat = Bio::SeqFeature::Generic->new( -primary_tag => "intron" );
	    if($f->{qLen}) {
		push @{$feats{query_mismatchL}},
		$qFeat = Bio::SeqFeature::Generic->new
		    (-primary_tag => "mismatch");
	    }
	}
	elsif ($stat == Q_TGAP) {
	    push @{$feats{target_gapL}},
	    $qFeat = Bio::SeqFeature::Generic->new(-primary_tag=>"target_gap");
	}
	elsif ($stat == Q_MISMATCH) {
	    push @{$feats{target_mismatchL}},
	    $tFeat = Bio::SeqFeature::Generic->new(-primary_tag => "mismatch");
	    push @{$feats{query_mismatchL}},
	    $qFeat = Bio::SeqFeature::Generic->new(-primary_tag => "mismatch");
	}

	# Set feature object data (i.e. positions, seq)
	if($tFeat) {
	    $tFeat->start($tStart); $tFeat->end($tEnd);
	    $tFeat->strand($tSeqObj->strand);
	    $tFeat->add_tag_value('abs_start', $tAbsStart);
	    $tFeat->add_tag_value('abs_end', $tAbsEnd);

#	    my $abs_loc = defined($fuzz) ?
#		Bio::Location::Fuzzy->new() : Bio::Location::Simple->new;
#	    $abs_loc->start($tAbsStart); $abs_loc->end($tAbsEnd);
#	    $abs_loc->strand($tSeqObj->strand);
#	    $tFeat->add_tag_value('abs_location', $abs_loc);

	    $tFeat->attach_seq($tSeqObj);
	}
	if($qFeat) {
	    $qFeat->start($qStart); $qFeat->end($qEnd);
	    $qFeat->attach_seq($qSeqObj);
	}

    }

    return %feats;
}


sub _all_queries_gap {
    my ($self, $chunk) = @_;

    return 0 if ($chunk->{tGap});
    foreach my $query_status (@{$chunk->{qStatus}}) {
	return 0 if($query_status != Q_GAP);
    }
    return 1;
}


sub _revcom {
    my ($self, $str) = @_;

    $str =~
	tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    $str = reverse $str;
    return $str;
}


#sub correct_strand {
#	my ($seq, $strand) = @_;
#	if ($strand == -1 or $strand eq "-")  {
#	    $seq = join('', reverse split('', $seq));
#	    $seq =~ tr/ACGTacgt/TGCAtgca/;
#	}
#	return $seq;
#}


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
	    push @target_exons, _merge_exons(@target_exon_buffer)
		if @target_exon_buffer;
	    push @query_exons,  _merge_exons(@query_exon_buffer)
		if @query_exon_buffer ;
	    @target_exon_buffer = ();
	    @query_exon_buffer  = ();
	    $intron_counter++;
	}
    }
    my @remaining_target_exons = (@target_exon_buffer, @$target_exon_ref);
    my @remaining_query_exons  = (@query_exon_buffer, @$query_exon_ref);

    push @target_exons, _merge_exons(@remaining_target_exons)
	if @remaining_target_exons;
    push @query_exons,  _merge_exons(@remaining_query_exons)
	if @remaining_query_exons ;

    return (\@query_exons, \@target_exons, \@introns);

}

sub _merge_exons  {
    my (@exons) = @_;
    return undef unless @exons;
    # print join(" ", ">>>", (map {$_->start."-". $_->end} @exons), "\n") if @exons>1;
    @exons = sort {$a->start <=> $b->start} @exons;
    my $merged_exon =
	Bio::SeqFeature::Gene::Exon->new(-start=>$exons[0]->start,
					 -end  =>$exons[-1]->end,
					 -strand => $exons[0]->strand);

    if ($exons[0]->has_tag("abs_start") and $exons[-1]->has_tag("abs_end"))  {
	$merged_exon->add_tag_value("abs_start",
				    $exons[0]->each_tag_value("abs_start"));
	$merged_exon->add_tag_value("abs_end",
				    $exons[-1]->each_tag_value("abs_end"));
    }
    $merged_exon->attach_seq($exons[0]->entire_seq);
    return $merged_exon;

}


#sub _get_mapping_sets
#{
#    my ($self, %args) = @_;
#
#    # Get sets of mappings (one alignment is made for each set)
#    my @mapping_sets;
#    if($args{'acc'})
#    {
#	# Check for misuse of args
#	if($args{'mappings'} or $args{'combine'}) {
#	    croak("The use of arg 'acc' together with ".
#		  "args 'mappings' or 'combine' is not supported")
#	}
#
#	my $acc = $args{'acc'};
#	my $query_db = $args{query_db} || $self->{query_db} ||
#	    croak ("No query database set");
#	@mapping_sets = map { [ $_ ] }
#	$query_db->get_mappings_for_acc($acc );#, undef, attach_query_seq => 1);
#    }
#    elsif($args{'accs'}) {
#	my $accsref = $args{'accs'};
#        unless (ref($accsref) eq "ARRAY")  {
#            croak ("accs argument not an array reference");
#        }
#	my $query_db = $args{query_db} || $self->{query_db} ||
#	    croak ("No query database set");
#	my @mappings = map { $query_db->get_mappings_for_acc($_) } @$accsref;
#	if($args{'combine'}) {
#	    @mapping_sets = (\@mappings);   
#	}
#	else {
#	    @mapping_sets = map { [ $_ ] } @mappings;
#	}
#    }
#    elsif($args{'mapping'}) {
#        @mapping_sets = [ $args{'mapping'} ];
#    }
#    elsif($args{'mappings'})  {
#	my $mappingsref = $args{'mappings'};
#        unless (ref($mappingsref) eq "ARRAY")  {
#            croak ("Mappings argument not an array reference");
#        }
#	if($args{'combine'}) {
#	    @mapping_sets = ( [ @$mappingsref ] );
#	}
#	else {
#	    @mapping_sets = map { [ $_ ] } @$mappingsref;
#	}
#    }
#    else  {
#        croak ("No acc, mappings or mapping argument");
#    }
#
#    return \@mapping_sets;
#}




1;
