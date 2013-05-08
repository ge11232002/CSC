# AT::DB::ChainedAlignment module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::DB::ChainedAlignment - interface to MySQL database with 
chained alignments in UCSC schema

=head1 SYNOPSIS

 my $hs_alndb = AT::DB::ChainedAlignment->connect(
    -dbname => "UCSC_hg18",
    -dbhost => "myhostn.mydomain",
    -dbuser => "myusername",
    -dbpass => "mypassword");

 my $hs_gendb = AT::DB::GenomeAssemblyTwoBit->new(...);
 my $mm_gendb = AT::DB::GenomeAssemblyTwoBit->new(...);

 $hs_alndb->target_assembly('hg18', $hs_gendb);
 $hs_alndb->query_assembly('mm8', $mm_gendb);

#get alignments by target region only
 my $alignments = $hs_alndb->get_chained_alignments_for_region(
    target_chr => 'chr22',
    target_start => 14472596,
    target_end => 14473215,
    query_id => 'mm8');

#get alignments by target region, confined to a given query region
 my $alignments = $db->get_alignments_for_regions(
    target_chr => 'chr22',             
    target_start => 14472596,
    target_end => 14473215,
    query_id => 'mm8',
    query_chr => 'chr22',             
    query_start => 14472596,
    query_end => 14473215);

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=head2 connect

 Title     : connect
 Usage     : my $db = $AT::DB::GenomeAlignment->connect(%args);
 Function  : Creates a new object of this class.
 Note      : For details, see docs for superclass AT::DB::MySQLdb.

=cut

package AT::DB::ChainedAlignment;

use strict;
use vars '@ISA';
use Bio::Seq;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use DBI;
use Carp;
use AT::DB::MySQLdb;
use AT::DB::Binner;
use AT::Tools::RangeHandler;


@ISA = qw(AT::DB::MySQLdb AT::DB::Binner);


sub _init
{
    my ($self) = @_;
}


sub target_assembly
{
    my $self = shift;
    if(@_) {
	if(@_ == 1) {
	    $self->{_target_assembly} = shift;
	    $self->{_target_assembly_id} = undef;
	}
	elsif(@_ == 2) {
	    $self->{_target_assembly_id} = shift;
	    $self->{_target_assembly} = shift;
	}
	else {
	    croak "Incorrect number of arguments to target_assembly()";
	}
    }
    return $self->{_target_assembly};
}


sub query_assembly
{
    my $self = shift;
    if(@_ == 1) {
	return $self->{_query_assembly_index}{$_[0]};
    }
    else {
	while(@_) {
	    my $id = shift;
	    my $asm = shift or croak "Incorrect number of arguments to query_assembly()";
	    $self->{_query_assembly_index}{$id} = $asm;
	}
    }
}


=head2 alignment_exists

 Title     : alignment_exists
 Usage     : $db->alignment_exists($query_asm_name);
 Function  : Check whether a particular alignment exists in the database
 Example   : $db->alignment_exists('danRer4');
 Returns   : True if the alignment exists, otherwise false.
 Arguments : Id of query genome assembly

=cut

sub alignment_exists
{
    my ($self, $qName) = @_;

    # Get names of all tables
    my %tables = map { $_ => 1 } $self->_get_table_names();

    # Set first char in query assembly id to uppercase
    $qName = ucfirst($qName);

    # First look for genome-wide chain tables
    my $table1 = 'chain'.$qName;
    my $table2 = $table1.'Link';
    unless($tables{$table1} and $tables{$table2}) {
		# If not found, look for per-chromosome chain tables
		my $target_asm = $self->target_assembly() or croak "Target assembly not set";   
		foreach my $chr ($target_asm->get_chr_names) {
	    	$table1 = $chr.'_chain'.$qName;
	    	$table2 = $table1.'Link';
	    	# Return failure if table not found
	    	return 0 unless($tables{$table1} and $tables{$table2});
		}
    }

    # Return success 
    return 1;
}

sub _get_chain_table_names
{
    my ($self, $tChr, $qName) = @_;

    # Set first char in query assembly id to uppercase
    $qName = ucfirst($qName);

    # First look for genome-wide chain table
    my $table1 = 'chain'.$qName;
    unless($self->_table_exists($table1)) {
	# If not found, assume we have per-chromosome chain table
	$table1 = $tChr.'_chain'.$qName;
    }

    return ($table1, $table1.'Link')
}


=head2 get_alignments_for_regions

 Title     : get_alignments_for_regions
 Usage     : my $alignments = $db->get_alignments_for_regions(
		target_chr => 'chr22',             
		target_start => 14472596,
		target_end => 14473215,
        query_id => 'mm8',
		query_chr => 'chr22',             
		query_start => 14472596,
		query_end => 14473215);
 Function  : Retrieves alignments for the genomic regions
	     in the target and query genomes.
 Returns   : Reference to array of Bio::SimpleAlign objects
 Arguments : strand   Orientation of alignment. Optional.
             max_gap  Max gap in alignment objects. Default 20.

=cut
sub get_alignments_for_regions
{
    my ($self, %args) = @_;

    # Parse arguments
    my $tChr = $args{target_chr} or croak "No target_chr argument";
    my $tStart_bound = $args{target_start} or croak "No target_start argument";
    my $tEnd_bound = $args{target_end} or croak "No target_end argument";
    my $qId = $args{query_id} or croak "No query_id argument";
    my $qChr = $args{query_chr};# or croak "No query_chr argument";
    my $qStart_bound = $args{query_start};# or croak "No query_start argument";
    my $qEnd_bound = $args{query_end};# or croak "No query_end argument";
    my $max_gap = $args{max_gap} || 20;
    my $req_strand = $args{strand};
    #my $confine = $args{confine};

    # Get sequence databases
    my $tSeq_db = $self->target_assembly() or croak "Target assembly not set";
    my $qSeq_db = $self->query_assembly($qId) or croak "Query assembly $qId not set";

    # Get all chain_info for the regions
    my ($chain_table, $chain_link_table) = $self->_get_chain_table_names($tChr, $qId);
    my $chain_info = $self->_get_chain_info_for_regions($chain_table, $tChr, $tStart_bound, $tEnd_bound, 
							$qId, $qChr, $qStart_bound, $qEnd_bound, $req_strand);

#    my $chain_info = $self->_get_chain_info_for_regions($tChr, $tStart_bound, $tEnd_bound, 
#							$qId, $qChr, $qStart_bound, $qEnd_bound, $req_strand);

    # Get alignment data for each chain and make alignments
    my @alignments;
    foreach my $chain (@$chain_info) {
	my ($chain_id, $qStrand, $qSize,$qChr_from_chain) = @$chain;
	# Get alignment data
	my $alndata = $self->_get_alndata_for_chain_part($chain_link_table, $chain_id, $tChr,
							 $tStart_bound, $tEnd_bound, 
							 $qId, $qStart_bound, $qEnd_bound,
							 $qStrand, $qSize, $max_gap);
	
	# Create an alignment for each alignment data entry
	foreach my $a (@$alndata) {    
	    my $min_tStart = $a->{tStart};
	    my $max_tEnd = $a->{tEnd};
	    my $min_qStart = $a->{qStart};
	    my $max_qEnd = $a->{qEnd};
	    
	    #if($confine) {
		$min_tStart = $tStart_bound if($min_tStart < $tStart_bound);
		$max_tEnd = $tEnd_bound if($max_tEnd > $tEnd_bound);

	    if(defined($qStart_bound) and defined($qEnd_bound)){
		$min_qStart = $qStart_bound if($min_qStart < $qStart_bound);
		$max_qEnd = $qEnd_bound if($max_qEnd > $qEnd_bound);
	    }
	    #}
	    my $aln = $self->_create_aln_from_alndata($a,
						      $tSeq_db, $qSeq_db,
						      $tChr, $min_tStart, $max_tEnd,
						      $qChr_from_chain, $min_qStart, $max_qEnd,
		                                      $qSize, $qStrand);
	    push @alignments, $aln if($aln);
	}
    }

    return \@alignments;
}


sub _get_chain_info_for_regions
{
    my ($self, $table, $tChr, $tStart, $tEnd, $qId, $qChr, $qStart, $qEnd, $strand) = @_;

    # Build query
#    my $table = $tChr.'_chain'.ucfirst($qId);
    my $bin_string = $self->bin_restriction_string($tStart, $tEnd);
    my $query = "SELECT id, qStrand, qSize, qName FROM $table
                 WHERE $bin_string AND tName = '$tChr' AND tEnd >= $tStart AND tStart < $tEnd";

    # WHERE $bin_string AND tName = '$tChr' AND qName = '$qChr' AND tEnd >= $tStart AND tStart < $tEnd ";
    if(defined($qChr)){
	$query .= "AND qName = '$qChr'";
    }
    if(defined($qStart) and defined($qEnd)){
	if(!$strand) {
	    $query .= "AND ((qEnd >= $qStart AND qStart < $qEnd AND qStrand = '+') OR
	                (qSize-qStart >= $qStart AND qSize-qEnd < $qEnd AND qStrand = '-'))";
	}
	elsif($strand eq '+') {
	    $query .= "AND qEnd >= $qStart AND qStart < $qEnd AND qStrand = '+'";
	}
	elsif($strand eq '-') {
	    $query .= "AND qSize-qStart >= $qStart AND qSize-qEnd < $qEnd AND qStrand = '-'";
	}
	else {
	    croak "Invalid strand $strand";
	}
    }
    # Run query and get ids
    my $sth = $self->dbh->prepare($query) or die "SQL query failed";
    $sth->execute() or die "SQL query failed";
    my @info;
    while (my $row = $sth->fetchrow_arrayref) {
	push @info, [@$row];
    }
    $sth->finish;

    # Return ids
    return \@info;
}



sub _get_alndata_for_chain_part
{
    my ($self, $table, $chain_id, $tChr, $tStart, $tEnd, $qId, $qStart, $qEnd, $qStrand, $qSize, $max_gap) = @_;

    # Build query
#    my $table = $tChr.'_chain'.ucfirst($qId).'Link';
    my $bin_string = $self->bin_restriction_string($tStart, $tEnd);
#    my $query = "SELECT tStart, tEnd, qStart, qStart+tEnd-tStart as qEnd
#                 FROM $table
#                 WHERE $bin_string AND chainId = ?
#                   AND tEnd >= ? AND tStart < ? 
#                   AND qStart+tEnd-tStart >= ? AND qStart < ?
#                 ORDER BY tStart";

    my $query = "SELECT tStart, tEnd, qStart, qStart+tEnd-tStart as qEnd
                 FROM $table
                 WHERE $bin_string AND chainId = ?
                   AND tEnd >= ? AND tStart < ? ";
    if(defined($qStart) and defined($qEnd)){
	$query .= "AND qStart+tEnd-tStart >= ? AND qStart < ? ";
    }
    $query .= "ORDER BY tStart";
    my $sth = $self->dbh->prepare($query);

    my @alndata;
    my ($prev_block_tEnd, $prev_block_qEnd) = (-$max_gap-1,-$max_gap-1);
    my $i = -1;

    # Run query and put results in array
    if($qStrand eq '+') {
	#$sth->execute($chain_id, $tStart, $tEnd, $qStart, $qEnd) or die "SQL query failed";
	if(defined($qStart) and defined($qEnd)){
	    $sth->execute($chain_id, $tStart, $tEnd, $qStart, $qEnd) or die "SQL query failed";
	}else{
	    $sth->execute($chain_id, $tStart, $tEnd) or die "SQL query failed";
	}
	while (my ($block_tStart, $block_tEnd, $block_qStart, $block_qEnd) = $sth->fetchrow_array) {
	    $i++ if($block_tStart > $prev_block_tEnd+$max_gap or
		    $block_qStart > $prev_block_qEnd+$max_gap);
	    push @{$alndata[$i]->{blocks}}, [$block_tStart+1, $block_qStart+1, $block_tEnd, $block_qEnd];
	    $prev_block_tEnd = $block_tEnd;
	    $prev_block_qEnd = $block_qEnd;
	}
	foreach my $a (@alndata) {
	    $a->{tStart} = $a->{blocks}[0][0];
	    $a->{tEnd} = $a->{blocks}[-1][2];
	    $a->{qStart} = $a->{blocks}[0][1];
	    $a->{qEnd} = $a->{blocks}[-1][3];
	}
    }
    elsif($qStrand eq '-') {
	my $qStart_rev = $qSize - $qEnd + 1;
#	$sth->execute($chain_id, $tStart, $tEnd, $qSize-$qEnd+1, $qSize-$qStart+1) or die "SQL query failed";
	if(defined($qStart) and defined($qEnd)){
	    $sth->execute($chain_id, $tStart, $tEnd, $qSize-$qEnd+1, $qSize-$qStart+1) or die "SQL query failed";
	}else{
	    $sth->execute($chain_id, $tStart, $tEnd) or die "SQL query failed";	
	}

	while (my ($block_tStart, $block_tEnd, $block_qStart, $block_qEnd) = $sth->fetchrow_array) {
	    $i++ if($block_tStart > $prev_block_tEnd+$max_gap or
		    $block_qStart > $prev_block_qEnd+$max_gap);
	    push @{$alndata[$i]->{blocks}}, [$block_tStart+1, $block_qStart+1, $block_tEnd, $block_qEnd];
	    $prev_block_tEnd = $block_tEnd;
	    $prev_block_qEnd = $block_qEnd;
	}
	foreach my $a (@alndata) {
	    $a->{tStart} = $a->{blocks}[0][0];
	    $a->{tEnd} = $a->{blocks}[-1][2];
	    $a->{qStart} = $qSize - $a->{blocks}[-1][3] + 1;
	    $a->{qEnd} = $qSize - $a->{blocks}[0][1] + 1;
	}
    }
    $sth->finish;

    # Return results
    return \@alndata;
}


sub _create_aln_from_alndata
{
    my ($self, $alndata, $tSeq_db, $qSeq_db, $tName, $tStart, $tEnd, $qName, $qStart, $qEnd, $qSize, $strand) = @_;

    # Get sequences
#print STDERR "getting tSeq $tName:$tStart-$tEnd\n";
    my $tSeq = $tSeq_db->get_genome_seq(chr => $tName,
					start => $tStart,
					end => $tEnd);
#print STDERR "getting qSeq $qName:$qStart-$qEnd\n";
    my $qSeq = $qSeq_db->get_genome_seq(chr => $qName,
					start => $qStart,
					end => $qEnd);

#print STDERR "  got seqs\n";

    # Create and return composite alignment
    return $self->_create_alignment($tSeq, $qSeq, $alndata, $strand, $qSize, 0, 0);
}

   
sub _create_alignment
{
    my ($self, $tSeq, $qSeq, $a, $strand, $qSize, $tOffset_major, $qOffset_major) = @_;

    # Get sequence lengths
    my $tLength = $tSeq->length;
    my $qLength = $qSeq->length;

    # Get some attributes into local variables
    my ($tStart, $tEnd) = ($a->{tStart}, $a->{tEnd});
    my ($qStart, $qEnd) = ($a->{qStart}, $a->{qEnd});
    my $tOffset = 1 - $tSeq->start;
    my $qOffset = $strand eq '+' ? 1 - $qSeq->start : $qSeq->end - $qSize;
    #print STDERR "Offsets: $tOffset $qOffset\n";

    # Get coords for ungapped blocks
    my $all_blocks = $a->{blocks};
    my @blocks;
    foreach my $block (@$all_blocks) {
	my ($s1, $s2, $e1, $e2) =
	    $self->_restrict_block($block->[0] + $tOffset, $block->[1] + $qOffset,
				   $block->[2] + $tOffset, $block->[3] + $qOffset,
				   $tLength, $qLength);
	#print "bef: $block->[0]-$block->[2] $block->[1]-$block->[3]\n";
	#print "aft: $s1-$e1 $s2-$e2\n";
	next if($s1 > $e1);
	push @blocks, ['', $s1, $s2, $e1, $e2]; # leave first empty for compatability w Malin's subs
    }
    #print STDERR "made ", scalar @blocks, " blocks from ",scalar @$all_blocks," alndata items\n";
    return unless(@blocks);

    # Create gapped sequences
    my ($alnseq1, $alnseq2) = $self->_create_gapped_seqs($tSeq, $qSeq, \@blocks, $strand);
    return unless($alnseq2);

    # Get start and end pos for alignment
    my @pos_seq1=($blocks[0][1] + $tOffset_major,
		  $blocks[-1][3] + $tOffset_major);
    my @pos_seq2=($blocks[0][2] + $qOffset_major,
		  $blocks[-1][4] + $qOffset_major);

    # Create SimpleAlign object
    my $simple_align_obj = $self->_create_subalignment_from_gapped_seqs
	($alnseq1, $alnseq2, $tSeq->id, $qSeq->id, \@pos_seq1, \@pos_seq2, $strand);

    # Return the alignment object
    return $simple_align_obj;
}


sub _create_gapped_seqs
{
  my ($self, $seq1, $seq2, $ungapped_blocks, $strand) = @_;

  #my $subseq1 = $seq1->subseq($tStart - $tSeq_start + 1, $tEnd - $tSeq_start + 1);
  #my $subseq2 = $seq2->subseq($qStart - $qSeq_start + 1, $qEnd - $qSeq_start + 1);

  # yes, we could save time by less revcoming
  # ie first get subseqs. then we need offsets o1 en o2 to subtract from all coords

  #If the second sequence was reverse complemented in the alignment by blastz, 
  #the positions above are calculaled on the reverse compemented sequence. The 
  #alignamble sequence should correspond to the reverse complementation, therefore
  #the reverse complementation of the original sequence will be used as a template.
  if($strand eq '-'){
    $seq2=$seq2->revcom;
  }

  #Creating two sequences of equal length that can be aligned using SimpleAlign
  #First ungapped block:
  my $alseq1=$seq1->subseq($ungapped_blocks->[0][1],$ungapped_blocks->[0][3]);
  my $alseq2=$seq2->subseq($ungapped_blocks->[0][2],$ungapped_blocks->[0][4]);

#print STDERR join("\t", @{$ungapped_blocks->[0]}), "\n";

  for my $i(1  .. (scalar (@$ungapped_blocks)-1)){
    
#print STDERR join("\t", @{$ungapped_blocks->[$i]}), "\n";
    
    #adding  regions in between ungapped blocks:
    if ($ungapped_blocks->[$i][1]>($ungapped_blocks->[($i-1)][3]+1)){
      $alseq1=$alseq1.$seq1->subseq($ungapped_blocks->[($i-1)][3]+1,$ungapped_blocks->[$i][1]-1);
    }
    if ($ungapped_blocks->[$i][2]>($ungapped_blocks->[($i-1)][4]+1)){
      $alseq2=$alseq2.$seq2->subseq($ungapped_blocks->[($i-1)][4]+1,$ungapped_blocks->[$i][2]-1);
    }
    my $difference=(($ungapped_blocks->[($i)][2]-$ungapped_blocks->[($i-1)][4])-($ungapped_blocks->[$i][1]-$ungapped_blocks->[($i-1)][3]));
    if($difference>0){
      $alseq1=$alseq1.'-' x $difference;
      $alseq2=$alseq2;
    }
    elsif($difference<0){
      $alseq2=$alseq2.'-' x (-$difference);
      $alseq1=$alseq1;
    }
    #adds aligned regions
  if($ungapped_blocks->[$i][1]<=$ungapped_blocks->[$i][3]){
      $alseq1=$alseq1.$seq1->subseq($ungapped_blocks->[$i][1],$ungapped_blocks->[$i][3]);
    }
    #In case of strange (wrong?) alignments where the start position of an ungapped block is larger than the end position, the region between the two adjacent ungapped blocks is added.(This problem has ocurred for at least one alignment, perhaps t is a bug in BlastZ???...)
    else{
      warn "Errors in the input alignment";
      return;
    }
    
    if($ungapped_blocks->[$i][2]<=$ungapped_blocks->[$i][4]){
      $alseq2=$alseq2.$seq2->subseq($ungapped_blocks->[$i][2],$ungapped_blocks->[$i][4]);
    }
    #In case of strange alignment where start position is larger than the end position
    #for an ungapped block.
    else{
      warn "Errors in the input alignment";
      return;
    }
  }
  return($alseq1, $alseq2);
}



sub _create_subalignment_from_gapped_seqs {
  my ($self, $locatableseq1, $locatableseq2, $id1, $id2, $pos_seq1, $pos_seq2, $strand) = @_;

  #LocatableSeq objects are used to make an alignment:
  $locatableseq1 = Bio::LocatableSeq -> new(-seq => $locatableseq1, -id => $id1, 
					    -start =>$pos_seq1->[0], -end => $pos_seq1->[1],
					    -strand => 1);
  $locatableseq2 = Bio::LocatableSeq -> new(-seq =>$locatableseq2, -id => $id2, 
					    -start=>$pos_seq2->[0], -end =>$pos_seq2->[1],
					    -strand => $strand eq '+' ? 1 : -1);
  
  #An alignment object is created:
  my $AlignObj= Bio::SimpleAlign->new();
  $AlignObj->add_seq($locatableseq1);
  $AlignObj->add_seq($locatableseq2);
  $AlignObj->gap_char('-');
  return ($AlignObj);
}
    

sub _restrict_block {
    my ($self, $s1, $s2, $e1, $e2, $l1, $l2) = @_;

    if($s1 < 1) {
	$s2 += (1-$s1);
	$s1 = 1;
    }
    if($s2 < 1) {
	$s1 += (1-$s2);
	$s2 = 1;
    }
    if($e1 > $l1) {
	$e2 -= ($e1-$l1);
	$e1 = $l1;
    }
    if($e2 > $l2) {
	$e1 -= ($e2-$l2);
	$e2 = $l2;
    }

    return ($s1, $s2, $e1, $e2);
}


sub _distance {
    my ($self, $s1, $e1, $s2, $e2) = @_;
    if($s1 < $s2) {
	my $d = $s2 - $e1;
	return $d > 0 ? $d : 0;
    }
    else {
	my $d = $s1 - $e2;
	return $d > 0 ? $d : 0;
    }
}


1;


