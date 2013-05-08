# AT::DB::GenomeAlignment module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::DB::GenomeAlignment - interface to MySQL database of
whole-genome alignments

=head1 SYNOPSIS

 my $hs_alndb = AT::DB::GenomeAligment->connect(
    -dbname => "AT_HS_JUL03",
    -dbhost => "myhostn.mydomain",
    -dbuser => "myusername",
    -dbpass => "mypassword");

 my $hs_gendb = AT::DB::GenomeAssembly->new(...);
 my $mm_gendb = AT::DB::GenomeAssembly->new(...);
 
 my $alnset = $hs_alndb->get_alignments_for_region(
    chr => 'chr22',             
    start => 14472596,
    end => 14473215,
    seq_dbs => [$hs_gendb, $mm_gendb],
    aln_type => 'net',
    confine => 1);

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=head2 connect

 Title     : connect
 Usage     : my $db = $AT::DB::GenomeAlignment->connect(%args);
 Function  : Creates a new object of this class.
 Note      : For details, see docs for superclass AT::DB::MySQLdb.

=cut

package AT::DB::GenomeAlignment;

use strict;
use vars '@ISA';
use Bio::Seq;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use DBI;
use Carp;
use AT::DB::MySQLdb;
use AT::DB::Binner;
use AT::DB::Importer;
use AT::Mapping;
use AT::mRNAInfo;
use AT::HSP;
use AT::BlatpslxIO;
use AT::AlignmentSet;
use AT::Alignment::Composite;
use AT::Alignment::Subalignment;
use AT::Tools::RangeHandler;


@ISA = qw(AT::DB::MySQLdb AT::DB::Binner);


my @_WGA_COLUMN_NAMES = qw(aln_id
			   bin	   score   strand
			   tName   tStart  tEnd
			   qName   qStart  qEnd
			   tStarts qStarts blockSizes);
my $_WGA_COLUMN_NAMES = join(',',@_WGA_COLUMN_NAMES);


sub _init
{
    my ($self) = @_;
}


sub _wga_column_values {
    my ($self, $aln, $qGenome_db) = @_;

    unless ($aln->isa("Bio::SimpleAlign")) { 
	croak "Wrong type of alignment object";
    }

    my ($tSeq, $qSeq) = $aln->each_seq;
    my $tName = $tSeq->id;
    my $tStart = $tSeq->start;
    my $tEnd = $tSeq->end;
    my $qName = $qSeq->id;
    my $qStart = $qSeq->start;
    my $qEnd = $qSeq->end;
    my $strand = $qSeq->strand == 1 ? '+' : '-';
    my $qSize;
    if($qName =~ /_rev$/) {
	$qName =~ s/_rev$//;
	$qSize = $qGenome_db->get_chr_size($qName);
	($qStart, $qEnd) = ($qSize - $qEnd + 1, $qSize - $qStart + 1);	
    }
    my $score = $aln->score();
    my $bin = $self->bin_from_coord_range($tStart, $tEnd);
    
    my @tStarts;
    my @qStarts;
    my @blockSizes;

    my $tSeqStr = $tSeq->seq;
    my $qSeqStr = $qSeq->seq;

    # Trim gaps at beginning or end of alignments (this should be rare, so the loops need not be optimized for speed)
    while(substr($tSeqStr,0,1) eq '-' or substr($qSeqStr,0,1) eq '-') {
	if(substr($tSeqStr,0,1,"") ne '-') { $tStart++ }
	if(substr($qSeqStr,0,1,"") ne '-') { if($strand eq '+') { $qStart++ } else { $qEnd-- } };
    }
    while(substr($tSeqStr,-1,1) eq '-' or substr($qSeqStr,-1,1) eq '-') {
	if(substr($tSeqStr,-1,1,"") ne '-') { $tEnd-- }
	if(substr($qSeqStr,-1,1,"") ne '-') { if($strand eq '+') { $qEnd-- } else { $qStart++ } };
    }

    my @tBlockLen = map {length($_)} split(/(\-+)/, $tSeqStr);
    my @qBlockLen = map {length($_)} split(/(\-+)/, $qSeqStr);
    #print join(' ', @tBlockLen), "\n";
    #print join(' ', @qBlockLen), "\n";
    my ($i, $j) = (0,0);
    my ($tRelStart, $qRelStart) = (1,1);
    while($tBlockLen[$i] != $qBlockLen[$j]) {
	#print "$tBlockLen[$i] $qBlockLen[$j] ";
	if($tBlockLen[$i] < $qBlockLen[$j]) {
	    # aligned block followed by target gap
	    my $block_size = $tBlockLen[$i];  # size of aligned block
	    push @tStarts, $tRelStart;
	    push @qStarts, $qRelStart;
	    push @blockSizes, $block_size;
	    $i++;
	    if($i >= @tBlockLen) {
		warn("Gapped sequences $tName:$tStart-$tEnd, $qName:$qStart-$qEnd have different lengths");
		return;
	    }
	    my $gap_size = $tBlockLen[$i];  # size of gap
	    $i++;
	    $tRelStart += $block_size;
	    $qRelStart += $block_size + $gap_size;
	    $qBlockLen[$j] -= $block_size + $gap_size;
	}
	else {
	    # aligned block followed by query gap
	    my $block_size = $qBlockLen[$j];  # size of aligned block
	    push @tStarts, $tRelStart;
	    push @qStarts, $qRelStart;
	    push @blockSizes, $block_size;
	    $j++;
	    if($j >= @qBlockLen) {
		warn("Gapped sequences $tName:$tStart-$tEnd, $qName:$qStart-$qEnd have different lengths");
		return;
	    }
	    my $gap_size = $qBlockLen[$j];  # size of gap
	    $j++;
	    $qRelStart += $block_size;
	    $tRelStart += $block_size + $gap_size;
	    $tBlockLen[$i] -= $block_size + $gap_size;
	}
        #print "$i $j $tRelStart $qRelStart\n";
    }
    # last block of alignment
    unless($i == @tBlockLen-1 and $j == @qBlockLen-1) {
	warn("Gapped sequences $tName:$tStart-$tEnd, $qName:$qStart-$qEnd have different lengths");
	return;
    }
    push @tStarts, $tRelStart;
    push @qStarts, $qRelStart;
    push @blockSizes, $tBlockLen[$i];

    my @values =   (undef,
		    $bin,
		    $score,
		    $strand,
		    $tName,
		    $tStart,
		    $tEnd,
		    $qName,
		    $qStart,
		    $qEnd,
		    join(',', @tStarts),
		    join(',', @qStarts),
		    join(',', @blockSizes));

    return \@values;
}


sub _create_wga_table
{
    my ($self, $table_name) = @_;
    my $sql = qq!
CREATE TABLE $table_name (
  aln_id mediumint(8) unsigned NOT NULL auto_increment,
  bin smallint(5) unsigned NOT NULL default '0',
  score mediumint(6) NOT NULL default '0',
  strand char(2) NOT NULL default '',
  tName varchar(40) NOT NULL default '',
  tStart int(10) unsigned NOT NULL default '0',
  tEnd int(10) unsigned NOT NULL default '0',
  qName varchar(40) NOT NULL default '',
  qStart int(10) unsigned NOT NULL default '0',
  qEnd int(10) unsigned NOT NULL default '0',
  tStarts text NOT NULL default '',
  qStarts text NOT NULL default '',
  blockSizes text NOT NULL default '',
  PRIMARY KEY(aln_id)
);!;

    my $rv = $self->dbh->do($sql);
    return defined($rv) ? 1 : undef;
}


=head2 create_alignment_importer

 Title     : create_aln_importer
 Usage     : my $importer = $db->create_alignment_importer($qGenome_db, 'net');
 Function  : Creates an AT::DB::Importer object that allows efficient
             import of a large number of alignments to the database.
             See the AT::DB::Importer module for further details.
 Example   : my $importer = $db->create_alignment_importer($qGenome_db, 'net');
             while (my $aln = $aln_stream->next_aln) {
		 $importer->store_object($aln);
	     }
             my @result = $importer->finish;
 Returns   : An AT::DB::Importer object
 Arguments : query_genome_db  AT genome assembly object for query genome
                              (i.e. an object of class AT::DB::GenomeAssembly
                              or AT::DB::GenomeAssemblyNibs)
             aln_type         Tag for alignment type (e.g. 'net' or 'tight')


=cut


sub create_alignment_importer {
    my ($self, $qGenome_db, $aln_type) = @_;
    
    # Figure out table name and create table if neccessary
    my $table = 'WGA_'.($qGenome_db->id).'_'.uc($aln_type);
    if($self->_table_exists($table)) {
	warn "Importing alignments into existing table $table\n";
    }
    else {
        $self->_create_wga_table($table) or die "Could not create table $table";
    }

    # define wrapper for storing one alignment
    my $storing_wrapper = sub {
	my ($importer, $aln) = @_;
	my $row = $self->_wga_column_values($aln, $qGenome_db);
	$importer->store_record($table, $row);
   };

    # prepare table info needed by importer
    my $tables_columns = {
	$table => [ @_WGA_COLUMN_NAMES ]
    };

    # create importer object
    my $importer = AT::DB::Importer->new( db => $self,
					  tables_columns => $tables_columns,
					  storing_wrapper => $storing_wrapper);

    # return the creation
    return $importer;
}


=head2 index_genome_alignment

 Title     : index_genome_alignment
 Usage     : $db->index_genome_alignment($query_asm_name, $aln_type);
 Function  : Indexes a whole-genome alignment stored in the database
 Example   : $db->index_genome_alignment('MM_MAY04', 'net');
 Returns   : -
 Arguments : query_asm_name  Name of query genome assembly
             aln_type        Tag for alignment type (e.g. 'net' or 'tight')

=cut


sub index_genome_alignment
{
    my ($self, $qName, $aln_type) = @_;

    # Get table name and check that table exists
    my $table = 'WGA_'.$qName.'_'.uc($aln_type);
    unless($self->_table_exists($table)) {
	die "Table $table does not exist\n";
    }
    
    # here: should check if there are any indices on the table already

    # Create indices
    $self->_do_or_die("alter table $table add unique chr_pos(tName, bin,tStart,tEnd)");
    $self->_do_or_die("alter table $table add unique chr_start(tName,tStart)");
    $self->_do_or_die("alter table $table add unique chr_end(tName,tEnd)");
}


=head2 alignment_exists

 Title     : alignment_exists
 Usage     : $db->alignment_exists($query_asm_name, $aln_type);
 Function  : Check whether a particular alignment exists in the database
 Example   : $db->alignment_exists('MM_MAY04', 'net');
 Returns   : True if the alignment exists, otherwise false.
 Arguments : query_asm_name  Name of query genome assembly
             aln_type        Tag for alignment type (e.g. 'net' or 'tight')

=cut

sub alignment_exists
{
    my ($self, $qName, $aln_type) = @_;
    my $table = 'WGA_'.$qName.'_'.uc($aln_type);
    return $self->_table_exists($table) ? 1 : 0;
}


=head2 get_alignments_of_seqs

 Title     : get_alignments_of_seqs
 Usage     : my $aln_set = $db->get_alignments_of_seqs
		(seqs => [$seq1, $seq2],
		 aln_type => 'net' );
 Function  : Retrieves alignments between two genomic
	     sequences.
 Returns   : AT::AlignmentSet, containing at most two  alignments
	     (one with $seq2 in forward orientation and one with
	     $seq2 reversed).
 Arguments : seqs  Sequences to align. Obtain using an
		   AT::DB::GenomeAssembly* object.
		   The first sequence should be from the target
		   assembly for which the alignment db was made.
	     aln_type Type of alignment to retrieve
	              (usually 'net' or 'tight')
		   
=cut

sub get_alignments_of_seqs
{
    my ($self, %args) = @_;

    # Get aln_type from arguments
    my $aln_type = $args{aln_type} or croak "No aln_type argument";

    # Get sequences from argument
    croak "No seqs argument" unless($args{seqs});
    my ($tSeq, $qSeq) = @{$args{seqs}};
    my ($tChr) = $tSeq->id =~ /(chr.*)/;
    my ($qChr) = $qSeq->id =~ /(chr.*)/;
    my ($tSeq_start, $tSeq_end) = ($tSeq->start, $tSeq->end);
    my ($qSeq_start, $qSeq_end) = ($qSeq->start, $qSeq->end);
    my ($query_id) = substr($qSeq->id,0,8);

    # Get alignment data
    my $alndata = $self->_get_alndata_for_regions($query_id, $aln_type,
						  $tChr, $tSeq_start, $tSeq_end,
						  $qChr, $qSeq_start, $qSeq_end);

    # Separate + and - alignment data
    my (@plus_alndata, @minus_alndata);
    foreach my $a (@$alndata) {
	if($a->{strand} eq '+') { push @plus_alndata, $a; }
	elsif($a->{strand} eq '-') { push @minus_alndata, $a; }
	else {
	    warn("Illegal strand ", $a->strand, " for WGA alignemnt", $a->aln_id);
	    return undef;
        }
    }
    
    # Create composite aligments, one for each strand
    my $plus_aln = $self->_create_aln_from_seqs_and_alndata($tSeq, $qSeq, \@plus_alndata);
    my $minus_aln = $self->_create_aln_from_seqs_and_alndata($tSeq, $qSeq, \@minus_alndata);

    # Create alignmentset
    my $alignmentset = AT::AlignmentSet->new(-seq1 => $tSeq, -seq2 => $qSeq);
    $alignmentset->add_alignment(-alignment => $plus_aln) if($plus_aln);
    $alignmentset->add_alignment(-alignment => $minus_aln) if($minus_aln);

    return $alignmentset;
}


=head2 get_alignments_for_region

 Title     : get_alignments_for_region
 Usage     : my $alnset = $db->get_alignments_for_region(
		chr => 'chr22',             
		start => 14472596,
		end => 14473215,
		seq_dbs => [$target_gendb, $query_gendb],
		aln_type => 'net',
		confine => 1);
 Function  : Retrieves alignments for a given genomic region
	     in the target genome.
 Returns   : AT::AlignmentSet
 Arguments : chr, start, end  Region in target genome
             seq_dbs  Genome sequence databases. Target first,
		      query second.
	     aln_type  Type of alignment to retrieve
	               (usually 'net' or 'tight')
	     confine  Set to 1 to slice alignments to
	              the given region. Optional.
		      
=cut

sub get_alignments_for_region
{
# args: chr, start, end
#       seqdb1, seqdb2
#       aln_type
# to be added: confinement args
    my ($self, %args) = @_;

    my $max_gap_in_composite_aln = 100000;

    # Parse arguments
    my $tChr = $args{chr} or croak "No chr argument";
    my $tStart_bound = $args{start} or croak "No start argument";
    my $tEnd_bound = $args{end} or croak "No end argument";
    my $aln_type = $args{aln_type} or croak "No aln_type argument";
    croak "No seq_dbs argument" unless($args{seq_dbs});
    my ($tSeq_db, $qSeq_db) = @{$args{seq_dbs}};
    my $confine = $args{confine};

    # Get alignment data
    my $alndata = $self->_get_alndata_for_region($qSeq_db->id, $aln_type,
						  $tChr,
						  $tStart_bound,
						  $tEnd_bound);

    # Create alignmentset
    my $alignmentset = AT::AlignmentSet->new();

    # Separate alignment data into sets and and create a composite  alignment
    # for each set.
    if(@$alndata) {
        my @alndata_set = ($alndata->[0]);
	my $min_tStart = $alndata->[0]->{tStart};
        my $max_tEnd = $alndata->[0]->{tEnd};
	my $min_qStart = $alndata->[0]->{qStart};
        my $max_qEnd = $alndata->[0]->{qEnd};
        foreach my $a (@$alndata[1..@$alndata-1]) {
	    if($a->{tStart} <= $max_tEnd+$max_gap_in_composite_aln and
	       $a->{qName} eq $alndata_set[0]->{qName} and
	       $a->{strand} eq $alndata_set[0]->{strand} and
	       $self->_distance($min_qStart, $max_qEnd, $a->{qStart}, $a->{qEnd}) < $max_gap_in_composite_aln) {
		push @alndata_set, $a;
		$max_tEnd = $a->{tEnd} if($a->{tEnd} > $max_tEnd);
		$min_qStart = $a->{qStart} if($a->{qStart} < $min_qStart);
		$max_qEnd = $a->{qEnd} if($a->{qEnd} > $max_qEnd);
	    }
	    else {
		$min_tStart = $tStart_bound if($confine and $min_tStart < $tStart_bound);
		$max_tEnd = $tEnd_bound if($confine and $max_tEnd > $tEnd_bound);
	        my $aln = $self->_create_aln_from_alndata(\@alndata_set,
	    					      $tSeq_db,
						      $qSeq_db,
						      $min_tStart,
						      $max_tEnd,
						      $min_qStart,
						      $max_qEnd);
	        $alignmentset->add_alignment(-alignment => $aln) if($aln);
	        @alndata_set = ($a);
		$min_tStart = $a->{tStart};
	        $max_tEnd = $a->{tEnd};
		$min_qStart = $a->{qStart};
	        $max_qEnd = $a->{qEnd};
	    }
	}
	$min_tStart = $tStart_bound if($confine and $min_tStart < $tStart_bound);
	$max_tEnd = $tEnd_bound if($confine and $max_tEnd > $tEnd_bound);
        my $aln = $self->_create_aln_from_alndata(\@alndata_set,  # create for last set
					      $tSeq_db,
					      $qSeq_db,
					      $min_tStart,
					      $max_tEnd,
					      $min_qStart,
					      $max_qEnd);
        $alignmentset->add_alignment(-alignment => $aln) if($aln);
    }
     
    return $alignmentset;
}


=head2 count_aligned_bp_for_region

 Title     : count_aligned_bp_for_region
 Usage     : my $alnset = $db->count_aligned_bp_for_region(
		chr => 'chr22',             
		start => 14472596,
		end => 14473215,
		query_assembly_name => "HS_JUL03",
		aln_type => 'net');
 Function  : Counts the number of aligned bp for a region of the
	     target genome. Matching and mismatching bp are counted,
	     but not bp in gaps.
 Returns   : Integer
 Arguments : chr, start, end   Region in target genome
	     query_assembly_name, aln_type   What alignment to use
		      
=cut


sub count_aligned_bp_for_region
{
    my ($self, %args) = @_;

    # Parse arguments
    my $chr = $args{chr} or croak "No chr argument";
    my $rgn_start = $args{start} or croak "No start argument";
    my $rgn_end = $args{end} or croak "No end argument";
    my $aln_type = $args{aln_type} or croak "No aln_type argument";
    my $query_assembly = $args{query_assembly_name} or croak "No query_assembly_name argument";

    # Get alignment data
    my $alndata = $self->_get_alndata_for_region($query_assembly, $aln_type,
						  $chr,
						  $rgn_start,
						  $rgn_end);

    my $aligned_bp = 0;
    
    foreach my $a (@$alndata) {
	my $block_starts = $a->{tStarts};
	my $block_sizes = $a->{blockSizes};
	my $offset = $a->{tStart}-1;
	my ($rel_rgn_start, $rel_rgn_end);
	for my $i (0..@$block_starts-1) {
	    my $block_start = $block_starts->[$i]+$offset;
	    last if $block_start > $rgn_end;
	    my $block_end = $block_start + $block_sizes->[$i] - 1;
	    my $ol_start = $rgn_start > $block_start ? $rgn_start : $block_start;
	    my $ol_end = $rgn_end < $block_end ? $rgn_end : $block_end;
	    my $ol = $ol_end - $ol_start + 1;
	    $aligned_bp += $ol if($ol > 0);

	}
    }

    return $aligned_bp;
}


=head2 get_grouped_segment_pairs_for_split_region

 Title     : get_grouped_segment_pairs_for_split_region
 Usage     : my $seg_pair_groups =
               $db->get_grouped_segment_pairs_for_split_region
	       (chr => 'chr2',
	        split_region => [[90929683, 90930047],
                                 [90930606, 90930944],
                                 [90935913, 90936071],
                                 [90937019, 90937117],
                                 [90937199, 90937321],
                                 [90937750, 90937803],
                                 [90938899, 90939098],
                                 [90939505, 90939788]],
		query_assembly_name => "HS_MAY04",
		aln_type => 'net');
	     # Above example is for Rapsn gene in MM_MAY04 assembly.
	     # Note that the coords in split_region are sorted
	     # in increasing order.
 Function  : Gets aligned segment pairs between two assemblies.
             Segment pairs are retrieved for the entire region
	     spanned by split_region, including intervening subregions.
	     The segment pairs are grouped by their position in
	     the target genome, so that each group should represent
	     one distinct locus in the target genome.
	     The groups are sorted, in decreasing order, by how much
	     they span the split region. The idea is that the most
	     "relevant" group should be the first in the returned
	     array.
	     The method is is implemented to serve an immediate
	     need. It can certainly be made more generic and
	     object-oriented in the future.
 Returns   : reference to array of arrays of arrays (3D array)
             Each element in the primary array represents one group.
	     Each element in a secondary array represents one segment
	     pair.
	     Each segment pair array has the following elements:
	      target chr, target start, target end,
	      query chr, query start, query end, strand
	     (Yes, this could be 'objectified'!)
 Arguments : chr           Chromosome in target genome
             split_region  2D array of start & end coords in target
	                   genome. Typically exon coords.
			   The coords should be sorted in increasing
			   order.
	     start         Start in target genome (optional)
	     end           End in target genome (optional)
	                   If start, end are not given, the bounds
			   of the split_region will be used.
	     query_assembly_name, aln_type   What alignment to use
		      
=cut

sub get_grouped_segment_pairs_for_split_region
{
    my ($self, %args) = @_;

    # Parse arguments
    my $tChr = $args{chr} or croak "No chr argument";
    my $split_region = $args{split_region} or croak "No split_region argument";

    my $aln_type = $args{aln_type} or croak "No aln_type argument";
    my $query_assembly = $args{query_assembly_name} or croak "No query_assembly_name argument";

    # Find start and end of the split region; assume blocks are sorted in incresing order
    my $tStart_bound =
    	($args{start} or $split_region->[0][0]) or
	 croak "Invalid split_region";
    my $tEnd_bound =
    	($args{end} or $split_region->[-1][1]) or
	croak "Invalid split_region";

    # Get alignment data
    my $alndata = $self->_get_alndata_for_region($query_assembly,
						 $aln_type,
						 $tChr,
						 $tStart_bound,
						 $tEnd_bound);

    # Group alignment data
    # Each group should represent a distinct locus & strand in the secondary genome
    my $max_d_in_group = ($tEnd_bound - $tStart_bound + 1) * 3;  # arbitarily set to 3 x size of primary region
    my $alndata_groups = $self->_group_alndata($alndata, $max_d_in_group);

    # Sort group by how much they span the split region in the primary genome
    $alndata_groups = $self->_sort_alndata_groups_by_span($alndata_groups, $split_region);

    my @segment_pair_groups =
	map { $self->_get_segment_pairs_for_alndata_group($_) }
	@$alndata_groups;

    return \@segment_pair_groups;    
}


sub get_first_aligned_region_after_pos
{
    my ($self, %args) = @_;

    # Parse arguments
    my $chr = $args{chr} or croak "No chr argument";
    my $pos = $args{pos} or croak "No pos argument";
    my $aln_type = $args{aln_type} or croak "No aln_type argument";
    my $query_assembly = $args{query_assembly_name} or croak "No query_assembly_name argument";

    # Get and return region
    my $alndata = $self->_get_first_alndata_after_pos($query_assembly, $aln_type, $chr, $pos);
    return $alndata ? $self->_get_region_from_alndata($alndata) : ();
}


sub get_first_aligned_region_before_pos
{
    my ($self, %args) = @_;

    # Parse arguments
    my $chr = $args{chr} or croak "No chr argument";
    my $pos = $args{pos} or croak "No pos argument";
    my $aln_type = $args{aln_type} or croak "No aln_type argument";
    my $query_assembly = $args{query_assembly_name} or croak "No query_assembly_name argument";

    # Get and return region
    my $alndata = $self->_get_first_alndata_before_pos($query_assembly, $aln_type, $chr, $pos);
    return $alndata ? $self->_get_region_from_alndata($alndata) : ();
}


sub _get_region_from_alndata
{
    my ($self, $alndata) = @_;
    return ($alndata->{tName},
	    $alndata->{tStart},
	    $alndata->{tEnd},
	    $alndata->{qName},
	    $alndata->{qStart},
	    $alndata->{qEnd},
	    $alndata->{strand});
}


sub _group_alndata
{
    my ($self, $alndata_ref, $max_d) = @_;
    
    # Sort by chr, strand, end in secondary genome
    my @alndata = sort { $a->{qName}  cmp $b->{qName} or
		         $a->{strand} cmp $b->{strand} or
		         $a->{qEnd} <=> $b->{qEnd} } @$alndata_ref;
		
    my @groups;      
    my (@group, $group_strand, $group_end);
    my $group_chr = "";
    foreach my $a (@alndata) {
	if($a->{qName} eq $group_chr and $a->{strand} eq $group_strand and $a->{qStart} < $group_end+$max_d) {
	    push @group, $a;
	    $group_end = $a->{qEnd} if($group_end < $a->{qEnd});
	}
	else {
	    push @groups, [@group] if(@group);
	    $group_chr = $a->{qName};
	    $group_strand = $a->{strand};
	    $group_end = $a->{qEnd};
	    @group = ($a);
        }
    }
    push @groups, [@group] if(@group);

#    #debug output
#    foreach my $g (@groups) {
#	print STDERR "Group\n";
#	foreach my $aln (@$g) {
#	    my $seq2 = $aln->get_seq2;
#            print STDERR $seq2->id, " ", $aln->get_strand, " ", $seq2->start, "-", $seq2->end, "\n";
#	}
#    }
    
    return \@groups;
}


sub _sort_alndata_groups_by_span
{
    my ($self, $groups, $range_to_span) = @_;

    my %span;    
    foreach my $group (@$groups) {
        my @aligned_range = map { [ $_->{tStart}, $_->{tEnd} ] } @$group;
        my $intersection =
	    AT::Tools::RangeHandler->compute_intersection($range_to_span, \@aligned_range);
	$span{$group} = AT::Tools::RangeHandler->sum($intersection);
    }

    my @sorted_groups = sort { $span{$b} <=> $span{$a} } @$groups;

    return \@sorted_groups;
}


sub _get_segment_pairs_for_alndata_group
{
    my ($self, $alndata_group_ref) = @_;

    my @alndata_group = sort { $a->{tStart} <=> $b->{tStart} } @$alndata_group_ref;

    my @segment_pairs = map {
	[ $_->{tName},
	  $_->{tStart},
	  $_->{tEnd},
	  $_->{qName},
	  $_->{qStart},
	  $_->{qEnd},
	  $_->{strand}
	]
    } @alndata_group;

    return \@segment_pairs;
}


sub _get_alndata_for_region
{
    my ($self, $query_id, $aln_type, $chr, $start, $end) = @_;
    my $table = 'WGA_'.$query_id.'_'.uc($aln_type);
    my $bin_string = $self->bin_restriction_string($start, $end);
    my $query = qq!
SELECT *
FROM $table
WHERE $bin_string
AND tName = ? AND tEnd >= ? AND tStart <= ?
ORDER BY tStart, qName, strand, qStart
!;
    my $sth = $self->dbh->prepare($query);
    $sth->execute($chr, $start, $end);
    my @alndata;
    while (my $row = $sth->fetchrow_hashref) {
	push @alndata, $self->_get_alndata_from_row($row);
    }
    $sth->finish;
    return \@alndata;
}


sub _get_alndata_for_regions
{
    my ($self, $query_id, $aln_type, $tChr, $tStart, $tEnd, $qChr, $qStart, $qEnd) = @_;
    my $table = 'WGA_'.$query_id.'_'.uc($aln_type);
    my $bin_string = $self->bin_restriction_string($tStart, $tEnd);
    my $query = qq!
SELECT *
FROM $table
WHERE $bin_string
AND tName = ? AND tEnd >= ? AND tStart <= ?
AND qName = ? AND qEnd >= ? AND qStart <= ?
ORDER BY tStart, qName, strand, qStart
!;
    my $sth = $self->dbh->prepare($query);
    $sth->execute($tChr, $tStart, $tEnd, $qChr, $qStart, $qEnd);
    my @alndata;
    while (my $row = $sth->fetchrow_hashref) {
	push @alndata, $self->_get_alndata_from_row($row);
    }
    $sth->finish;
    return \@alndata;
}


# before
# 1. SELECT tStart FROM WGA_HS_MAY04_NET WHERE tName = 'chr4' AND tStart < '43248561' order by tStart DESC limit 1;
# 2. _get_alndata_for_region(start, pos)

# after
# 1. SELECT tEnd FROM WGA_HS_MAY04_NET WHERE tName = 'chr4' AND tEnd > '43248561' order by tEnd limit 1;
# 2. _get_alndata_for_region(pos,end)

sub _get_first_alndata_before_pos
{
    my ($self, $query_id, $aln_type, $chr, $pos) = @_;
    my $table = 'WGA_'.$query_id.'_'.uc($aln_type);
    my ($query, $sth);

    $query = "SELECT tStart FROM $table WHERE tName = ? AND tStart < ? ORDER BY tStart DESC LIMIT 1";
    $sth = $self->dbh->prepare_cached($query);
    $sth->execute($chr, $pos);
    my ($start) = $sth->fetchrow_array();
    $sth->finish;
    return undef unless($start);

    my $alndata_set = $self->_get_alndata_for_region($query_id, $aln_type, $chr, $start, $pos);

    my $rightmost_alndata;
    my $max_tEnd = 0;
    foreach my $alndata (@$alndata_set) {
	my $tEnd =  $alndata->{tEnd};
	if($tEnd > $max_tEnd) {
	    $rightmost_alndata = $alndata;
	    $max_tEnd = $tEnd;
	}
    }
    
    return $rightmost_alndata;

#    $query = "SELECT * FROM $table WHERE tName = ? AND tEnd = ?";
#    $sth = $self->dbh->prepare_cached($query);
#    $sth->execute($chr, $end);
#    my $row = $sth->fetchrow_hashref();
#    my $alndata = $self->_get_alndata_from_row($row);
#    $sth->finish;
#    return $alndata;
}


sub _get_first_alndata_after_pos
{
    my ($self, $query_id, $aln_type, $chr, $pos) = @_;
    my $table = 'WGA_'.$query_id.'_'.uc($aln_type);
    my ($query, $sth);

    $query = "SELECT tEnd FROM $table WHERE tName = ? AND tEnd > ? ORDER BY tEnd LIMIT 1";
    $sth = $self->dbh->prepare_cached($query);
    $sth->execute($chr, $pos);
    my ($end) = $sth->fetchrow_array();
    $sth->finish;
    return undef unless($end);

    my $alndata_set = $self->_get_alndata_for_region($query_id, $aln_type, $chr, $pos, $end);
    return $alndata_set->[0];

#    $query = "SELECT * FROM $table WHERE tName = ? AND tStart = ?";
#    $sth = $self->dbh->prepare_cached($query);
#    $sth->execute($chr, $start);
#    my $row = $sth->fetchrow_hashref();
#    my $alndata = $self->_get_alndata_from_row($row);
#    $sth->finish;
#    return $alndata;
}


sub _get_alndata_from_row
{
    my ($self, $row) = @_;
    my %alndata = (%$row); # could use arrayref for row but less stable
    $alndata{tStarts} = [split(/,/,$alndata{tStarts})];
    $alndata{qStarts} = [split(/,/,$alndata{qStarts})];
    $alndata{blockSizes} = [split(/,/,$alndata{blockSizes})];
    return \%alndata;
}


sub _create_aln_from_alndata
{
    my ($self, $alndata, $tSeq_db, $qSeq_db, $tStart, $tEnd, $qStart, $qEnd) = @_;

    # If no data, no work to do
    return undef unless(@$alndata);

#print STDERR "creating alignment for ", scalar(@$alndata), " alndata\n";

    # tName and qName should be same for all subalignments
    my ($tName, $qName) = ($alndata->[0]->{tName},
			   $alndata->[0]->{qName});

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


    # Sometimes it would be more efficient only to retrieve sequences for the
    # subalignments; however Alignment::Composite stores the entire sequences
    # so we get them here for compatability.
    # If in the future we want to get subsequences only, this can be accomodated
    # by copying the _create_aln_from_seqs_and_alndata code into this method
    # and giving offsets for the subsequencesin arg 4 & 5 when calling _create_subaln

    # Create and return composite alignment
    return $self->_create_aln_from_seqs_and_alndata($tSeq, $qSeq, $alndata);        
}

   
sub _create_aln_from_seqs_and_alndata
{
    my ($self, $tSeq, $qSeq, $alndata) = @_;

    # If no data, no work to do
    return undef unless(@$alndata);

    # Get strand of query; should be same for all subalignments
    my $strand = $alndata->[0]->{strand};

    # Create a composite alignment
    my $alignment = AT::Alignment::Composite->new(-seq1 => $tSeq,
						  -seq2 => $qSeq,
						  -strand => $strand);

    # Create and add subalignments to composite alignment
    my $nr_subalignments = 0;    
    foreach my $a (@$alndata) {
	my $sa = $self->_create_subaln($tSeq, $qSeq, $a, 0, 0);
	if($sa) {
	    $alignment ->add_subalignment(-sub => $sa);
	    $nr_subalignments++;
	}
    }
 
    return $nr_subalignments ? $alignment : undef;
}


sub _create_subaln
{
    my ($self, $tSeq, $qSeq, $a, $tOffset_major, $qOffset_major) = @_;

    # Get sequence lengths
    my $tLength = $tSeq->length;
    my $qLength = $qSeq->length;

    # Get strand of query
    my $strand = $a->{strand};
    
    # Get some attributes into local variables
    my ($tStart, $tEnd) = ($a->{tStart}, $a->{tEnd});
    my ($qStart, $qEnd) = ($a->{qStart}, $a->{qEnd});
    my $tOffset = $tStart - $tSeq->start;
    my $qOffset = ($strand eq '+') ? ($qStart - $qSeq->start)
			           : ($qSeq->end - $qEnd);
    # Get coords for ungapped blocks
    my @blocks;
    my $tStarts = $a->{tStarts};
    my $qStarts = $a->{qStarts};
    my $blockSizes = $a->{blockSizes};
    foreach my $i (0..@$tStarts-1) {
	my ($s1, $s2, $e1, $e2) =
	    $self->_restrict_block
    		($tStarts->[$i] + $tOffset,
		 $qStarts->[$i] + $qOffset,
		 $tStarts->[$i] + $tOffset + $blockSizes->[$i]-1,
		 $qStarts->[$i] + $qOffset + $blockSizes->[$i]-1,
		 $tLength,
		 $qLength);
	    next if($s1 > $e1);
	    push @blocks, ['', $s1, $s2, $e1, $e2]; # leave first empty for compatability w Malin's subs
    }
    return unless(@blocks);

    # Create gapped sequences
    my ($alnseq1, $alnseq2) = $self->_create_gapped_seqs($tSeq, $qSeq, \@blocks, $strand);
    return unless($alnseq2);

    # Get start and end pos for alignment
    my @pos_seq1=($blocks[0][1] + $tOffset_major,
		  $blocks[-1][3] + $tOffset_major);
    my @pos_seq2=($blocks[0][2] + $qOffset_major,
		  $blocks[-1][4] + $qOffset_major);
    my @lengths = ($pos_seq1[1] - $pos_seq1[0] + 1,
    	           $pos_seq2[1] - $pos_seq2[0] + 1);

    # Create SimpleAlign object
    my $simple_align_obj = $self->_create_subalignment_from_gapped_seqs
	($alnseq1, $alnseq2,
	 $tSeq->id, $qSeq->id,
	\@pos_seq1, \@pos_seq2);

    # Create Subaligment object
    my $subalignment = AT::Alignment::Subalignment->new
	        (-ps1 => \@pos_seq1,
	         -ps2 => \@pos_seq2,
	         -score => $a->{score}, # note: this might not be the actual score we might not use the entire alignment
	         -ao => $simple_align_obj,
	         -strand => $strand,
	         -lengths => \@lengths,
	         -ub => \@blocks);

    return $subalignment;
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
  my ($self, $locatableseq1, $locatableseq2, $id1, $id2, $pos_seq1, $pos_seq2) = @_;

  #LocatableSeq objects are used to make an alignment:
  $locatableseq1 = Bio::LocatableSeq -> new(-seq => $locatableseq1, -id => $id1, -start =>$pos_seq1->[0], -end => $pos_seq1->[1]);
  $locatableseq2 = Bio::LocatableSeq -> new(-seq =>$locatableseq2, -id => $id2, -start=>$pos_seq2->[0], -end =>$pos_seq2->[1]);
  
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


