# AT::DB::GenomeMapping2 module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::DB::GenomeMapping2 - interface to MySQL relational database of genome
mappings


=head1 DESCRIPTION

This is a temporary, experimental module to test a new revision of the
mapping db schema. This schema is used for the FANTOM3 mapping db (F3_AT_MM_MAY04).
This module will eventually be removed, as AT::DB::GenomeMapping is upgreaded
to support a newer, more efficient schema.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

THE MODULE DOCUMENTATION IS INCOMPLETE. SEE SOURCE CODE.

=head2 connect

 Title     : connect
 Usage     : my $db = $AT::DB::GenomeAssembly->connect(%args);
 Function  : Creates a new object of this class.
 Note      : For details, see docs for superclass AT::DB::MySQLdb.

=cut

package AT::DB::GenomeMapping2;

use strict;
use vars '@ISA';
use Bio::Seq;
use Bio::LocatableSeq;
use DBI;
use Carp;
use AT::DB::MySQLdb;
use AT::DB::Binner;
use AT::DB::Importer;
use AT::DB::ResultStream;
use AT::Mapping;
use AT::mRNAInfo;
use AT::HSP;
use AT::BlatpslxIO;


@ISA = qw(AT::DB::MySQLdb AT::DB::Binner);


sub _init
{
    my ($self, %args) = @_;
    my $data = $self->dbh->selectall_arrayref("SELECT id, value from CONFIG") or
	croak "Could not read CONFIG table";
    $self->{'_config'} = { map {@$_} @$data };
    $self->{skip_if_no_mRNA_info} = $self->{_config}{use_mrnainfo} ? 1 : 0;
    $self->{query_seq_db} = $args{-query_seq_db} || $args{query_seq_db};
}


=head2 get_query_seq

 Title     : get_query_seq
 Usage     : my $seq = $db->get_query_seq("NM_001243");
 Function  : Gets query sequence using associated query seq db.
             The call is simply passed on to the query seq db:
             $self->query_seq_db->get_seq(@_);
 Returns   : A Bio::Seq object with sequence and description on success.
             Undef if no query seq db is associated or no sequence was
             found.
 Args      : 1. Accession number
             2. Version number (optional)

=cut


sub get_query_seq {
    my $self = shift;
    if(my $query_seq_db = $self->query_seq_db) {
	return $query_seq_db->get_seq(@_);
    }
    else {
	return undef;
    }
}


sub get_libs_with_accurate_ori_annotation
{
    my ($self, %args) = (@_);
    unless($self->_table_exists('LIB_ORIANN_ACCURACY')) {
	warn "No annotation accuracy assessment for ".($self->dbname);
	return [];
    }
#    my $min_nr_ests_annori = $args{'min_nr_ests_for_estimate'} or croak "Missing argument";
    my $confidence_level = $args{'confidence_level'} or croak "Missing argument confidence_level";
    my $min_ppt_correct = $args{'min_ppt_correct'} or croak "Missing argument min_ppt_correct";
    my $query = "SELECT libid FROM LIB_ORIANN_ACCURACY
		 WHERE correct_ppt_CI${confidence_level}lo >= $min_ppt_correct";
    my $id_list = $self->dbh->selectcol_arrayref($query);
    return [@$id_list];
}



=head2 get_mapping_by_id

 Title     : get_mapping_by_id
 Usage     : my $mapping = $db->get_mapping_by_id(112);
 Function  : Gets the mapping with a particular mapping id from the
             database.
 Returns   : An AT::Mapping object.
 Args      : A mapping id.
 Note      : This may not be a good way to get mappings, as mapping ids are
             primarily indented as internal db keys and will change if the
             database is rebuilt. 

=cut


sub get_mapping_by_id {
    my ($self, $id) = @_;
    return $self->_get_mappings("a.mapping_id = ?", [$id])->[0];
}


sub get_mappings_by_cluster_id
{
    my ($self, $id, %args) = @_;
    return $self->_get_mappings("a.cluster_id = ?", [$id], %args);
}


#sub get_mRNA_mappings_by_cluster_id
#{
#    my ($self, $id, %args) = @_;
#    return $self->_get_mappings("a.cluster_id = ? AND qType = \"mRNA\"", [$id], %args);
#}


sub get_mappings_for_acc {
    my ($self, $acc, %args) = @_;
    my $mappings = $self->_get_mappings("a.qName = ?", [$acc], %args);
    return $mappings ? @$mappings : ();
}


sub get_mappings_in_region
{
    my ($self, %args) = @_;

    my $chr = $args{chr} || die "no chr arg";
    my $start = $args{start} || die "no start arg";
    my $end = $args{end} || die "no end arg";

    my $bin_string = $self->bin_restriction_string($start,$end);
    my $where_string = " (tName = '$chr' AND $bin_string AND tEnd >= $start and tStart <= $end) ";

    my $mappings_ref = $self->_get_mappings($where_string,
					    [],
					    dont_cache_query => 1,
					    %args);

    return @$mappings_ref;
}


sub get_mappings_in_region_as_cluster_stream
{
    my ($self, $chr, $start, $end, %args) = @_;

    my $cluster_query = "SELECT cluster_id FROM MAP_CLUSTER
			 WHERE tName = ? AND tEnd >= ? AND tStart <= ?
			 ORDER BY tName, tStart";
    my $cluster_sth = $self->dbh->prepare($cluster_query);
    $cluster_sth->execute($chr, $start, $end);

    my $next_method = sub {
	my ($cluster_id) = $cluster_sth->fetchrow_array();
	return unless($cluster_id);
	$self->get_mappings_by_cluster_id($cluster_id, %args);
    };

    return AT::DB::ResultStream->new(next_method => $next_method);
}


sub get_mRNA_info {
    my ($self, $id) = @_;
    my $query1 = qq!
	SELECT
	    acc,
	    version,
	    type,
	    library,
	    mrnaClone,
	    cds
	FROM mrna
	WHERE acc = ?
        !;
    my $sth1 = $self->dbh->prepare_cached($query1);     
    $sth1->execute($id);
    my $data1 = $sth1->fetchrow_hashref();
    $sth1->finish;
    my $refSeqStatus;
    if($id =~ /^NM_/) {
	my $query2 = "select status from refSeqStatus where mrnaAcc = ?";
	my $sth2 = $self->dbh->prepare_cached($query2);
	$sth2->execute($id);
	($refSeqStatus) = $sth2->fetchrow_array();
	$sth2->finish;
    }
    my $info = AT::mRNAInfo->new(%$data1,
				 refSeqStatus => $refSeqStatus);
    return $info;


}



sub _get_mappings {
    my ($self, $condition, $bindvars, %args) = @_;

    my ($query, $sth);
    my $mapping_table = "MAPPING2";

    if($args{table}) {
	$mapping_table = $args{table};
    }

    my $strand_req = "";
    if($args{strand}) {
	my $strand = $args{strand};
	if($strand ne '+' and $strand ne '-') { die "Strand argument is $strand; must be + or -"; }
	$strand_req = " AND strand = '$strand' ";
    }

    my $type_req = "";
    if($args{types} and $args{types} eq 'mRNA') { $type_req = "AND qType = 'mRNA'" }
    my $types = $args{types};
    if(defined($types) and @$types) {
	my $tmp_str = join " OR a.qType = ",
	    map {$self->dbh->quote($_)} @$types;
	$type_req = " AND (a.qType = $tmp_str)";
    }

    my ($acc_table, $acc_req) = ("","");
    if($args{acc_table}) {
	$acc_table = $args{acc_table};
	$acc_req = " AND a.qName = $acc_table.acc";
	$acc_table = ",$acc_table";
    }
    my ($id_table, $id_req) = ("","");
    if($args{id_table}) {
	$id_table = $args{id_table};
	$acc_req = " AND a.mapping_id = $id_table.mapping_id";
	$acc_table = ",$id_table";
    }

    my $status_req;
    my $status_option = $args{status_requirement} || '';
    if($status_option eq 'none') {
	$status_req = "";
    }
    elsif($status_option eq 'one_best') {
	$status_req = "AND (MAP_STATUS.status = 'P' OR (MAP_STATUS.status = 'A' AND MAP_STATUS.note < 2))";
    }
    elsif($status_option eq '' or $status_option eq 'max_two_best') {
	$status_req = "AND (MAP_STATUS.status = 'P' OR (MAP_STATUS.status = 'A' AND MAP_STATUS.note < 3))";
    }
    else {
	croak "Unknown status_requirement $status_option";
    }

    $query =
    "SELECT a.*,
	    MAP_STATUS.status,
	    MAP_STATUS.note,
	    mrna.version,
	    mrna.library,
	    mrna.mrnaClone,
	    mrna.cds,
	    refSeqStatus.status
	FROM $mapping_table a, MAP_STATUS, mrna
             $acc_table
             $id_table
	LEFT JOIN refSeqStatus ON a.qName = refSeqStatus.mrnaAcc
	WHERE a.mapping_id = MAP_STATUS.mapping_id
	  AND a.qName = mrna.acc
          $strand_req
          $type_req
	  $status_req
          $acc_req
          $id_req
	  AND ($condition)";

#    $query =
#    "SELECT a.*,
#	    MAP_STATUS.status,
#	    MAP_STATUS.note,
#	    mrna.version,
#	    mrna.library,
#	    mrna.mrnaClone,
#	    cds.name,
#	    refSeqStatus.status
#	FROM $mapping_table a, MAP_STATUS, mrna, cds
#	LEFT JOIN refSeqStatus ON a.qName = refSeqStatus.mrnaAcc
#	WHERE a.mapping_id = MAP_STATUS.mapping_id
#	  AND a.qName = mrna.acc AND mrna.cds = cds.id
#	  $status_req
#	  AND ($condition)";

    $sth = $args{dont_cache_query} ? $self->dbh->prepare($query) : $self->dbh->prepare_cached($query);
    $sth->execute(@$bindvars);
    my @mappings;
    while(my $row = $sth->fetchrow_arrayref()) {
        push @mappings, $self->_mapping_from_row($row);
    }
    $sth->finish();

    return \@mappings;
}


sub _mapping_from_row
{
    my ($self, $row) = @_;

   my ($mapping_id,
       $cluster_id,
       $qType,    
       $bin,
       $matches,
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
	$changed,
	$qReadDir,
	$percentId,
	$score,
	$qInitTs,
	$qTermAs,
	$status,
	$status_note,
	$version,
	$library,
	$mrnaClone,
	$cds,
	$refSeqStatus
	) = @$row;
   
    $qStart++;
    $tStart++;

    my @HSPs;
    $blockSizes =~ s/,\s*$//;
    $qStarts    =~ s/,\s*$//;
    $tStarts    =~ s/,\s*$//;
    my @blockSizeList = split /,\s*/, $blockSizes;
    my @qStartList    = split /,\s*/, $qStarts;
    my @tStartList    = split /,\s*/, $tStarts;
    for my $i (0..@blockSizeList-1)  {
	push @HSPs, AT::HSP->new
	    (qStart    => $qStartList[$i] + 1,
	     qEnd      => $qStartList[$i] + $blockSizeList[$i],
	     tStart    => $tStartList[$i] + 1,
	     tEnd      => $tStartList[$i] + $blockSizeList[$i],
	     blockSize => $blockSizeList[$i],
	     );
    }

#    if($library == 0) {
#	if(length($qName) == 10) {
#	    $library = "ri".substr($qName,0,5);
#	}
#	elsif(length($qName) > 10) {
#	    $library = "ri.".substr($qName,3,5);
#	}
#   }

    my $mRNAInfo = AT::mRNAInfo->new(_mapping_id => $mapping_id,
				     acc => $qName,
				     version => $version,
				     type => $qType,
				     library => $library,
				     mrnaClone => $mrnaClone,
				     cds => $cds || '',
				     refSeqStatus => $refSeqStatus);

   return AT::Mapping->new(mapping_id	 => $mapping_id,
			   cluster_id	 => $cluster_id,
			   qType	 => $qType,
			   bin		 => $bin,
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
			    changed	 => $changed,
			    qReadDir	 => $qReadDir,
			    percentId	 => $percentId,
			    score	 => $score,
			    qInitTs	 => $qInitTs,
			    qTermAs	 => $qTermAs,
			    status	 => $status,
			    status_note  => $status_note,
			    HSPs         => [@HSPs],
			    db		 => $self,
			    mRNAInfo	 => $mRNAInfo);
}


=head2 _create_mapping

 Title     : _create_mapping
 Usage     : my $mappings = $db->_create_mapping($mapping_data);
 Function  : Creates an AT::Mapping object from a database record
 Returns   : An AT::Mapping object
 Args      : A reference to a hash containing the data from a row
             in the MAPPING table.

=cut


sub _create_mapping {
    my ($self, $data, $hsps, $mRNAInfo) = @_;

    if(!$mRNAInfo and $self->skip_if_no_mRNA_info) {
	warn("Skipping mapping ".$data->{qName}." due to missing mRNA info\n");
	return;
    }

    return AT::Mapping->new
	(%$data,
	 db => $self,
	 HSPs => $hsps || [$self->_get_HSPs_for_mapping_id($data->{'mapping_id'})],
	 mRNAInfo => $mRNAInfo
    );
}


sub DESTROY  {
    my $self = shift;
    $self->dbh->disconnect() if $self->dbh();
}

1;


