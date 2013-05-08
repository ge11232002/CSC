# AT::DB::GenomeMapping module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::DB::GenomeMapping - interface to MySQL relational database of genome
mappings


=head1 SYNOPSIS

=over 4

=item * creating a database object by connecting to an existing AT-type database

my $db = AT::DB::GenomeMapping->connect(
    -dbname => "AT",
    -dbhost => "myhostn.mydomain",
    -dbuser => "myusername",
    -dbpass => "mypassword");

=item * retrieving previously stored genome mappings of the sequence AF100234

my @mappings = $db->get_mappings_for_acc("AF100234");

=item * storing the contents of an AT::Mapping object into the database

$db->store_mapping($mapping);

=back

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=head2 connect

 Title     : connect
 Usage     : my $db = $AT::DB::GenomeMapping->connect(%args);
 Function  : Creates a new object of this class.
 Note      : For details, see docs for superclass AT::DB::MySQLdb.
    
=cut


package AT::DB::GenomeMapping;

use strict;
use vars '@ISA';
use Bio::Seq;
use Bio::LocatableSeq;
use DBI;
use Carp;
use AT::DB::MySQLdb;
use AT::DB::Binner;
use AT::DB::Importer;
use AT::Mapping;
use AT::mRNAInfo;
use AT::HSP;
use AT::BlatpslxIO;


@ISA = qw(AT::DB::MySQLdb AT::DB::Binner);


my %SIMPLE_FEATURE_INFO =
    ('CpG island' => { table => 'cpgIslandExt',
		       is_binned => 0 },
     'gap' => {table => 'gap',
	       is_binned => 1 },
     'assembly' => {table => 'gold',
		    is_binned => 1 }
    );


# The following three array-subroutine pairs define the column-method
# assocations between three MySQL tables (MAPPING, HSP, QUERY_SEQUENCES)
# and their corresponding Perl objects.
# Used by store_x and create_x_importer methods.

my @_MAPPING_COLUMN_NAMES = qw(bin
			       mapping_id     target_db
			       matches        misMatches
			       repMatches     nCount
			       qNumInsert     qBaseInsert
			       tNumInsert     tBaseInsert
			       strand
			       qName          qSize
			       qStart         qEnd
			       qType
			       tName          tSize
			       tStart         tEnd);

sub _mapping_column_values {
    my ($self, $mapping) = @_;

    return ($mapping->bin,
	    $mapping->mapping_id,  $mapping->target_db,
	    $mapping->matches,     $mapping->misMatches,
	    $mapping->repMatches,  $mapping->nCount,
	    $mapping->qNumInsert,  $mapping->qBaseInsert,
	    $mapping->tNumInsert,  $mapping->tBaseInsert,
	    $mapping->strand,
	    $mapping->qName,       $mapping->qSize,
	    $mapping->qStart,      $mapping->qEnd,
	    $mapping->qType,
	    $mapping->tName,       $mapping->tSize,
	    $mapping->tStart,      $mapping->tEnd );
}

my @_HSP_COLUMN_NAMES =  qw(hsp_id       mapping_id
			    qStart       qEnd
			    tStart       tEnd
			    blockSize
			    qSeq         tSeq);

sub _hsp_column_values {
    my ($self, $hsp) = @_;

    return ($hsp->hsp_id,    $hsp->mapping_id,
	    $hsp->qStart,    $hsp->qEnd,
	    $hsp->tStart,    $hsp->tEnd,
	    $hsp->blockSize,
	    $hsp->qSeq,      $hsp->tSeq);
}

my @_QUERY_SEQ_COLUMN_NAMES = qw(acc version description seq);

sub _query_seq_column_values {
    my ($self, $seqobj) = @_;

    return ($seqobj->can('acc') ? $seqobj->acc : $seqobj->id,
	    $seqobj->version,
	    $seqobj->desc, 
	    $seqobj->seq);    
}


sub _init
{
    my ($self,%args) = @_;
    if($self->_table_exists("CONFIG")) {
	my $data = $self->dbh->selectall_arrayref("SELECT id, value from CONFIG") or
	    croak "Could not read CONFIG table";
	$self->{'_config'} = { map {@$_} @$data };
    }
    $self->{skip_if_no_mRNA_info} = $self->{_config}{use_mrnainfo} ? 1 : 0;
    $self->{query_seq_db} = $args{-query_seq_db} || $args{query_seq_db};
    $self->{_table_formats} = {};
}


=head2 store_mapping

 Title     : store_mapping
 Usage     : $db->store_mapping($mapping);
 Function  : Stores a mapping in the database
 Returns   : 1 on success
 Args      : An AT::Mapping object

=cut


sub store_mapping  {
    my ($self, $mapping) = @_;

    # Check argument
    unless ($mapping->isa("AT::Mapping")) { 
        croak "Wrong type of mapping object";
    }
    
    # Insert data into MAPPING table
    $mapping->mapping_id(0);   # Make sure auto_increment is used to set id
    my $mapping_query =
        'INSERT INTO MAPPING(' . (join ',', @_MAPPING_COLUMN_NAMES) .
        ') VALUES(' . '?,'x(scalar(@_MAPPING_COLUMN_NAMES)-1) . '?)';
    my $mapping_sth = $self->dbh->prepare($mapping_query);
    $mapping_sth->execute($self->_mapping_column_values($mapping))
        or do { warn "Error executing query.";
                return undef; };
    $mapping_sth->finish;
    my $mapping_id = $self->_last_insert_id;
    $mapping->mapping_id($mapping_id);

    # Insert data into HSP table
   my $hsp_query =
        'INSERT INTO HSP(' . (join ',', @_HSP_COLUMN_NAMES) .
        ') VALUES(' . '?,'x(scalar(@_HSP_COLUMN_NAMES)-1) . '?)';
    my $hsp_sth = $self->dbh->prepare($hsp_query);
    foreach my $hsp ($mapping->all_HSPs)  {
	$hsp->hsp_id(0);  # Make sure auto_increment is used to set id
	$hsp->mapping_id($mapping_id);
        $hsp_sth->execute($self->_hsp_column_values($hsp));
        $hsp->hsp_id($self->_last_insert_id);
    }
    $hsp_sth->finish;   
    return 1;
}


#=head2 store_query_seq
#
# Title     : store_query_seq
# Usage     : $db->store_query_seq($seq);
# Function  : Stores a query (mRNA or EST) sequence in the database
# Returns   : 1 on success
# Args      : A Bio::SeqI-compliant object
#
#=cut
#
#
#sub store_query_seq  {
#    my ($self, $seqobj) = @_;
#
#    # Check argument
#    croak "Wrong sequence object type" 
#        unless (ref($seqobj) and $seqobj->isa("Bio::SeqI"));
#
#    # Store sequence
#    my $query =
#        'INSERT INTO QUERY_SEQUENCES(' . (join ',', @_QUERY_SEQ_COLUMN_NAMES) .
#        ') VALUES(' . '?,'x(scalar(@_QUERY_SEQ_COLUMN_NAMES)-1) . '?)';
#    my $sth = $self->dbh->prepare($query);
#    $sth->execute($self->_query_seq_column_values($seqobj));
#    $sth->finish;
#    return 1;
#}


=head2 create_mapping_importer

 Title     : create_mapping_importer
 Usage     : my $importer = $db->create_mapping_importer();
 Function  : Creates an AT::DB::Importer object that allows effective
             import of a large number of mappings to the database.
             See the AT::DB::Importer module for further details.
 Example   : my $importer = $db->create_mapping_importer();
             while (my $mapping = $mapping_stream->next_mapping) {
		 $importer->store_object($mapping);
	     }
             my @result = $importer->finish;
 Returns   : An AT::DB::Importer object

=cut


sub create_mapping_importer {
    my $self = shift;

    # define wrapper for storing one mapping
    my $storing_wrapper = sub {
	my ($importer, $mapping) = @_;
	unless ($mapping->isa("AT::Mapping")) { 
	    croak "Wrong type of mapping object";
	}
	unless (defined $mapping->qType) {
	    croak "Type undefined for mapping object".
		($mapping->mapping_id || '');
	}
	my $mapping_id = $importer->counter;
	$importer->counter($importer->counter+1);
	$mapping->mapping_id($mapping_id);
	$mapping->bin($self->bin_from_coord_range($mapping->tStart,
						  $mapping->tEnd));
	$importer->store_record('MAPPING',
				[$self->_mapping_column_values($mapping)]);
	foreach my $hsp ($mapping->all_HSPs)  {
	    $hsp->hsp_id(0);  # Make sure auto_increment is used to set id
	    $hsp->mapping_id($mapping_id);
	    $importer->store_record('HSP',
				    [$self->_hsp_column_values($hsp)]);
	    # Note: $hsp->hsp_id is not set.
	}
   };

    # prepare table info needed by importer
    my $tables_columns = {
	MAPPING => [ @_MAPPING_COLUMN_NAMES ],
	HSP => [ @_HSP_COLUMN_NAMES ]
    };

    # create importer object; tell it to lock database
    my $importer = AT::DB::Importer->new( db => $self,
					  tables_columns => $tables_columns,
					  storing_wrapper => $storing_wrapper,
					  lock => 1);

    # set mapping_id counter
    #   (mapping_ids have to be set "manually" since we can't use
    #    last_insert_id in batch storage)
    if ($self->num_mappings) {
	my $id_sth = $self->dbh->prepare('SELECT MAX(mapping_id) FROM MAPPING');
	$id_sth->execute or do { warn "Error executing query.";
			     return undef; };
	my ($highest_id) = $id_sth->fetchrow_array;
	$id_sth->finish;
	$importer->counter($highest_id+1);
    }
    else {
	$importer->counter(1);
    }

    # return the creation
    return $importer;
}


#=head2 create_query_seq_importer
#
# Title     : create_query_seq_importer
# Usage     : my $importer = $db->create_query_seq_importer();
# Function  : Creates an AT::DB::Importer object that allows effective
#             import of a large number of query sequences to the database.
#             See the AT::DB::Importer module for further details.
# Example   : my $importer = $db->create_query_seq_importer();
#             while (my $seq = $seq_stream->next_seq) {
#		 $importer->store_object($seq);
#	     }
#             my @result = $importer->finish;
# Returns   : An AT::DB::Importer object
#
#=cut
#
#
#sub create_query_seq_importer {
#    my $self = shift;
#
#    # define wrapper for storing one seq
#    my $storing_wrapper = sub {
#	my ($importer, $seqobj) = @_;
#	unless (ref($seqobj) and $seqobj->isa("Bio::SeqI")) {
#	    croak "Wrong sequence object type";
#	}
#	$importer->store_record('QUERY_SEQUENCES',
#				[ $self->_query_seq_column_values($seqobj) ]);
#   };
#
#    # prepare table info needed by importer
#    my $tables_columns = {
#	QUERY_SEQUENCES => [ @_QUERY_SEQ_COLUMN_NAMES ],
#    };
#
#    # create importer object
#    my $importer = AT::DB::Importer->new( db => $self,
#					  tables_columns => $tables_columns,
#					  storing_wrapper => $storing_wrapper,
#					  lock => 0); # no lock required
#
#    # return the creation
#    return $importer;
#}


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


=head2 get_features_by_id

 Title     : get_features_by_id
 Usage     : my @ft = $db->get_features_by_id(id => "NM_001243",
					      types => ['mRNA',
					                'EST']
					     );
 Function  : Gets genomic features by id from the database.
 Returns   : An array of AT::Mapping objects.
 Args      : id - typically an accession number
             types - types of features to look for; array reference.
 Notes:
 Types 'mRNA' and 'EST' are accepted for all databases.
 Some databases have features of other types, e.g.:
 'genscan', 'ensGene', 'flyBaseGene', 'knownGene'.
 The names of these other feature types should correspond to table names
 in the database (these table names are usually the same as in the
 original UCSC database).
 For mRNA/EST mappings, this method simply calls get_mappings_by_acc,
 and the id argument maps to the qName column in the db table.
 For other features, the id argument maps to the name column in the db table.
 Currently all features are returned as AT::Mapping objects. This will
 change, as e.g. gene annotations are not mappings and may be better
 represented by a different class (with a compatible interface).

=cut


sub get_features_by_id {
    my ($self, %args) = @_;

    my $types = $args{types};
    croak "Missing or invalid types argument" unless($types and ref($types) eq 'ARRAY');
    my $id = $args{id} || croak "No id argument";

    my @features;

    my $get_mappings = 0;
  
    foreach my $type (@$types) {
	my $lc_type = lc $type;
	if($lc_type eq 'mrna' or $lc_type eq 'est') {
	    $get_mappings = 1;
        }
	else {
	    my $genes = $self->_get_genes_by_id($id, $type);
	    my $mappings = $self->_genes_to_mappings($genes);
	    push @features, @$mappings;
	}
    }
    push @features, $self->get_mappings_for_acc($id) if($get_mappings);

    return @features;
}


=head2 get_features_in_region

 Title     : get_features_in_region
 Usage     : my @ft = $db->get_features_in_region(chr => "chr10",
						  start => "2134000",
						  end => "2140000",
						  types => ['mRNA',
						            'ensGene']
						);
 Function  : Gets all the features that, completely or partially, are
             in a particlar genomic region.
 Returns   : An array of AT::Mapping objects.
 Args      : chr, start, end - specifies the region; all required
             strand - set to '+' or '-' select features from
               one strand only; optional
             types - types of features to look for; array reference;
	       required
 Notes:
 See get_features_by_id() above.
 For types 'mRNA' and 'EST', this method simply calls get_mappings_in_region.
    
=cut


sub get_features_in_region
{
    my ($self, %args) = @_;
    my $types = delete $args{types};
    croak "Missing or invalid types argument" unless($types and ref($types) eq 'ARRAY');
    my $chr = $args{'chr'} or croak "No chr argument";
    my $start = $args{'start'} or croak "No start argument";
    my $end = $args{'end'} or croak "No end argument";
    my $strand = $args{'strand'};

    my @features;
    my @mapping_types;

    foreach my $type (@$types) {
	my $lc_type = lc $type;
	if($lc_type eq 'mrna' or $lc_type eq 'est') {
	    push @mapping_types, $type;
        }
	else {
	    my $genes = $self->_get_genes_in_region($chr, $start, $end, $strand, $type);
	    my $mappings = $self->_genes_to_mappings($genes);
	    push @features, @$mappings;
	}
    }
    push @features, $self->get_mappings_in_region(%args, types => \@mapping_types)
	if(@mapping_types);

    return @features;
}


=head2 get_mappings_for_acc

 Title     : get_mappings_for_acc
 Usage     : my @mappings = $db->get_mappings_for_acc("NM_001243");
 Function  : Gets all the mappings for a particular query sequence from 
             the database.
 Returns   : An array of AT::Mapping objects.
 Args      : An accession number for a query sequence.

=cut


sub get_mappings_for_acc {
    my ($self, $acc, $target_db, %args) = @_;

    # Only get mappings to a certain assembly? 
    my $target_db_condition = defined $target_db ? 
	"and target_db = ".$self->dbh->quote($target_db) : "";

    # Get mRNA info
    my $mRNAInfo = $self->get_mRNA_info($acc) || 0;

    # Get mappings
    my $sth = $self->dbh->prepare(qq!SELECT DISTINCT *
				  FROM MAPPING WHERE MAPPING.qName = ?
				  $target_db_condition
				  ORDER BY matches desc!);
    $sth->execute($acc);
    my @mappings;
    while (my $hashref = $sth->fetchrow_hashref)  {
        push @mappings, $self->_create_mapping($hashref, undef, $mRNAInfo);
    }

    # Attach query sequences to mappings if asked to
    if($args{attach_query_seq}) {
	warn 'argument "attach_query_seq" currently unsupported';
	#my $query_seq = $self->get_query_seq($acc) ||
	 #   croak "Could not get query sequence $acc";
	#foreach my $mapping (@mappings) {
	 #   $mapping->attach_query_seq($query_seq);
	#}
    }

    return @mappings;
}


sub get_first_mapping_starting_in_region
{
# experimental method
    my ($self, %args) = @_;
    my ($chr) = ($args{'chr'} or $args{'-chr'} or
		 croak "No chr argument");
    my ($start) = ($args{'start'} or $args{'-start'} or 
		   croak "No start argument");
    my ($end) = ($args{'end'} or $args{'-end'} or $start);
    my ($strand) = ($args{'strand'} or $args{'-strand'});
    my ($target_db) = ($args{'target_db'} or $args{'-target_db'});
    my ($types) = ($args{'types'} or $args{'-types'});

    # Determine options to select statement
    my $map_table = 'MAPPING';
    my $where_options = "";
    # Position conditions (WHERE clause)
    if($self->{'_config'}{'use_binning'}) {
	$where_options .= ' and '.$self->bin_restriction_string($start, $end);
    }
    if(defined($target_db)) {
	$where_options .= " and $map_table.target_db = ".
	    $self->dbh->quote($target_db);
    }
    if(defined($strand)) {
	if($strand eq '-1') {$strand = '-';}
	elsif($strand eq '1') {$strand = '+';}
	$where_options .= " and $map_table.strand = ".$self->dbh->quote($strand);
    }
    if(defined($types) and @$types) {
	my $tmp_str = join " or $map_table.qType = ",
	    map {$self->dbh->quote($_)} @$types;
	$where_options .= " and ($map_table.qType = $tmp_str)";
    }

    my $query = qq!
	SELECT min(tStart)
	FROM $map_table
	WHERE tName = ?
	AND tStart >= ?
	AND tStart <= ?
	$where_options
    !;
    my $sth = $self->dbh->prepare($query);
    $sth->execute($chr, $start, $end);
    $start = $sth->fetchrow_array;
    $sth->finish();
    return unless ($start);
    return $self->get_mappings_in_region(%args, start => $start, end => $start, last => 1);
}


=head2 get_mappings_in_region

 Title     : get_mappings_in_region
 Usage     : my @mappings = $db->get_mappings_in_region(chr => "chr10",
						        start => "2134000",
						        end => "2140000");
 Function  : Gets all the mappings to a particlar genomic region.
 Returns   : An array of AT::Mapping objects.
 Args      : chr, start, end - specifies the region; all required
             Optional args:
             strand - set to '+','-',1 or -1 to select mappings from
               one strand only
             confine - set to 1 to only retrieve mappings entirely
               within the region
             confine_start - set to 1 to only retrieve mappings that
               start within the region
             confine_end - set to 1 to only retrieve mappings that
               end within the region
             first - get only the mapping with the lowest start value
             types - reference to array of strings. Currently 'mRNA'
               and 'EST' are the only types. E.g. types => ['mRNA']
               only retrieves mRNA mappings. By default mappings of
	       all types are retrieved.
    
=cut


sub get_mappings_in_region  {
    my ($self, %args) = @_;
    my ($chr) = ($args{'chr'} or $args{'-chr'} or
		 croak "No chr argument");
    my ($start) = ($args{'start'} or $args{'-start'} or 
		   croak "No start argument");
    my ($end) = ($args{'end'} or $args{'-end'} or $start);
    my ($strand) = ($args{'strand'} or $args{'-strand'});
    my ($target_db) = ($args{'target_db'} or $args{'-target_db'});
    my ($types) = ($args{'types'} or $args{'-types'});

    # Determine options to select statement
    my $map_table = 'MAPPING';
    my $from = $map_table;
    my $where_options = "";
    my $order_and_limit_options = "";
    # Position conditions (WHERE clause)
    if($self->{'_config'}{'use_binning'} and
       !defined($args{'no_bin'}) and !defined($args{'-no_bin'})) {
	$where_options .= ' and '.$self->bin_restriction_string($start, $end);
    }
    $where_options .= ($args{'confine'} || $args{'-confine'} ||
		       $args{'confine_start'} || $args{'-confine_start'}) ?
		       " and $map_table.tStart >= ? " : " and $map_table.tEnd >= ? ";
    $where_options .= ($args{'confine'} || $args{'-confine'} ||
		       $args{'confine_end'} || $args{'-confine_end'}) ?
		       " and $map_table.tEnd <= ? " : " and $map_table.tStart <= ? ";
    # Optional conditions (WHERE clause)
    if(defined($target_db)) {
	$where_options .= " and $map_table.target_db = ".
	    $self->dbh->quote($target_db);
    }
    if(defined($strand)) {
	if($strand eq '-1') {$strand = '-';}
	elsif($strand eq '1') {$strand = '+';}
	$where_options .= " and $map_table.strand = ".$self->dbh->quote($strand);
    }
    if(defined($types) and @$types) {
	my $tmp_str = join " or $map_table.qType = ",
	    map {$self->dbh->quote($_)} @$types;
	$where_options .= " and ($map_table.qType = $tmp_str)";
    }
    # Ordering, limits
    #  (we need to join in the mrna table for the first and last arguments
    #   to this function to work properly - otw in case we get a mapping without
    #   an mrna entry first, this mapping will be tossed later, and will not
    #   return any mapping; this is getting messy; maybe a separate method is
    #   needed for first and last)
    if($args{'first'} or $args{'-first'}) {
	$from .= ', mrna';
	$where_options .= " and $map_table.qName = mrna.acc ";
	$order_and_limit_options .= " order by $map_table.tStart limit 1";
    }
    elsif($args{'last'} or $args{'-last'}) {
	$from .= ', mrna';
	$where_options .= " and $map_table.qName = mrna.acc ";
	$order_and_limit_options .= " order by $map_table.tEnd desc limit 1";
    }
    else {
	$order_and_limit_options .= " order by $map_table.matches*".
	    "($map_table.matches+$map_table.misMatches)/$map_table.qSize desc";
    }

    # Create tmp table with mapping data
#print STDERR "making tmp table\n";
    my $query1 = qq!
		CREATE TEMPORARY TABLE
		_AT_tmp._mapping
		TYPE=HEAP
		SELECT
		$map_table.*
		FROM $from
		WHERE $map_table.tName = ?
		$where_options
		$order_and_limit_options
		!;
    my $sth1 = $self->dbh->prepare($query1);
    $sth1->execute($chr, $start, $end);

    # Get HSP data and create HSP objects
#print STDERR "getting HSPs\n";
    my $query3 = qq!
		SELECT HSP.* from
		_AT_tmp._mapping as tmp, HSP
		where tmp.mapping_id = HSP.mapping_id
		!;
    my $sth3 = $self->dbh->prepare_cached($query3);
    $sth3->execute;
    my %HSPs;
    while(my $hashref = $sth3->fetchrow_hashref) {
	push @{$HSPs{$hashref->{mapping_id}}}, $self->_create_HSP($hashref);
    }

    # Get mrna data and create mrnaInfo objects
#print STDERR "getting mrna info\n";
    my $mRNAInfo = $self->_get_mRNA_info_batch('_AT_tmp._mapping')
	if($self->{'_config'}{'use_mrnainfo'});

    # Get mapping data and create mapping objects
#print STDERR "creating mappings\n";
    my $sth2 = $self->dbh->prepare_cached("SELECT * FROM _AT_tmp._mapping");
    $sth2->execute();
    my @mappings;
    while (my $hashref = $sth2->fetchrow_hashref)  {
	my $mapping_id = $hashref->{mapping_id};
        push @mappings, $self->_create_mapping($hashref,
					       $HSPs{$mapping_id},
					       $mRNAInfo->{$mapping_id} || 0);
    }
#print STDERR "done\n";

    # Attach query sequences to mappings if asked to
    if(my $seqdb = $args{attach_query_seq}) {
	my $seqhash = $seqdb->get_seqhash_using_acctable($self->dbh,'_AT_tmp._mapping',
							 'qName'); #ugly
	foreach my $mapping (@mappings) {
	    my $seq = $seqhash->{$mapping->qName};
	    $mapping->attach_query_seq($seq) if($seq);
	}
    }

    # Drop temporary table
    $self->dbh->do("DROP TABLE _AT_tmp._mapping");

    return @mappings;
}


#sub get_mappings_in_region_old  {
#    my ($self, %args) = @_;
#    my ($chr) = ($args{'chr'} or $args{'-chr'} or
#		 croak "No chr argument");
#    my ($start) = ($args{'start'} or $args{'-start'} or 
#		   croak "No start argument");
#    my ($end) = ($args{'end'} or $args{'-end'} or $start);
#    my ($strand) = ($args{'strand'} or $args{'-strand'});
#    my ($target_db) = ($args{'target_db'} or $args{'-target_db'});
#    my ($types) = ($args{'types'} or $args{'-types'});
#
#    # Determine options to select statement
#    my $map_table = 'MAPPING';
#    my $where_options = "";
#    my $order_and_limit_options = "";
#    my $use_query_seq_table = 0;
#    # Position conditions (WHERE clause)
#    if($self->{'_config'}->{'use_binning'} and
#       !defined($args{'no_bin'}) and !defined($args{'-no_bin'})) {
#	$where_options .= ' and '.$self->bin_restriction_string($start, $end);
#    }
#    $where_options .= ($args{'confine'} || $args{'-confine'} ||
#		       $args{'confine_start'} || $args{'-confine_start'}) ?
#		       " and $map_table.tStart >= ? " : " and $map_table.tEnd >= ? ";
#    $where_options .= ($args{'confine'} || $args{'-confine'} ||
#		       $args{'confine_end'} || $args{'-confine_end'}) ?
#		       " and $map_table.tEnd <= ? " : " and $map_table.tStart <= ? ";
#    # Optional conditions (WHERE clause)
#    if(defined($target_db)) {
#	$where_options .= " and $map_table.target_db = ".
#	    $self->dbh->quote($target_db);
#    }
#    if(defined($strand)) {
#	if($strand eq '-1') {$strand = '-';}
#	elsif($strand eq '1') {$strand = '+';}
#	$where_options .= " and $map_table.strand = ".$self->dbh->quote($strand);
#    }
#    if(defined($types) and @$types) {
#	#my $tmp_str = join " or QUERY_SEQUENCES.type = ",
#	#    map {$self->dbh->quote($_)} @$types;
#	#$where_options .= " and (QUERY_SEQUENCES.type = $tmp_str)";
#	#$use_query_seq_table = 1;
#	my $tmp_str = join " or $map_table.qType = ",
#	    map {$self->dbh->quote($_)} @$types;
#	$where_options .= " and ($map_table.qType = $tmp_str)";
#    }
#    # Ordering, limits
#    if($args{'first'} or $args{'-first'}) {
#	$order_and_limit_options .= " order by $map_table.tStart limit 1";
#    }
#    elsif($args{'last'} or $args{'-last'}) {
#	$order_and_limit_options .= " order by $map_table.tEnd desc limit 1";
#    }
#    else {
#	$order_and_limit_options .= " order by $map_table.matches*".
#	    "($map_table.matches+$map_table.misMatches)/$map_table.qSize desc";
#    }
#    my $tables = $map_table;
#    if($use_query_seq_table) {
#	$tables .= ", QUERY_SEQUENCES";
#	$where_options .= " and $map_table.qName = QUERY_SEQUENCES.acc";
#    }
#
#    # Get mappings
##    my $query = qq!
##		SELECT
##		$map_table.*, mrna.library, mrna.mrnaClone
##		FROM $tables, mrna
##		WHERE $map_table.tName = ? AND MAPPING.qName = mrna.acc
##		$where_options
##		$order_and_limit_options
##		!;
#    my $query = qq!
#		SELECT
#		$map_table.*
#		FROM $tables
#		WHERE $map_table.tName = ?
#		$where_options
#		$order_and_limit_options
#		!;
#    #my $print_q = $query; $print_q =~ s/\?/\"$chr\"/; $print_q =~ s/\?/\"$start\"/;
#    #$print_q =~s/\?/\"$end\"/; print STDERR "$print_q\n";
#    my $sth = $self->dbh->prepare($query);
#    $sth->execute($chr, $start, $end);
#    #print STDERR "query done\n";
#    my @mappings;
#    while (my $hashref = $sth->fetchrow_hashref)  {
#        push @mappings, $self->_create_mapping($hashref);
#    }
#
#    # Attach query sequences to mappings if asked to
#    if($args{attach_query_seq}) {
#	warn 'argument "attach_query_seq" currently unsupported';
#	#foreach my $mapping (@mappings) {
#	#    my $query_seq = $self->get_query_seq($mapping->qName) ||
#	#	croak("Could not get query sequence ".$mapping->qName);
#	#    $mapping->attach_query_seq($query_seq);
#	#}
#    }
#
#    return @mappings;
#}




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
    my $sth = $self->dbh->prepare(qq!SELECT *
				  FROM MAPPING WHERE mapping_id = ?!);
    $sth->execute($id);
    my $hashref = $sth->fetchrow_hashref;
    my $mapping = AT::Mapping->new
	(db => $self,
	 %$hashref,
	 HSPs=>[$self->_get_HSPs_for_mapping_id($hashref->{'mapping_id'})]
	 );
    return $mapping;
}


=head2 get_distance_to_HSP

 Title     : get_distance_to_HSP
 Usage     : my $d = $db->get_distance_to_HSP
		(chr => "chr22",
	         pos => 17007486,
		 max_dist => -5000,
		 types => ['mRNA']);							  
 Function  : Get distance from a given position to the nearest HSP
 Returns   : Undef if there are no HSPs within the given region.
	     0 if some HSP overlaps pos.
	     Otherwise positive or negative integer.
 Args      : chr - chromosome
	     pos - position on chromosome
	       Typically this is set to the start of some mapping - 1.
	       Note: if pos is set to the exact start of some mapping,
	       the result will always be 0, since the fist HSP of that
	       mapping overlaps pos.
	     max_dist - the maximum distance to consider
	       If >0, the region [pos, pos+max_dist] will be searched.
	       If <0, the region [pos+max_dist, pos] will be searched.
	       Typically this is set to -5000.
	     Optional args:
             types - reference to array of strings. Currently 'mRNA'
               and 'EST' are the only types. E.g. types => ['mRNA']
               to only consider mRNA mappings. By default mappings of
	       all types are considered.

=cut

sub get_distance_to_HSP {
    my ($self, %args) = @_;

    my $chr = $args{'chr'} || croak "No chr arg";
    my $pos = $args{'pos'} || croak "No start arg";
    my $max_d = $args{'max_dist'} || croak "No max_dist arg";
    my $types = $args{'types'};

	my ($start, $end, $border, $order);
	if($max_d > 0) {
		($start, $end) = ($pos, $pos + $max_d);
		$border = 'tStart';
		$order = 'ASC';
	}
	else {
		($start, $end) = ($pos + $max_d, $pos);
		$border = 'tEnd';
		$order = 'DESC';
	}
    
    my $bin_arg = "";
    if($self->{'_config'}->{'use_binning'} and
       !defined($args{'no_bin'}) and !defined($args{'-no_bin'})) {
		$bin_arg = $self->bin_restriction_string($start, $end). " AND ";
    }
    my $type_arg = "";
    if(defined($types) and @$types) {
		my $tmp_str = join " or MAPPING.type = ",
		    map {$self->dbh->quote($_)} @$types;
		$type_arg = " (MAPPING.qType = $tmp_str) AND ";
    }

    my $sth = $self->dbh->prepare(
		  "SELECT HSP.$border FROM MAPPING, HSP
		   WHERE MAPPING.mapping_id = HSP.mapping_id AND
	       MAPPING.tName = ? AND
	       $bin_arg
	       $type_arg
	       MAPPING.tStart <= ? AND MAPPING.tEnd >= ? AND
	       HSP.tStart <= ? AND HSP.tEnd >= ?
		   ORDER BY $border $order");
    unless(defined $sth) { croak "Error querying database"; }

	$sth->execute($chr, $end, $start, $end, $start);
	my ($nearest) = $sth->fetchrow_array();
	$sth->finish;

	# If there is no HSP in the range, return undef
    return undef unless ($nearest);

	# Calculate relative distance. Set to 0 if the nearest HSP overlaps $pos
	my $d = $nearest - $pos;
	$d = 0 if(($max_d > 0 and $d < 0) or ($max_d < 0 and $d > 0));

	# Return the distance
    return $d; 
}


#sub OLD_get_distance_to_HSP {
#    my ($self, %args) = @_;
#
#    my $chr = ($args{'chr'} or $args{'chr'} or croak "No chr arg");
#    my $pos = ($args{'pos'} or $args{'pos'} or croak "No start arg");
#    my $dir = $args{'dir'} || croak "No dir arg";
#    my $max_d = $args{'max_dist'};
#    my ($types) = ($args{'types'} or $args{'-types'});
#    
#    my $bin_arg = "";
#    if($self->{'_config'}->{'use_binning'} and
#       !defined($args{'no_bin'}) and !defined($args{'-no_bin'})) {
#		$bin_arg = $self->bin_restriction_string($pos, $pos). " AND ";
#    }
#    my $type_arg = "";
#    if(defined($types) and @$types) {
#		my $tmp_str = join " or MAPPING.type = ",
#		    map {$self->dbh->quote($_)} @$types;
#		$type_arg = " (MAPPING.qType = $tmp_str) AND ";
#    }
#
#    # First check if any HSPs overlap with the position
#    my ($nr_overlapping) = $self->dbh->selectrow_array(
#		"SELECT count(*) FROM MAPPING, HSP
#		 WHERE MAPPING.mapping_id = HSP.mapping_id AND
#	       MAPPING.tName = '$chr' AND
#	       $bin_arg
#	       $type_arg
#	       MAPPING.tStart <= $pos AND
#	       MAPPING.tEnd >= $pos AND
#	       HSP.tStart <= $pos AND
#	       HSP.tEnd >= $pos");
#    unless(defined $nr_overlapping) { croak "Error querying database"; }
#    if($nr_overlapping) { return 0; }  # we have overlapping HSPs -> dist is 0
#
#    # If we get here, there are no overlapping HSPs
#    # Find distance to the closest
#    # Note: we do a straight join to make sure that the index HSP.tStart
#    #  or HSP.tEnd) is used
#    my $query;
#    if($dir == 1) {
#	my $max_d_arg = $max_d ?
#	    "HSP.tStart <= ".($pos+$max_d)." AND " : "";
#	$query = "SELECT HSP.tStart FROM HSP STRAIGHT_JOIN MAPPING
#		    WHERE MAPPING.mapping_id = HSP.mapping_id AND
#			  MAPPING.tName = '$chr' AND
#			  $max_d_arg
#			  $type_arg
#			  HSP.tStart > $pos
#		    ORDER BY HSP.tStart
#		    LIMIT 1";
#    }
#    elsif($dir == -1) {
#	my $max_d_arg = $max_d ?
#	    "HSP.tEnd >= ".($pos-$max_d)." AND " : "";
#	$query = "SELECT HSP.tEnd FROM HSP STRAIGHT_JOIN MAPPING
#		    WHERE MAPPING.mapping_id = HSP.mapping_id AND
#			  MAPPING.tName = '$chr' AND
#			  $max_d_arg
#			  $type_arg
#			  HSP.tEnd < $pos	  
#		    ORDER BY HSP.tEnd DESC
#		    LIMIT 1";
#    }
#    else { croak "Illegal dir argument value $dir"; }
#    my ($min_d) = $self->dbh->selectrow_array($query);
#
#    return undef unless ($min_d); # should check for error
#    return abs($pos - $min_d); 
#}


sub get_nearest_HSPs {
    my ($self, %args) = @_;

    # call get_distance_to_HSP to get distance
    # get all HSP objects with that distance

}


sub get_refSeqStatus {
    my ($self, $id) = @_;
    print STDERR "GenomeMapping::get_refSeqStatus is deprecated\n";
    my ($status) = $self->dbh->selectrow_array
	("select status from refSeqStatus where mrnaAcc = \"$id\"");
    return $status;
}


sub get_library_and_clone_ids {
    my ($self, $id) = @_;
    print STDERR "GenomeMapping::get_library_and_clone_ids is deprecated\n";
    my ($library, $clone) = $self->dbh->selectrow_array
	("select library, mrnaClone from mrna where acc = \"$id\"");
    return ($library, $clone);
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
	    cds.name as cds
	FROM mrna, cds
	WHERE acc = ? AND mrna.cds = cds.id
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


sub get_libs_with_accurate_ori_annotation
{
    my ($self, %args) = (@_);
    unless($self->_table_exists('LIB_ORIANN_ACCURACY')) {
	warn "No annotation accuracy assessment for ".($self->dbname);
	return [];
    }
    my $min_nr_ests_annori = $args{'min_nr_ests_for_estimate'} or croak "Missing argument";
    my $min_ppt_correct = $args{'min_ppt_correct'} or croak "Missing argument";
    my $query = "SELECT libid FROM LIB_ORIANN_ACCURACY
		 WHERE nr_ests_annori >= $min_nr_ests_annori
		 AND correct_ppt >= $min_ppt_correct";
    my $id_list = $self->dbh->selectcol_arrayref($query);
    return [@$id_list];
}


sub _get_mRNA_info_batch {
    my ($self, $mapping_table) = @_;
    my %info;
    my $query1 = qq!
	SELECT
	    mapping_id as _mapping_id,
	    mrna.acc,
	    mrna.version,
	    mrna.type,
	    mrna.library,
	    mrna.mrnaClone,
	    cds.name as cds,
	    refSeqStatus.status as refSeqStatus
	FROM $mapping_table, mrna, cds
	LEFT JOIN refSeqStatus ON $mapping_table.qName = refSeqStatus.mrnaAcc
	WHERE $mapping_table.qName = mrna.acc AND mrna.cds = cds.id
        !;
    my $sth1 = $self->dbh->prepare_cached($query1);     
    $sth1->execute();
    while(my $hashref = $sth1->fetchrow_hashref) {
	$info{$hashref->{_mapping_id}} = AT::mRNAInfo->new(%$hashref);
    }
    return \%info;
}


=head2 num_mappings

 Title     : num_mappings
 Usage     : my $n = $db->num_mappings();
 Returns   : The number of mappings in the database.

=cut


sub num_mappings {
    my $self = shift;
    return $self->_num_records("MAPPING");
}


=head2 num_query_seqs

 Title     : num_mappings
 Usage     : my $n = $db->num_mappings();
 Returns   : The number of query sequences in the database.

=cut


sub num_query_seqs {
    my $self = shift;
    return $self->_num_records("QUERY_SEQUENCES");
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


sub _create_HSP {
    my ($self, $data) = @_;
    if($data->{blockSize} == 0) {
        warn ("blockSize is 0 for hsp_id ".$data->{hsp_id});
        return;
    }
    return AT::HSP->new(%$data);
}


=head2 _get_HSPs_for_mapping_id

 Title     : _get_HSPs_for_mapping_id
 Usage     : my @mappings = $db->_get_HSPs_for_mapping_id(1145);
 Function  : Gets all the HSPs for a particlar mapping
 Returns   : An array of AT::HSP objects.
 Args      : A mapping id.

=cut


sub _get_HSPs_for_mapping_id  {
    my ($self, $mapping_id) = @_;
    
    my $sth = $self->dbh->prepare(q!SELECT DISTINCT *
				  FROM HSP WHERE mapping_id =? 
				  ORDER BY qStart!);
    $sth->execute($mapping_id);
    my @HSPs;    
    while (my $hashref = $sth->fetchrow_hashref)  {
        push @HSPs, $self->_create_HSP($hashref);
    }
    return @HSPs;
}



sub _get_genes_by_id
{
    my ($self, $id, $table) = @_;

    my $sth = $self->dbh->prepare("select * from $table where name = ?");
    $sth->execute($id);
    my $genes = $self->_get_genes_from_query($sth);
    return $genes;
}


sub _get_genes_in_region
{
    my ($self, $chr, $start, $end, $strand, $table) = @_;

    my $table_format = $self->_get_table_format($table);
    my $query;
    if($table_format eq 'genePred') {
        $query = "SELECT * FROM $table WHERE chrom = ? AND txEnd >= ? AND txStart < ?";
    }
    elsif($table_format eq 'bed') {
        $query = "SELECT * FROM $table WHERE chrom = ? AND chromEnd >= ? AND chromStart < ?";
    }
    else {
        croak "Unsupported gene table format $table_format";
    }
    if(defined $strand) {
	$query .= " AND strand = '$strand'";
    }
    my $sth = $self->dbh->prepare($query);
    $sth->execute($chr, $start, $end);
    my $genes = $self->_get_genes_from_query($sth);
    return $genes;
}


sub _get_genes_from_query
{
    my ($self, $sth) = @_;
    my @genes;
    while(my $genedata_ref = $sth->fetchrow_hashref) {
        my %gene = %$genedata_ref; # create new hash since DBI may reuse references
	# Split comma-separated fields to arrays
	foreach my $field qw (exonStarts exonEnds blockSizes chromStarts) {
	    $gene{$field} = [ split ',',$gene{$field} ] if(defined $gene{$field});
	}
	# Convert from zero-based, half-open to one-based coords
    	foreach my $field qw(txStart cdsStart chromStart thickStart) {
	    $gene{$field}++ if(defined $gene{$field});
        }
	if(my $starts = $gene{exonStarts}) {
	    for my $i (0..@$starts-1) {
		$starts->[$i]++;
	    }
	}
	push @genes, \%gene;
    }
    return \@genes;
}


sub _genes_to_mappings
{
# note: this method is destructive: data is removed from the gene hashes that are input
    my ($self, $genes) = @_;
    my @mappings;
    foreach my $gene (@$genes) {
	my $strand = delete $gene->{strand};
	my $qName = delete $gene->{name};
	my $tName = delete $gene->{chrom};
	my $tStart = delete $gene->{txStart} || delete $gene->{chromStart};
	my $tEnd = delete $gene->{txEnd} || delete $gene->{chromEnd};
	my $blockCount = delete $gene->{exonCount} || delete $gene->{blockCount};
	my @hsps;
	my $starts;
	if($starts = delete $gene->{exonStarts}) {
	    my $ends = delete $gene->{exonEnds};
	    for my $i (0..$blockCount-1) {
		my ($s, $e) = ($starts->[$i], $ends->[$i]);
		push @hsps, AT::HSP->new(tStart => $s, tEnd => $e, blockSize => $e-$s+1);
	    }
	}
	elsif($starts = delete $gene->{chromStarts}) {
	    my $sizes = delete $gene->{blockSizes};
	    for my $i (0..$blockCount-1) {
		my $start = $starts->[$i] + $tStart;
		my $size = $sizes->[$i];
		push @hsps, AT::HSP->new(tStart => $start,
					 tEnd => $start+$size-1,
					 blockSize => $size);
	    }
        }
	else { croak "Missing exon starts for gene" }
	my $mapping = AT::Mapping->new(
		       db => $self,
		       strand => $strand,
		       qName => $qName,
		       tName => $tName,
		       tStart => $tStart,
		       tEnd => $tEnd,
		       HSPs => \@hsps,
		       %$gene);
	push @mappings, $mapping;

    }
    return \@mappings;
}


sub _get_table_format
{
# "smart" subroutine that guesses the format of a table
# currently only works for genepred and bed formats
    my ($self, $table) = @_;
    my $formats = $self->{_table_formats};
    my $format = $formats->{$table};
    unless($format) {
	my @columns = $self->_get_column_names($table);
	foreach my $c (@columns) {
	    if($c eq 'txStart') {
		$format = "genePred";
		last;
	    }
	    elsif($c eq 'chromStart') {
		$format = "bed";
		last;
	    }
	}
	croak "unrecognized table or table format for table $table" unless($format);
	$formats->{table} = $format;
    }
    return $format;    
}


sub get_simple_features_in_region
{
    my ($self, %args) = @_;
    my $sth = $self->_execute_simple_feature_query("*", %args);
    my @all_ft_data;
    while(my $data_ref = $sth->fetchrow_hashref) {
	my %feature_data = %$data_ref;
	$feature_data{chromStart}++;
	push @all_ft_data, \%feature_data;
    }
    return @all_ft_data;
}


sub get_simple_feature_locations_in_region
{
    my ($self, %args) = @_;
    my $sth = $self->_execute_simple_feature_query("chromStart, chromEnd", %args);
    my @locations;
    while(my ($start,$end) = $sth->fetchrow_array) {
	push @locations, [$start+1, $end];
    }
    return @locations;
}


sub _execute_simple_feature_query
{
    my ($self, $select_fields, %args) = @_;
    my $chr = $args{'chr'} or croak "No chr argument";
    my $start = $args{'start'} or croak "No start argument";
    my $end = $args{'end'} or croak "No end argument";
    my ($table, $is_binned);
    if($args{'feature_type'}) {
	my $feature = $args{'feature_type'};
	my $feature_info = $SIMPLE_FEATURE_INFO{$feature} or croak "Unknown feature type $feature";
	$table = $feature_info->{table};
	$is_binned = $feature_info->{is_binned};
    }
    elsif($args{'table'}) {
	$table = $args{'table'};
	$is_binned = $args{'table_is_binned'};
    }

    my $query = "SELECT $select_fields
	         FROM $table
		 WHERE chrom = ? AND chromEnd >= ? AND chromStart <= ?";
    if($is_binned) {
	my $bin_string = $self->bin_restriction_string($start, $end);
	$query .= " AND ($bin_string)";
    }
    my $sth = $self->dbh->prepare($query); 
    croak "Failed to prepare query [$query]" unless($sth);
    $sth->execute($chr, $start, $end);

    return $sth;
}


sub DESTROY  {
    my $self = shift;
    $self->dbh->disconnect() if $self->dbh();
}

1;


