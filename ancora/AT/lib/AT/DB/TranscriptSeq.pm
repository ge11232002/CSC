# AT::DB::TranscriptSeq module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::DB::TranscriptSeq

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::DB::TranscriptSeq;

use strict;
use vars '@ISA';
use Carp;
use AT::DB::MySQLdb;
use AT::DB::Importer;
use AT::DB::DataSourceHolder;
use AT::Seq::Generic;

@ISA = qw(AT::DB::MySQLdb AT::DB::DataSourceHolder);

# The following three array-subroutine pairs define the column-method
# assocations between the SEQ table and Bio::Seq object
# Used by store_x and create_x_importer methods.

my @_SEQ_COLUMN_NAMES = qw(
    acc
    version
    division
    replaced
    locus
    type
    organism
    data_source_id
    description
    seq);

# Missing: replaced, type, source_id
# To provide for this extra information, we should subclass Bio::Seq
# (or use annotation)

sub _seq_column_values {
    my ($self, $seqobj) = @_;

    return ($seqobj->can('accession_number') ?
		$seqobj->accession_number : $seqobj->id,
	    $seqobj->version,
	    $seqobj->division,
	    $seqobj->replaced,
	    $seqobj->id,
	    $seqobj->type,
	    $seqobj->species->binomial,
	    $seqobj->data_source_id,
	    $seqobj->desc, 
	    $seqobj->seq);    
}


sub _init
{
    my ($self) = @_;
    my $data = $self->dbh->selectall_arrayref("SELECT id, value from CONFIG") or
	croak "Could not read CONFIG table";

    my $config = { map {@$_} @$data };
    $self->{'_config'} = $config;

    $config->{'schema_version'} ||= 100;

    if($config->{'schema_version'} >= 200) {
	my $org_ss_string = $config->{'organism_subsections'};
	unless($org_ss_string) {
	    warn "No organism subsections defined in CONFIG table for db ",$self->db_name,
	    ". Using deafult (human,mouse)\n";
	    $org_ss_string = 'human,mouse';
	}
	$self->{'_organism_subsection_hash'} =
	    { map { $_ => 1 } $self->_parse_comma_separated_config_string($org_ss_string) };
    }

    $self->{'lazy_load_seq'} = 0;
}


=head2 store_seq

 Title     : store_seq
 Usage     : $db->store_seq($seq);
 Function  : Stores a (mRNA or EST) sequence in the database
 Returns   : 1 on success
 Args      : A Bio::SeqI-compliant object

=cut


sub store_query_seq  {
    my ($self, $seqobj) = @_;

    croak "store_query_seq not updated to current schema";

    # Check argument
    croak "Wrong sequence object type" 
        unless (ref($seqobj) and $seqobj->isa("Bio::SeqI"));

    # Store sequence
    my $query =
        'INSERT INTO QUERY_SEQUENCES(' . (join ',', @_SEQ_COLUMN_NAMES) .
        ') VALUES(' . '?,'x(scalar(@_SEQ_COLUMN_NAMES)-1) . '?)';
    my $sth = $self->dbh->prepare($query);
    $sth->execute($self->_seq_column_values($seqobj));
    $sth->finish;
    return 1;
}


=head2 create_seq_importer

 Title     : create_seq_importer
 Usage     : my $importer = $db->create_seq_importer();
 Function  : Creates an AT::DB::Importer object that allows effective
             import of a large number of sequences into the database.
             See the AT::DB::Importer module for further details.
 Example   : my $importer = $db->create_seq_importer();
             while (my $seq = $seq_stream->next_seq) {
		 $importer->store_object($seq);
	     }
             my @result = $importer->finish;
 Returns   : An AT::DB::Importer object

=cut


sub create_seq_importer {
    my $self = shift;

    croak "create_seq_importer not updated to schema 2.00" if($self->{_config}{schema_version} >= 200);

    # define wrapper for storing one seq
    my $storing_wrapper = sub {
	my ($importer, $seqobj) = @_;
	unless (ref($seqobj) and $seqobj->isa("Bio::SeqI")) {
	    croak "Wrong sequence object type";
	}
	$importer->store_record('SEQ',
				[ $self->_seq_column_values($seqobj) ]);
   };

    # prepare table info needed by importer
    my $tables_columns = {
	'SEQ' => [ @_SEQ_COLUMN_NAMES ],
    };

    # create importer object
    my $importer = AT::DB::Importer->new( db => $self,
					  tables_columns => $tables_columns,
					  storing_wrapper => $storing_wrapper,
					  lock => 0); # no lock required

    # return the creation
    return $importer;
}


=head2 get_seq

 Title     : get_seq
 Usage     : my $seq = $db->get_seq("NM_001243", 1);
 Function  : Gets a sequence from the database by accession and
	     version number. If a version is not given, the most
	     recent sequence will be retrieved.
 Returns   : A Bio::Seq object with sequence and description on success.
             Undef if the accession number is not in the database.
 Args      : 1. Accession number
	     2. Version (optional)

=cut


sub get_seq {
    my $self = shift;

    if($self->{_config}{schema_version} < 200) {
    	return $self->_get_seq_v100(@_);
    }
    else {
	return $self->_get_seq_v200(@_);
    }
}


# add perldocs for these two

sub organism_subsection_list { keys %{shift->{_organism_subsection_hash}} }
sub organism_subsection_listref { [ keys %{shift->{_organism_subsection_hash}} ]; }

# add perldoc for lazy_load_seq() here : autoloaded method

sub _get_seq_v100 {
    my ($self, $acc, $req_version) = @_;

    my $version_string = $req_version ?
	" AND version = ?" : " AND replaced = 0";
    my $sth = $self->dbh->prepare_cached
	("SELECT DISTINCT ".(join ',', @_SEQ_COLUMN_NAMES)." ".
	 "FROM SEQ
	  WHERE acc = ?
	  $version_string");
    if($req_version) { $sth->execute($acc, $req_version); }
    else { $sth->execute($acc); }
    my $row = $sth->fetchrow_arrayref();
    $sth->finish;
    return undef unless($row);
    return $self->_create_seq_from_row_v100($row);
}


sub _get_seq_v200 {
    my ($self, $acc, $req_version, %options) = @_;

#options can be:
# seq_sets [mrna,est,refseq,...]
# organisms  [as configured (see CONFIG table): usually mouse, human or other]
# if no org groups given -> scan all org groups conf'd + other
# if some given -> translate those not conf'd to other 
# DO:
# Foreach org_group until seq found
#  If seq_types is subset of mrna,est,refseq)
#    If acc is refseq: query refseq if refseq in set
#    Else: query mrna if mrna in set; query est if est in set
#  Else    
#    Query all seq_types in turn until seq found

# Q: new sequence object needed to hold extra data?

# table name format:
# seq_section_organism
# (in mapping db we have map_mapset)
# 

    my $tables = $self->_get_seq_table_names($acc, $options{seq_sets}, $options{organisms});
    my $seq_column = $self->lazy_load_seq() ? ", seq" : "";
    my $version_string = $req_version ?	" AND version = ?"  : "";
    foreach my $table (@$tables) {
	my $sth = $self->dbh->prepare_cached
	("SELECT acc, version, mod_date, type,
                 organism, gene, product, length, cds,
                 library, mrna_clone, refseq_status,
                 read_dir, init_Ts, term_As,
                 description
                 $seq_column
	  FROM $table
	  WHERE acc = ?
	  $version_string
	  ORDER BY VERSION DESC");
	if($req_version) { $sth->execute($acc, $req_version); }
	else { $sth->execute($acc); }
	my $row = $sth->fetchrow_arrayref();
	$sth->finish;
	return $self->_create_seq_from_row_v200($row) if($row);  # sequence found
    }

    return undef;  # no sequence found
}


sub _get_seq_table_names
{
    my ($self, $acc, $seq_sets, $organisms) = @_;

    # Figure out which sequence sets to use
    my ($mrna_ss, $rs_ss, $est_ss, $custom_ss);
    foreach my $seq_set (@$seq_sets) {
	if($seq_set eq 'mrna') {
	    $mrna_ss = 1;
	}
	elsif($seq_set eq 'refseq') {
	    $rs_ss = 1;
	}
	elsif($seq_set eq 'est') {
	    $est_ss = 1;
	}
	else {
	    $custom_ss = 1;
	    last;
	}
    }
    unless($custom_ss) {
	if($acc =~ /^NM_/) {
	    $seq_sets = $rs_ss ? [ 'refseq' ] : [ ];
	}
	else {
	    my @seq_sets;
	    push @seq_sets, 'mrna' if($mrna_ss);
	    push @seq_sets, 'est' if($est_ss);
	    $seq_sets = \@seq_sets;
	}
    }

    # Figure out which organism subsections to use
    my $organism_subsections;
    if($organisms) {
	my $all_organism_subsections = $self->{_organism_subsection_hash};
	my %organism_subsection_hash;
	foreach my $o (@$organisms) {
	    if($all_organism_subsections->{$o}) {
		$organism_subsection_hash{$o} = 1;
	    }
	    else {
		$organism_subsection_hash{'other'} = 1;
	    }
	}
	$organism_subsections = [keys %organism_subsection_hash];
    }
    else {
	$organism_subsections = [keys %{$self->{_organism_subsection_hash}}];
    }
	
    # Produce table names
    my @table_names;
    foreach my $s (@$seq_sets) {
	foreach my $o (@$organism_subsections) {
	    push @table_names, "seq_${s}_${o}";
	}
    }

    # Return table names
    return \@table_names;
}


sub seq_table_name
{
    my ($self, $seq_set, $organism) = @_;
    $organism = 'other' unless($self->{_organism_subsection_hash}{$organism});
    my $table = "seq_${seq_set}_${organism}";
    return $table;
}


sub get_seqhash_using_acctable
{
# we use a foreign dbh, which is very ugly; TEMPORARY
    my ($self, $dbh, $acc_table, $acc_col, $ver_col) = @_;
    my $version_string = defined($ver_col) ?
	" AND SEQ.version = $acc_table.$ver_col" : " AND SEQ.replaced = 0";
    my $sth = $dbh->prepare_cached
	("SELECT DISTINCT ".(join ',', @_SEQ_COLUMN_NAMES)." ".
	 "FROM ".$self->dbname.".SEQ as SEQ, $acc_table
	  WHERE acc = $acc_table.$acc_col
	  $version_string");
    $sth->execute() or return undef;
    my %seqs;
    while(my $row = $sth->fetchrow_arrayref) {
	my $seq = $self->_create_seq_from_row_v100($row);
	$seqs{$seq->accession_number} = $seq;
    }
    return \%seqs;
}


sub _create_seq_from_row_v100 {
    my ($self, $row) = @_;
    my ($acc,
	$version,
	$division,
	$replaced,
	$locus,
	$type,
	$organism,
	$data_source_id,
	$desc,
	$seqstring) = @$row;
    my $seqobj = AT::Seq::Generic->new(-accession_number => $acc,
				       -version => $version,
				       -division => $division,
				       -replaced => $replaced,
				       -id => $locus,
				       -type => $type,
				       -desc=> $desc,
				       -seq=> $seqstring,
				       -alphabet => 'dna');
    if($organism ne 'unknown') {
	#my $speciesobj = Bio::Species->new(-classification =>
	#				   [ reverse split /\s+/, $organism ]); 
	#$seqobj->species($speciesobj);
	# Bio::Species does not work with current perl
    }
    return $seqobj;
}


sub _create_seq_from_row_v200 {
    my ($self, $row) = @_;

    my ($acc,
	$version,
	$mod_date,
	$type,
	$organism,
	$gene,
	$product,
	$length,
	$cds,
	$library,
	$clone,
	$refseq_status,
	$read_dir,
	$init_Ts,
	$term_As,
	$desc,
	$seqstring) = @$row;

    my $seqobj = AT::Seq::Transcript->new(-id => $acc,
					  -accession_number => $acc,
					  -version => $version,
					  -alphabet => 'dna',
					  -desc=> $desc,
					  -length => $length,
					  -seq=> $seqstring,
					  # non-standard attributes follow
					  -mod_date => $mod_date,
					  -type => $type,
					  -gene_name => $gene,
					  -product_name => $product,
					  -cds_string => $cds,
					  -library_id => $library,
					  -clone_id => $clone,
					  -refseq_status => $refseq_status,
					  -read_dir => $read_dir,
					  -init_Ts => $init_Ts,
					  -term_As => $term_As
					  );
    if($organism and $organism ne 'unknown') {
	#my $speciesobj = Bio::Species->new(-classification =>
	#				   [ reverse split /\s+/, $organism ]); 
	#$seqobj->species($speciesobj);
	# Bio::Species does not work with current perl
    }
    return $seqobj;
}


=head2 riken_clone_id_to_acc

 Title     : riken_clone_id_to_acc
 Usage     : my $acc = $db->riken_clone_id_to_acc("E130016A22");
 Function  : Gets the accession number corresponding to a riken clone id.
 Returns   : String
 Args      : String

=cut


sub riken_clone_id_to_acc
{
    my ($self, $clone_id) = @_;
    my ($acc) = $self->dbh->selectrow_array(
	"SELECT acc FROM RIKEN_ID WHERE clone_id = \"$clone_id\"");
    return $acc;
}


=head2 acc_to_riken_clone_id

 Title     : acc_to_riken_clone_id
 Usage     : my $acc = $db->acc_to_riken_clone_id("AK123456");
 Function  : Gets the riken clone id corresponding to an accession number.
 Returns   : String
 Args      : String

=cut


sub acc_to_riken_clone_id
{
    my ($self, $acc) = @_;
    my ($clone_id) = $self->dbh->selectrow_array(
	"SELECT clone_id FROM RIKEN_ID WHERE acc = \"$acc\"");
    return $clone_id;
}


=head2 nr_seqs

 Title     : nr_seqs
 Usage     : my $n = $db->nr_seqs();
 Returns   : The number of sequences in the database.

=cut


sub nr_seqs {
    my $self = shift;
    return $self->_num_records('SEQ');
}


=head2 is_wanted

 Title     : is_wanted
 Usage     : $self->is_wanted($acc, $version);
 Function  : Checks whether an identifier is on the wanted list.
	     Used by transcript_seqs2sql.pl
 Returns   : true or false
 Args      : accession, version (both required)

=cut

sub is_wanted {
    my ($self, $acc, $version) = @_;
    $self->_get_wanted_hash unless ($self->{_wanted_hash});
    return $self->{_wanted_hash}{$acc.'.'.$version};
}


=head2 _get_wanted_hash

 Title     : _get_wanted_hash
 Usage     : my $hash = $self->_get_wanted_hash();
 Function  : Sets $self->{_wanted_hash}.
 Returns   : -
 Args      : -

=cut

sub _get_wanted_hash {
    my ($self, $acc, $req_version) = @_;
    my $sth = $self->dbh->prepare
	("SELECT CONCAT(acc,'.',version) FROM WANTED");
    $sth->execute();
    my %hash;
    while(my $acc_ver = $sth->fetchrow_array()) {
	$hash{$acc_ver} = 1;
    }
    $sth->finish;
    $self->{_wanted_hash} = \%hash;
}


sub DESTROY  {
    my $self = shift;
    $self->dbh->disconnect() if $self->dbh();
}

1;


