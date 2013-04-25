# CNE::DB
#
# Copyright 2007 Lenhard Group, BCCS, University of Bergen
# You may distribute this module under the same terms as perl itself
#

=head1 NAME

CNE::DB - interface to MySQL relational database of CNEs

=head1 SYNOPSIS

=over 4

=item * creating a database object by connecting to an existing CNE-type database

my $db = CNE::DB->connect(
    dbname => "cne",
    dbhost => "myhostn.mydomain",
    dbuser => "myusername",
    dbpass => "mypassword");

=item * retrieving CNEs for region...


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


package CNE::DB;

use strict;
use vars '@ISA';
use DBI;
use Carp;
use AT::DB::MySQLdb;
use AT::DB::Binner;
use AT::DB::Importer;
use CNE::CNE;


@ISA = qw(AT::DB::MySQLdb AT::DB::Binner);


sub _init
{
    my ($self,%args) = @_;
    # put any initialization here
}


=head2 get_cne_table_name

 Title     : get_cne_table_name
 Usage     : my $table = $db->get_cne_table_name(assembly1 => $asm1,
						 assembly2 => $asm2,
						 min_length => $length,
						 min_identity => $min_id,
						 version => $version);
 Returns   : Name of CNE table. The table need not exist.
 Arge      : ...

=cut

sub get_cne_table_name {
    my ($self, %args) = @_;
    my $asm1 = $args{assembly1} or croak "Missing argument: assembly1";
    my $asm2 = $args{assembly2} or croak "Missing argument: assembly2";
    my $len = $args{min_length} or croak "Missing argument: min_length";
    my $id = $args{min_identity} or croak "Missing argument: min_identity";
    my $version = $args{version} || 1;
    return $self->_get_cne_table_name($asm1, $asm2, $len, $id, $version);
}


sub get_all_cne_table_names {
    my ($self) = @_;
    my @cne_tables = grep { $_ =~ /^cne_/ } $self->_get_table_names();
    return \@cne_tables;
}


sub get_cne_table_names_for_assembly {
# this method has not been tested
    my ($self, $asm) = @_;
    my $all_tables = $self->get_all_cne_table_names();
    my @some_tables;
    foreach my $t (@$all_tables) {
	my ($asm1, $asm2) = $self->_get_assemblies_from_table_name($t);
	push @some_tables, $t if($asm1 eq $asm1 or $asm2 eq $asm);
    }
    return \@some_tables;
}


sub get_cne_table_names_for_assemblies {
    my ($self, $asm1, $asm2) = @_;
    ($asm1, $asm2) = ($asm2, $asm1) if($asm1 gt $asm2);
    my @cne_tables = grep { $_ =~ /^cne_twoWay_${asm1}_${asm2}_/ } $self->_get_table_names();
    return \@cne_tables;
}


sub get_cne_table_info {
    my ($self, $table) = @_;
    my ($type, $asm1, $asm2, $len, $id, $version) = $table =~ /^cne_([^_]+)_([^_]+)_([^_]+)_len(\d+)_id(\d+)_v(\d+)$/;
    croak "Illegal table name $table" unless(defined $version);
    return { assembly1 => $asm1,
	     assembly2 => $asm2,
	     min_length => $len,
	     min_identity => $id / 1000,
	     version => $version };
}


sub _get_cne_table_name {
    my ($self, $asm1, $asm2, $len, $id, $version) = @_;
    ($asm1, $asm2) = ($asm2, $asm1) if($asm1 gt $asm2);
    $id = int($id * 1000);
    return "cne_twoWay_${asm1}_${asm2}_len${len}_id${id}_v${version}";
}


sub _get_assemblies_from_table_name {
    my ($self, $table) = @_;
    my ($asm1, $asm2) = $table =~ /^cne_twoWay_([^_]+)_([^_]+)_/;
    croak "Illegal table name $table" unless($asm1 and $asm2);
    return ($asm1, $asm2);
}


=head2 get_cnes_in_region

 Title     : get_cnes_in_region
 Usage     : my $cnes = $db->get_cnes_in_region(table_name => $table_name,
						assembly => "hg18",
						chr => "chr10",
						start => 2134000,
						end => 2140000,
						);
 Function  : Gets all the CNEs that, completely or partially, are
             in a particlar genomic region.
 Returns   : Reference to an array of CNE::CNE objects.
 Args      : table_name - name of CNE table, obtained through a call
               to get_cne_table_name()
             assembly - name of the assembly coordinates refer to
             chr, start, end - specifies the region
             All arguments are required.

=cut

sub get_cnes_in_region
{
    my $self = shift;

    # Parse args
    my ($table, $assembly, $chr, $start, $end, $min_length) = $self->_parse_region_args(@_);

    # Build query
    my $query;
    my ($asm1, $asm2) = $self->_get_assemblies_from_table_name($table);
    if($assembly eq $asm1) {
	my $bin_string = $self->bin_restriction_string($start, $end, 'bin1');
	$query = "SELECT id, chr1, start1, end1, chr2, start2, end2, strand, similarity
	          FROM $table
		  WHERE chr1 = ? AND end1 >= ? AND start1 <= ? AND ($bin_string)";
    }
    elsif($assembly eq $asm2) {
	($asm1, $asm2) = ($asm2, $asm1);
	my $bin_string = $self->bin_restriction_string($start, $end, 'bin2');
	$query = "SELECT id, chr2, start2, end2, chr1, start1, end1, strand, similarity
	          FROM $table
		  WHERE chr2 = ? AND end2 >= ? AND start2 <= ? AND ($bin_string)";
    }
    else {
	croak "Table $table does not have CNEs for assembly $assembly";
    }
    $query .= " AND end1-start1 >= $min_length AND end2-start2 >= $min_length" if($min_length);

    # Prepare and execute query
    my $sth = $self->dbh->prepare($query) or croak "Failed to prepare query [$query]";
    $sth->execute($chr, $start-1, $end) or croak "Failed to execute query [$query]";

    # Fetch data and create CNE objects
    my $cnes = $self->_get_cnes_from_sth($asm1, $asm2, $table, $sth);

    return $cnes;
}


sub get_cne_ranges_in_region
{
    my $self = shift;

    # Parse args
    my ($table, $assembly, $chr, $start, $end, $min_length) = $self->_parse_region_args(@_);

    # Build query
    my $query;
    my ($asm1, $asm2) = $self->_get_assemblies_from_table_name($table);
    if($assembly eq $asm1) {
	my $bin_string = $self->bin_restriction_string($start, $end, 'bin1');
	$query = "SELECT start1, end1
	          FROM $table
		  WHERE chr1 = ? AND end1 >= ? AND start1 <= ? AND ($bin_string)";
    }
    elsif($assembly eq $asm2) {
	my $bin_string = $self->bin_restriction_string($start, $end, 'bin2');
	$query = "SELECT start2, end2
	          FROM $table
		  WHERE chr2 = ? AND end2 >= ? AND start2 <= ? AND ($bin_string)";
    }
    else {
	croak "Table $table does not have CNEs for assembly $assembly";
    }
    $query .= " AND end1-start1 >= $min_length AND end2-start2 >= $min_length" if($min_length);

    # Prepare and execute query
    my $sth = $self->dbh->prepare($query) or croak "Failed to prepare query [$query]";
    $sth->execute($chr, $start-1, $end) or croak "Failed to execute query [$query]";

    # Fetch data and create CNE objects
    my @cnes;
    while(my ($start1, $end1) = $sth->fetchrow_array) {
	push @cnes, [$start1+1, $end1];
	next;
    }

    return \@cnes;
}


sub get_cne_ranges_in_region_partitioned_by_other_chr
{
    my $self = shift;

    # Parse args
    my ($table, $assembly, $chr, $start, $end, $min_length) = $self->_parse_region_args(@_);

    # Build query
    my $query;
    my ($asm1, $asm2) = $self->_get_assemblies_from_table_name($table);
    if($assembly eq $asm1) {
	my $bin_string = $self->bin_restriction_string($start, $end, 'bin1');
	$query = "SELECT start1, end1, chr2
	          FROM $table
		  WHERE chr1 = ? AND end1 >= ? AND start1 <= ? AND ($bin_string)";
    }
    elsif($assembly eq $asm2) {
	my $bin_string = $self->bin_restriction_string($start, $end, 'bin2');
	$query = "SELECT start2, end2, chr1
	          FROM $table
		  WHERE chr2 = ? AND end2 >= ? AND start2 <= ? AND ($bin_string)";
    }
    else {
	croak "Table $table does not have CNEs for assembly $assembly";
    }
    $query .= " AND end1-start1 >= $min_length AND end2-start2 >= $min_length" if($min_length);

    # Prepare and execute query
    my $sth = $self->dbh->prepare($query) or croak "Failed to prepare query [$query]";
    $sth->execute($chr, $start-1, $end) or croak "Failed to execute query [$query]";

    # Fetch data and create CNE objects
    my %cnes;
    while(my ($start1, $end1, $chr2) = $sth->fetchrow_array) {
	push @{$cnes{$chr2}}, [$start1+1, $end1];
	next;
    }

    return \%cnes;
}


sub _parse_region_args {
    my ($self, %args) = @_;
    my $table = $args{table_name} or croak "No table name argument";
    my $assembly = $args{assembly} or croak "No assembly argument";
    my $chr = $args{chr} or croak "No chr argument";
    my $start = $args{start} or croak "No start argument";
    my $end = $args{end} or croak "No end argument";
    my $min_length = $args{min_length};
    return ($table, $assembly, $chr, $start, $end, $min_length);
}


=head2 get_cne_by_id

 Title     : get_cne_by_id
 Usage     : my $cne = $db->get_cne_by_id(table_name => $table_name,
					  assembly => "hg18",
					  id => 110);
 Function  : Gets the CNE with the specified id
 Returns   : A CNE::CNE object
 Args      : table_name - name of CNE table, obtained through a call
               to get_cne_table_name(); required
             assembly - name of the assembly; optional
             id - cne id; required

=cut

sub get_cne_by_id
{
    my ($self, %args) = @_;

    # Parse args
    my $table = $args{table_name} or croak "No table name argument";
    my $assembly = $args{assembly};
    my $id = $args{id} or croak "No id argument";

    # Build query
    my ($asm1, $asm2) = $self->_get_assemblies_from_table_name($table);
    my $columns;
    if(!defined($assembly) or $assembly eq $asm1) {
	$columns = "id, chr1, start1, end1, chr2, start2, end2, strand, similarity";
    }
    elsif($assembly eq $asm2) {
	($asm1, $asm2) = ($asm2, $asm1);
	$columns = "id, chr2, start2, end2, chr1, start1, end1, strand, similarity";
    }
    else {
	croak "Table $table does not have CNEs for assembly $assembly";
    }
    my $query = "SELECT $columns FROM $table WHERE id = ?";

    # Prepare and execute query
    my $sth = $self->dbh->prepare($query) or croak "Failed to prepare query [$query]";
    $sth->execute($id) or croak "Failed to execute query [$query]";

    # Fetch data and create CNE object
    my $cnes = $self->_get_cnes_from_sth($asm1, $asm2, $table, $sth);

    return $cnes->[0];
}


=head2 get_all_cnes_in_table

 Title     : get_all_cnes_in_table
 Usage     : my $cne = $db->get_all_cnes_in_table(table_name => $table_name,
					          assembly => "hg18");
 Function  : Gets all CNEs in specified table
 Returns   : A CNE::CNE object
 Args      : table_name - name of CNE table, obtained through a call
               to get_cne_table_name(); required
             assembly - name of the assembly; optional

=cut

sub get_all_cnes_in_table
{
    my ($self, %args) = @_;

    # Parse args
    my $table = $args{table_name} or croak "No table name argument";
    my $assembly = $args{assembly};

    # Build query
    my ($asm1, $asm2) = $self->_get_assemblies_from_table_name($table);
    my $columns;
    if(!defined($assembly) or $assembly eq $asm1) {
	$columns = "id, chr1, start1, end1, chr2, start2, end2, strand, similarity";
    }
    elsif($assembly eq $asm2) {
	($asm1, $asm2) = ($asm2, $asm1);
	$columns = "id, chr2, start2, end2, chr1, start1, end1, strand, similarity";
    }
    else {
	croak "Table $table does not have CNEs for assembly $assembly";
    }
    my $query = "SELECT $columns FROM $table";

    # Prepare and execute query
    my $sth = $self->dbh->prepare($query) or croak "Failed to prepare query [$query]";
    $sth->execute() or croak "Failed to execute query [$query]";

    # Fetch data and create CNE object
    my $cnes = $self->_get_cnes_from_sth($asm1, $asm2, $table, $sth);

    return $cnes;
}


sub load_cigar_for_cne
{
    my ($self, $cne) = @_;
    my $table = $cne->table_name or croak "No table name for cne";
    my $id = $cne->id or croak "No id for cne";
    my $query = "select cigar from $table where id = ?";
    my $sth = $self->dbh->prepare($query) or croak "Failed to prepare query [$query]";
    $sth->execute($id) or croak "Failed to execute query [$query]";
    my ($cigar) = $sth->fetchrow_array();
    $sth->finish;
    croak "No cigar for cne $id in table $table" unless($cigar);
    $cne->cigar($cigar);
    $cne->flip_cigar if($cne->assembly1 gt $cne->assembly2);
}


sub _get_cnes_from_sth
{
    my ($self, $asm1, $asm2, $table, $sth) = @_;
    my @cnes;
    while(my ($id, $chr1, $start1, $end1, $chr2, $start2, $end2, $strand, $sim) = $sth->fetchrow_array) {
	my $cne = CNE::CNE->new(db => $self,
				table_name => $table,
				id => $id,
				assembly1 => $asm1,
				assembly2 => $asm2,
				chr1 => $chr1, 
				chr2 => $chr2,
				start1 => $start1+1,
				start2 => $start2+1,
				end1 => $end1,
				end2 => $end2,
				strand => $strand,
	                        similarity => $sim);
	push @cnes, $cne;
    }
    return \@cnes;
}


=head2 count_cnes_in_table

 Title     : count_cnes_in_table
 Usage     : my $n = $db->count_cnes_in_table($table_name);
 Returns   : The number of CNEs in the table.
 Args      : Name of CNE table.    

=cut

sub count_cnes_in_table {
    my $self = shift;
    my $table_name = shift;
    croak "need table name" unless($table_name);
    return $self->_num_records($table_name);
}


sub cne_table_exists {
    my ($self, $table_name) = @_;
    return $self->_table_exists($table_name);
}


sub drop_cne_table {
    my ($self, $table_name) = @_;
    $self->dbh->do("drop table $table_name") or croak "Failed to remove table $table_name\n";
}

=head2 create_cne_table

 Title     : create_cne_table
 Usage     : $db->create_cne_table($table_name);
 Returns   : True
 Args      : Name of CNE table to create.      

=cut

sub create_cne_table {
    my ($self, $table_name) = @_;
    croak "Table $table_name already exists" if($self->_table_exists($table_name));
    my $st = "create table $table_name (
              id mediumint unsigned not null auto_increment primary key,
              bin1 smallint unsigned not null,
              chr1 varchar(255) not null,
              start1 int unsigned not null,
              end1 int unsigned not null,
              bin2 smallint unsigned not null,
              chr2 varchar(255) not null,
              start2 int unsigned not null,
              end2 int unsigned not null,
              strand char(1) not null,
              similarity float unsigned not null,
              cigar text not null)";
    $self->dbh->do($st) or croak "create table statement $st failed\n";
    return 1;
}


=head2 index_cne_table

 Title     : index_cne_table
 Usage     : $db->index_cne_table($table_name);
 Returns   : True
 Args      : Name of CNE table to index.      

=cut

sub index_cne_table {
    my ($self, $table_name) = @_;
    my $dbh = $self->dbh;
    unless($dbh->do("alter table $table_name add key (chr1, bin1, start1)") and
	   $dbh->do("alter table $table_name add key (chr1, bin1, end1)") and
	   $dbh->do("alter table $table_name add key (chr2, bin2, start2)") and
	   $dbh->do("alter table $table_name add key (chr2, bin2, end2)")) {
	die "Alter table add key statement failed\n";
    }
    return 1;
}


=head2 create_cne_importer

 Title     : create_cne_importer
 Usage     : my $importer = $db->create_cne_importer($table_name);
 Function  : Creates an AT::DB::Importer object that allows effective
             import of a large number of CNEs to the database.
             See the AT::DB::Importer module for further details.
 Example   : my $importer = $db->create_cne_importer();
             while (my $cne = $cne_stream->next_cne) {
		 $importer->store_object($cne);
	     }
             my @result = $importer->finish;
 Returns   : An AT::DB::Importer object
 Args      : Name of CNE table to import into.      

=cut

sub create_cne_importer {
    my ($self, $table_name) = @_;

    # define wrapper for storing one cne
    my $storing_wrapper = sub {
	my ($importer, $cne) = @_; 
	croak "Wrong type of CNE object" unless $cne->isa("CNE::CNE");
	croak "Missing data in CNE object" unless $cne->validate;
	$cne = $cne->clone->swap_locations if($cne->assembly1 gt $cne->assembly2);
	my $start1 = $cne->start1;
	my $end1 = $cne->end1;
	my $start2 = $cne->start2;
	my $end2 = $cne->end2;
	my $bin1 = $self->bin_from_coord_range($start1, $end1);
	my $bin2 = $self->bin_from_coord_range($start2, $end2);
	$importer->store_record($table_name,
				[$bin1, $cne->chr1, $start1-1, $end1,
				 $bin2, $cne->chr2, $start2-1, $end2,
				 $cne->strand, $cne->similarity, $cne->cigar]);
   };

    # prepare table info needed by importer
    my $tables_columns = {
	$table_name => [ qw(bin1 chr1 start1 end1 bin2 chr2 start2 end2 strand similarity cigar) ]
    };

    # create importer object; tell it to lock database
    my $importer = AT::DB::Importer->new( db => $self,
					  tables_columns => $tables_columns,
					  storing_wrapper => $storing_wrapper);

    # return the creation
    return $importer;    
}


=head2 store_cne

 Title     : store_cne
 Usage     : $db->store_cne($cne);
 Function  : Stores a CNE in the database
 Returns   : 1 on success
 Args      : A CNE::CNE object
 Note      : Not yet implemented.

=cut


sub store_cne {
    my ($self, $cne) = @_;

    # Check argument
    unless ($cne->isa("CNE::CNE")) { 
        croak "Wrong type of CNE object";
    }
    
    die "not implemented";

    # Insert data into table

    return 1;
}


sub get_assembly_info {
    my ($self, $asm_id) = @_;
    croak "Need assembly id as first argument" unless($asm_id);
    my $query = "select * from assembly where assembly_id = ?";
    my $sth = $self->dbh->prepare($query) or croak "Failed to prepare query [$query]";
    $sth->execute($asm_id) or croak "Failed to execute query [$query]";
    my $asm_info = $sth->fetchrow_hashref;
    $sth->finish;
    return $asm_info ? { %$asm_info } : undef;
    # ^ return a reference to a copy of the hash created by DBI in case DBI reuses that hash
}


sub DESTROY  {
    my $self = shift;
    $self->dbh->disconnect() if $self->dbh();
}

1;


