# AT::DB::Importer module
#
# Copyright Boris Lenhard
# You may distrubute this code under the same terms as perl itself
#
# POD documentation:

=head1 NAME

AT::DB::Importer - handles batch import into a MySQL database.

=head2 SYNOPSIS

# Ask your favourite database object to create an importer for some object it can handle

my $importer = $mapping_db->create_mapping_importer;

# Feed the importer some objects; they will be stored in a temp file

while(my $mapping = $mapping_stream->next_mapping) {
    $importer->store_object($mapping);
}

# Have your objects loaded into the database in one batch

$importer->finish;

=head2 DESCRIPTION

Inserting many records into a database by successive INSERT queries is inefficient due to query overhead, repeated reindexing etc. This module provides an alternative and often more efficient approach. The records to be inserted are first written to temporary files. One file is used for each database table. A LOAD DATA query is then issued for each table, asking the database to import the corresponding temporary file. The main drawback is that when data is imported to several tables, the tables may have to be locked for the lifespan of the importer.

=cut

# The code begins HERE

package AT::DB::Importer;

use strict;
use vars '@ISA';
use Carp;
use File::Temp qw( tempfile );
use AT::Root;

@ISA = qw(AT::Root);



sub new  {
    my ($caller, %args) = @_;

    # Check required parameters
    my $tables_columns =  $args{tables_columns} or
	croak ("No tables specified");
    my $db = $args{db} or
	croak ("No db specified.");
    my $storing_wrapper = $args{storing_wrapper} or
	croak("No storing wrapper specified");

    # Prepare table info data structure
    my %table_info;
    while(my ($table_name, $columns) = each %$tables_columns) {
	$table_info{$table_name} = { columns => $columns,
				     num_records => 0 };
    }

    # Make self
    my $self = bless {_table_info => { %table_info },
		      _db => $db,
		      _storing_wrapper => $storing_wrapper,
		      _lock_on => 0,
		      _cleanup => 1,
		      counter => $args{counter} || 1
		      }, ref $caller || $caller;

    $self->_init(lock => $args{lock} || 0);

    return $self;
}


sub _init {
    my ($self, %args) = @_;

    # Open a tmp file for each table
    while(my ($table_name, $table_info) = each %{$self->{_table_info}}) {
	my ($fh, $fn) = tempfile();
	$table_info->{fh} = $fh;
	$table_info->{filename} = $fn;
	print STDERR "Importer: writing data for import into $table_name to temporary file $fn\n";
    }
    
    # Lock tables if told to
    # Note: there are some safety flaws here. If the same db object makes any
    # db queries or creates another Importer before this one is finished,
    # the locking might be cancelled (haven't tested this).
    if ($args{lock}) {
	my @table_names = (keys %{$self->{_table_info}});
	$self->_lock_tables( \@table_names,
			     [ ("WRITE") x scalar(@table_names) ] )
	    or do { warn "Importer: could not lock tables.";
		    return undef; }
    }
}


# Store one object. The storing wrapper arranges the information in the object
# as records and calls store_record (below) to write each record to a
# temp file.

sub store_object {
    my ($self, $object) = @_;
    my $storing_wrapper = $self->{_storing_wrapper};
    &$storing_wrapper($self, $object);
}

# Write one record (one line) to a temp file.
# This method should only be called by the storing wrapper

sub store_record {
    my ($self, $table, $values) = @_;
    
    my $table_info = $self->{_table_info}->{$table};

    # Replace undefined values with MySQL NULL ('\N')
    my @my_values;
    foreach my $value (@$values) {
	push @my_values, defined($value) ? $value : '\N';
    }

    # Write line to tempfile
    my $fh = $table_info->{fh};
    print $fh (join "\t", @my_values), "\n";

    # Increment record counter
    $table_info->{num_records}++;
}


# Load files > database

sub finish {
    my $self = shift;
    my @num_records_loaded;

    # Since parsing can take long, we must check that we are still connected
    my $db = $self->{_db};
    unless($db->dbh->ping) {
	warn "Importer: disconnected from database; trying to reconnect...\n";
	unless($db->reconnect()) {
	    warn "Importer: failed to reconnect to database\n";	
	    return undef;
	}
	warn "Importer: reconnected to database.\n";
    }

    # For each table/tmpfile
    while(my ($table_name, $table_info) = each %{$self->{_table_info}}) {

	# Close the file
	close $table_info->{fh};
	undef $table_info->{fh};

	# Get the filename, and undef it to make sure we don't do anything further with the file.
	# This means the file will be kept unless it is unlinked in this scope.
	my $filename = $table_info->{filename};
	undef $table_info->{filename};

	# Import the file > database
	my $num_records_loaded = $self->_load_data($filename,
						   $table_name,
						   $table_info->{columns});
	if (!defined($num_records_loaded)) {
	    warn "Importer: error executing LOAD DATA query with file $filename (file kept)";
	    return undef;
	}
	elsif ($num_records_loaded != $table_info->{num_records}) {
	    warn ("Importer: loaded ".(0+$num_records_loaded)." out of ".
		  $table_info->{num_records}." records from file ".
		  $filename." into table ".
		  $table_name." (skipped records possibly duplicate values ".
		  "on a unique index; file kept)");
	}
	else {
	    # Remove the file if import went well
	    unlink $filename;
	}

	push @num_records_loaded, $num_records_loaded;
    }

    # Unlock tables
    $self->_unlock_tables if $self->{_lock_on};

    return (@num_records_loaded);
}


sub disable_cleanup
{
    my $self = shift;
    $self->{_cleanup} = 0;
}


sub _load_data
{
    my ($self, $filename, $table, $columns) = @_;
    my $query = "LOAD DATA LOCAL INFILE \"$filename\" INTO TABLE $table (".
	    (join ',', @$columns) . ')';
    my $dbh = $self->{_db}->dbh;
    my $nr_rows = $dbh->do($query);
    print STDERR "Importer: ", $dbh->errstr, "\n" unless(defined $nr_rows);
    return $nr_rows;
}


sub _lock_tables {
    my ($self, $tables, $modes) = @_;
    
    my @lock_items;
    foreach my $i (0 .. @$tables-1) {
	my $table = $tables->[$i];
	my $mode = $modes->[$i];
	unless ($mode =~ /^(read|write)$/i) {
	    warn "Importer: must specify read or write mode for each table";
	    return undef; }
	push @lock_items, "$table $mode";
    }

    my $query = "LOCK TABLES " . (join ',', @lock_items);

    $self->{_db}->dbh->do($query) or return undef;

    $self->{_lock_on} = 1;
    return 1;
}


sub _unlock_tables {
    my $self = shift;
    $self->{_db}->dbh->do("UNLOCK TABLES") or return undef;
    $self->{_lock_on} = 0;
    return 1;
}


sub DESTROY
{
    my $self = shift;

    my $cleanup = $self->{_cleanup};

    # Close and remove temporary files
    while(my ($table_name, $table_info) = each %{$self->{_table_info}}) {
	close $table_info->{fh} if defined ($table_info->{fh});
	unlink $table_info->{filename} if($cleanup and defined($table_info->{filename}));
    }

    # Unlock tables
    $self->_unlock_tables if($cleanup and $self->{_lock_on});
}

1;
