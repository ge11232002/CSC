# AT::DB::DataSourceHolder module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::DB::DataSourceHolder

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::DB::DataSourceHolder;

use strict;
use vars '@ISA';
use Carp;
use AT::DB::MySQLdb;
use AT::DB::DataSource;

@ISA = qw(AT::DB::MySQLdb);


=head2 add_data_source

 Title     : add_data_source
 Usage     : 
 Function  : 
 Returns   : 
 Args      : 

=cut


sub add_data_source  {
    my ($self, $source) = @_;

    # Check argument
    croak "Wrong data source object type" 
        unless (ref($source) and $source->isa('AT::DB::DataSource'));

    # Store object
    my $query = "INSERT DATA_SOURCE (name, version, component, ".
		"load_start_time) VALUES (?, ?, ?, NOW())";
    $self->dbh->do($query, undef, $source->name, $source->version,
		   $source->component)
	or return undef;
    $source->data_source_id($self->_last_insert_id);
}


=head2 get_data_source

 Title     : get_seq
 Usage     : my $seq = $db->get_seq("NM_001243");
 Function  : Gets a sequence from the database by accession number.
 Returns   : A Bio::Seq object with sequence and description on success.
             Undef if the accession number is not in the database.
 Args      : An accession number.

=cut


sub get_data_source {
    my ($self, $name, $version, $component) = @_;

    unless(defined($name) and defined($version) and defined($component)) {
	croak "Data source name, version and/or component undefined";
    }
    
    my $data = $self->dbh->selectrow_hashref(
	"SELECT * FROM DATA_SOURCE WHERE name = ? AND version = ? AND ".
	"component = ?", undef, $name, $version, $component);
    return undef unless ($data);
    return AT::DB::DataSource->new(%$data);
}



sub set_data_source_load_end_time {
    my ($self, $ds, $name, $version, $component);
    if(@_ == 2) {
	($self, $ds) = @_;
	$name = $ds->name,
	$version = $ds->version;
	$component = $ds->component;
    }
    elsif(@_ == 4)  {
	($self, $name, $version, $component) = @_;
    }
    else { croak 'Invalid number of arguments'; }

    unless(defined($name) and defined($version) and defined($component)) {
	croak "Data source name, version and/or component undefined";
    }
    
    return $self->dbh->do(
	"UPDATE DATA_SOURCE SET load_end_time = NOW() ".
	"WHERE name = ? AND version = ? AND component = ?",
	undef, $name, $version, $component);
}


1;


