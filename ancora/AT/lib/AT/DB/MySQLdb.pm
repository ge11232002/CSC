# AT (Alternative Transcripts) module for AT::MySQLdb
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::DB::MySQLdb - a superclass for the interface to MySQL relational databases 

=head1 SYNOPSIS

You should not use this class directly. Use the subclasses
AT::DB::GenomeAssembly and AT::DB::GenomeMapping

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::DB::MySQLdb;

use strict;
use vars '@ISA';
use DBI;
use Carp;
use AT::Root;
use Term::ReadKey;

@ISA = qw(AT::Root);


=head2 connect

 Title     : connect
 Usage     : my $db = AT::DB::GenomeMapping->connect(-dbname => "AT",
						     -dbhost => "myhostn.mydomain",
						     -dbuser => "myusername",
						     -dbpass => "mypassword",
						     -dbport => 3306);
 Function  : create a database object by connecting to an existing AT-type
             database
 Returns   : A database object
 Args      : All arguments are optional.
             dbname  Database name (default: test)
             dbhost  Host running the MySQL server (default: localhost)
             dbuser  MySQL username (no default)
             dbpass  MySQL password (no default);
                     If set to '?', the password is prompted for on console
                     If set to '??', emtpy inputs cause repeated prompting
             dbport  MySQL port number (no default)

=cut


sub connect  {
    my ($caller, %args) = @_;

    # Get parameters
    my $dbhost = $args{'-dbhost'} || $args{'dbhost'} || 'mysql-cbu.bccs.uib.no'; 
    my $dbname = $args{'-dbname'} || $args{'dbname'} || 'test'; 
    my $dbuser = $args{'-dbuser'} || $args{'dbuser'} || undef; 
    my $dbpass = $args{'-dbpass'} || $args{'dbpass'} || undef; 
    my $dbport = $args{'-dbport'} || $args{'dbpost'} || undef;
    ($dbhost, $dbport) = split(/\:/, $dbhost) unless(defined $dbport);
 
    # Create object
    my $self = bless {
	dbhost => $dbhost,
	dbname => $dbname,
	dbuser => $dbuser,
	dbpass => $dbpass,
	dbport => $dbport
    }, ref $caller || $caller;

    # Ask for password if password was given as '?'
    if(defined $dbpass) {
	if($dbpass eq '?') {
	    $self->_ask_password();
	}
	elsif($dbpass eq '??') {
	    $self->_ask_password("",1);
	}
    }

    # Connect to db
    $self->reconnect() or croak "Could not connect to mysql database.";

    # Initialize
    $self->_init(%args);

    return $self;
}


=head2 reconnect

 Title     : reconnect
 Usage     : $db->reconnect();
 Function  : Reconnects to database using parameters from earlier
             call to connect(). Can be used if db connection was lost.
 Returns   : True on success; false on failure.
 Args      : None

=cut


sub reconnect  {
    my ($self) = @_;

    # Get connection parameters
    my $dbhost = $self->{dbhost} or croak "No dbhost set";
    my $dbname = $self->{dbname} or croak "No dbname set";
    my $dbuser = $self->{dbuser};
    my $dbpass = $self->{dbpass};
    my $dbport = $self->{dbport};
    my $data_source = "dbi:mysql:database=$dbname;host=$dbhost"; 
    $data_source .= ";port=$dbport" if($dbport);

    # Connect to db
    my $dbh = DBI->connect($data_source, $dbuser, $dbpass);
    $self->{_dbh} = $dbh;

    # Return
    return $dbh ? 1 : 0;
}


=head2 ping

 Title     : ping
 Usage     : $db->reconnect unless($db->ping())
 Function  : Check if database connection is working.
 Returns   : True if there seems to be a working connection, otherwise false
 Args      : None

=cut


sub ping {
    my ($self) = @_;
    return $self->{_dbh}->ping;
}


=head2 dbh

 Title     : dbh
 Usage     : my $dbh = $db->dbh();
 Function  : Retrieve DBI database handle, which can be used for
             direct queries the database (not recommended).
 Returns   : A database handle object.
 Args      : None.
 See also  : perldoc DBI

=cut


sub dbh {
    # method should be read-only
    return $_[0]->{'_dbh'}
}


sub _init  {
    # intentionally left blank - to be overriden in subclasses if necessary   
}


sub _ask_password {
    my ($self, $question, $repeat) = @_;
    $question = "Password for database ".($self->{dbname}).": "
	unless($question);
    my $answer;
    ReadMode('noecho');
    do {
	print $question;
	$answer = ReadLine(0);
	chomp $answer;
	print "\n";
    } until ($answer or !$repeat);
    ReadMode('restore');
    $self->{dbpass} = $answer;
    return $answer;
}


sub _do_or_die {
    my $self = shift;
    $self->dbh->do(@_) or die($self->dbh->errstr);
}

sub _last_insert_id  {
    my ($self) = @_;
    my $sth = $self->dbh->prepare_cached("SELECT LAST_INSERT_ID()");
    $sth->execute;
    my ($liid) = $sth->fetchrow_array;
    $sth->finish();
    return $liid;
}


sub _num_records {
    my ($self, $table) = @_;
    my $sth = $self->dbh->prepare
	("SELECT COUNT(*) FROM $table");
    $sth->execute;
    my ($num) = $sth->fetchrow_array;
    $sth->finish();
    return $num;
}


sub _get_table_names {
    my ($self, $table) = @_;
    my @tables;
    foreach my $t ($self->dbh->tables(undef, undef, undef, 'TABLE')) {
	$t =~ s/.*?\.//; # MySQL 5 returns `DBNAME`.`TABLENAME` strip everything but tablename
	push @tables, substr($t,1,length($t)-2); # chop qoutes around table name
    }
    return @tables;
}


sub _table_exists {
    my ($self, $t1) = @_;
    foreach my $t2 ($self->_get_table_names) {
	return 1 if ($t1 eq $t2);
    }
    return 0;
    #my $result = grep '^.?$table.?', $self->dbh->tables(undef, undef, undef, 'TABLE');
    #return $result;
}


sub _get_column_names {
    my ($self, $table) = @_;
    my $sth = $self->dbh->column_info(undef, undef, $table, '%');
    my @columns;
    while (my $row = $sth->fetchrow_arrayref) {
	push @columns, $row->[3];
    }
    return @columns;
}


sub _column_exists {
    my ($self, $table, $column) = @_;
    my $sth = $self->dbh->column_info(undef, undef, $table, '%');
    while (my $row = $sth->fetchrow_arrayref) {
	return 1 if($row->[3] eq $column);
    }
    return 0;
}


sub _parse_comma_separated_config_string
{
    my ($self,$string) = @_;
    my @list;
    foreach my $item (split ",", $string) {
	($item) = $item =~ /^\s*(\w+)\s*$/;
	if(!defined $item) {
	    die "Error in string $string in configuration for db ", $self->dbname;
	}
	elsif($item ne '') {
	    push @list, $item;
	}
    }
    return \@list;
}


sub DESTROY  {
    my $self = shift;
    $self->dbh->disconnect() if $self->dbh();
}

1;


