# AT::GFX::LocusConfig module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::GFX::LocusConfig - configuration for drawing loci

=head1 SYNOPSIS


=cut


package AT::GFX::LocusConfig;

use strict;
use vars '@ISA';
use Carp;
use Class::Struct;
use AT::GFX::LocusConfig::Assembly;


@ISA = qw(AT::Root);


struct my_conf_sql_conn => [
    id => '$',
    host => '$',
    user => '$',
    pass => '$'
];


sub new
{
    my ($caller, %args) = @_;

    my $self = bless {}, ref $caller || $caller;
    my $filename = $args{file} || croak("No file argument");
    $self->_read_config($filename);

    return $self;   
}


sub _read_config
{
    my ($self, $fn) = @_;
    my (%conn_index, %asm_index);
    open IN, $fn or croak("could not open config file $fn");
    my $block_type = 0;
    my $block_id;
    while(my $line = <IN>) {
	chomp $line;
	my ($param, @data) = split /\s+/, $line;
	next if(!$param or $param =~ /^\#/);
	if($param eq 'sqlConnection') {
	    my ($id, $host, $user, $pass) = @data;
	    croak("No id for sqlConnection in $fn") unless($id);
	    croak("No host for sqlConnection in $fn") unless($host);
	    my $conn = my_conf_sql_conn->new(id => $id,
					host => $host,
					user => $user,
					pass => $pass);
	    $conn_index{$id} = $conn;
	    $block_type = 0;
	}
	elsif($param eq 'assembly') {
	    $block_type = 'assembly';
	    $block_id = $data[0];
	    croak("No id for assembly in $fn") unless($block_id);
	    $asm_index{$block_id} = AT::GFX::LocusConfig::Assembly->new(id => $block_id);
	}
	elsif($param eq 'name') {
	    croak("$param line outside assembly block in $fn") unless($block_type eq 'assembly');
	    my $name = $data[0];
	    $name =~ s/_/ /g;
	    die "No data for $param in $fn" unless($name);
	    $asm_index{$block_id}->name($name);
	}
	elsif($param eq 'seq_nib') {
	    croak("$param line outside assembly block in $fn") unless($block_type eq 'assembly');
	    my $id = $data[0];
	    die "No data for $param in $fn" unless($id);
	    $asm_index{$block_id}->seq_nib($id);
	}
	elsif($param eq 'seq_2bit') {
	    croak("$param line outside assembly block in $fn") unless($block_type eq 'assembly');
	    my $id = $data[0];
	    croak("No data for $param in $fn") unless($id);
	    $asm_index{$block_id}->seq_2bit($id);
	}
	elsif($param eq 'seq_db') {
	    croak("$param line outside assembly block in $fn") unless($block_type eq 'assembly');
	    my ($db_name, $conn_id) = @data;
	    croak("No database name for $param in $fn") unless($db_name);
	    croak("No connection id for $param in $fn") unless($conn_id);
	    my $conn = $conn_index{$conn_id};
	    croak("Undefined connection $conn_id for $param $db_name in $fn") unless($conn);
	    $asm_index{$block_id}->seq_db_name($db_name);
	    $asm_index{$block_id}->seq_db_conn($conn);
	}
	elsif($param eq 'trackDb') {
	    croak("$param line outside assembly block in $fn") unless($block_type eq 'assembly');
	    my ($db_name, $conn_id) = @data;
	    croak("No database name for $param in $fn") unless($db_name);
	    croak("No connection id for $param in $fn") unless($conn_id);
	    my $conn = $conn_index{$conn_id};
	    croak("Undefined connection $conn_id for $param $db_name in $fn") unless($conn);
	    $asm_index{$block_id}->track_db_name($db_name);
	    $asm_index{$block_id}->track_db_conn($conn);
	}
	elsif($param eq 'AT') {
	    croak("$param line outside assembly block in $fn") unless($block_type eq 'assembly');
	    my ($db_name, $conn_id) = @data;
	    croak("No database name for $param in $fn") unless($db_name);
	    croak("No connection id for $param in $fn") unless($conn_id);
	    my $conn = $conn_index{$conn_id};
	    croak("Undefined connection $conn_id for $param $db_name in $fn") unless($conn);
	    $asm_index{$block_id}->AT_db_name($db_name);
	    $asm_index{$block_id}->AT_db_conn($conn);
	}
	elsif($param eq 'trackTable') {
	    croak("$param line outside assembly block in $fn") unless($block_type eq 'assembly');
	    my ($id, $table, $format, $glyph, $bgcolor, $fgcolor, $enabled) = @data;
	    croak("Missing data for line $line in $fn") unless($id and $table and $format and $glyph and $bgcolor and $fgcolor);
	    $id =~ s/_/ /g;
	    my $tracks = $asm_index{$block_id}->tracks;
	    push @$tracks, { 'name' => $id,
			     'type' => 'table',
			     'tables' => [split /,/, $table],
			     'format' => $format,
			     'glyph' => $glyph,
			     'bgcolor' => $bgcolor,
			     'fgcolor' => $fgcolor
			     };
	    $asm_index{$block_id}->enabled_tracks->{$id} = 1 if($enabled);
	}
	elsif($param eq 'trackFile') {
	    croak("$param line outside assembly block in $fn") unless($block_type eq 'assembly');
	    my ($id, $trackfile, $format, $glyph, $bgcolor, $fgcolor, $enabled) = @data;
	    croak("Missing data for line $line in $fn") unless($id and $trackfile and $format and $glyph and $bgcolor and $fgcolor);
	    $id =~ s/_/ /g;
	    my $tracks = $asm_index{$block_id}->tracks;
	    push @$tracks, { 'name' => $id,
			     'type' => 'file',
			     'filename' => $trackfile,
			     'format' => $format,
			     'glyph' => $glyph,
			     'bgcolor' => $bgcolor,
			     'fgcolor' => $fgcolor
			     };
	    $asm_index{$block_id}->enabled_tracks->{$id} = 1 if($enabled);
	}
	else {
	    croak("unknown parameter $param in $fn");
	}
    }
    
    $self->{sql_conn} = \%conn_index;
    $self->{assemblies} = \%asm_index;
}


1;
