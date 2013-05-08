# AT::GFX::LocusConfig::Assembly module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::GFX::LocusConfig::Assembly - assembly configuration for drawing loci

=head1 SYNOPSIS


=cut


package AT::GFX::LocusConfig::Assembly;

use strict;
use vars '@ISA';
use Carp;
use Class::Struct;
use AT::Root;
use AT::DB::GenomeMapping;
use AT::DB::GenomeAlignment;
use AT::DB::GenomeAssembly;
use AT::DB::GenomeAssemblyNibs;
use AT::DB::GenomeAssemblyTwoBit;


@ISA = qw(AT::Root);


sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
	id => $args{id},
	name => $args{name} || $args{id},
	seq_2bit => $args{seq_2bit},
	seq_nib => $args{seq_nib},
	seq_db_name => $args{seq_db_name},
	seq_db_conn => $args{seq_db_conn},
	AT_db_name => $args{AT_db_name},
	AT_db_conn => $args{AT_db_conn},
	track_db_name => $args{AT_db_name},
	track_db_conn => $args{AT_db_conn},
	tracks => $args{tracks} || [],
	enabled_tracks => $args{enabled_tracks} || {}
    }, ref $caller || $caller;

    return $self;   
}


sub get_drawing_parameters
{
    my ($self) = @_;

    my %result;

    # connect to AT database
    my $at_db_param = $self->AT_db_conn or die "No AT database connection for assembly ".($self->id);
    my $aln_db = AT::DB::GenomeAlignment->connect(-dbname => $self->AT_db_name,
						  -dbhost => $at_db_param->host,
						  -dbuser => $at_db_param->user,
						  -dbpass => $at_db_param->pass);
    $result{"alignment_db"} = $aln_db;

    # connect to UCSC database
    my $ucsc_db_param = $self->track_db_conn or die "No track database connection for assembly ".($self->id);
    my $ucsc_db = AT::DB::GenomeMapping->connect(-dbname => $self->track_db_name,
						 -dbhost => $at_db_param->host,
						 -dbuser => $at_db_param->user,
						 -dbpass => $at_db_param->pass);
    $result{"track_db"} = $ucsc_db;

    # connect to sequence database, if any
    my $seq_db;
    if($self->seq_db_name) {
	my $seq_db_param = $self->seq_db_conn;
	$seq_db = AT::DB::GenomeAssembly->connect(-dbname => $self->seq_db_name,
						  -dbhost => $seq_db_param->host,
						  -dbuser => $seq_db_param->user,
						  -dbpass => $seq_db_param->pass);
    }
    elsif($self->seq_2bit) {
	$seq_db = AT::DB::GenomeAssemblyTwoBit->new(file => $self->seq_2bit,
						    id => $self->id);
    }
    elsif($self->seq_nib) {
	$seq_db = AT::DB::GenomeAssemblyNibs->new(dir => $self->seq_2bit,
						  id => $self->id);
    }
    $result{assembly_db} = $seq_db;

    # Set organism and tracks parameters
    $result{organism} = $self->name,
    my $enabled_tracks = $self->enabled_tracks;
    $result{tracks} = [grep { $enabled_tracks->{$_->{name}} } @{$self->tracks}];
    
    return \%result;
}


1;
