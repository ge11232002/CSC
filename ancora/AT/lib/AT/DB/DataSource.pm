# AT::DB::DataSource module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::DB::DataSource

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::DB::DataSource;

use strict;
use vars '@ISA';
use DBI;
use Carp;

@ISA = qw(AT::Root);


=head2 new

 Title     : new
 Usage     : 
 Function  : Constructor
 Returns   : 
 Args      : 

=cut

sub new {
    my ($caller, %args) = @_;
    my $self = bless { data_source_id => undef,
		       name => undef,
		       version => undef,
		       component => undef,
		       load_start_time => undef,
		       load_end_time => undef,
		       %args }, ref $caller || $caller;
    return $self;
}


sub id_str {
    my ($self) = @_;
    return $self->name.':'.$self->version.':'.$self->component;
}

1;


