# AT::DB::Collection module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::DB::Collection

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::DB::Collection;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;

@ISA = qw(AT::Root);

sub new
{
    my ($class, %args) = @_;
    my $self = bless {
	_dbL => $args{db_list} || [],
	db_selection_subs => $args{db_selection_subs} || {}
    }, ref($class) || $class;
    $self->_init();
    return $self;
}


sub _init
{
}

sub add_db_selection_sub
{
    my ($self, $method, $sub) = @_;
    $self->{db_selection_subs}{$method} = $sub;
}


sub _select_db {
    my ($self, $method, @args) = @_;
    my $sub = $self->{db_selection_subs}->{$method};
    return $sub ? $sub->(@args) : undef;
}


sub _generic_db_call {
    my ($self, $method, @args) = @_;

    my $db = $self->_select_db($method, @args);
    if($db) {
	unless(ref $db) {
	    die "Cannot handle db names yet";
	}
        return $db->$method(@args);
    }
    else {
        foreach my $db ($self->db_list) {
	    my @result = $db->$method(@args);
	    return @result if(@result);
        }
    }
    return;
}


1;


