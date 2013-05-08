# AT::DB::ResultStream module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::DB::ResultStream

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::DB::ResultStream;

use strict;
use vars '@ISA';
use Carp;

@ISA = qw(AT::Root);


sub new
{
    my ($caller, %args) = @_;   
    my $self = bless { _next_method => $args{next_method} }, ref($caller) || $caller;
    return $self;
}

sub next
{
    my ($self) = @_;
    my $next_method = $self->{_next_method};
    &$next_method();
}

1;


