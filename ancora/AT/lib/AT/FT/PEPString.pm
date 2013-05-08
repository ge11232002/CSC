# AT::FT::PEPString module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME


=head1 SYNOPSIS



=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::FT::PEPString;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;


@ISA = qw/AT::Root/;


sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
	%args
    }, ref $caller || $caller;
   
    return $self;
}


sub delete
{
    my ($string) = @_;
    my $l = $string->left_pep;
    my $r = $string->right_pep;
    delete $l->right_strings->{$r};
    delete $r->left_strings->{$l};
}


1;
