#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::mRNAInfo - supplementary information about a transcript sequence

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::mRNAInfo;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;

@ISA = qw(AT::Root);

=head2 new

 Title     : new
 Usage     : 
 Function  : Constructor
 Returns   : AT::mRNAInfo
 Args      : 

=cut

sub new  {
    my ($caller, %args) = @_;
    my $self = bless {
	acc => $args{acc},
	version => $args{version},
	type => $args{type},
	library => $args{library},
	mrnaClone => $args{mrnaClone},
	cds => $args{cds},
	refSeqStatus => ($args{refSeqStatus} || '')
    }, ref $caller || $caller;

    return $self;
}


1;

