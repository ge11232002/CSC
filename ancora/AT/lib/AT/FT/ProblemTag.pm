#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::FT::GFMappingHSP

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::FT::GFMappingHSP;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;

@ISA = qw/AT::Root/;


sub new  {
    my ($caller, %args) = @_;
    my $self = bless {
	start => ($args{start} || croak "No start arg"),
	end => ($args{end} || croak "No end arg"),
	left_gap => ($args{left_gap} || undef),
	right_gap => ($args{right_gap} || undef),
	next => ($args{'next'} || 0),
	#open_pep => ($args{open_pep} || 0),  # <-- keep?
	#close_pep => ($args{close_pep} || 0) # <-- kepp?
    }, ref $caller || $caller;
    return $self;
}


sub clone {
    my ($self) = @_;
    return $self->new(start => $self->start,
		      end => $self->end);
}


sub length {
    my ($self) = @_;
    return $self->end - $self->start + 1;
}

1;

