#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::Mapping - mapping of a transcript sequence to genomic sequence

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::FT::CloneMapping;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;


@ISA = qw/AT::Root/;

use Class::Struct '_clone_mapping_part' => [
    start => '$',
    end => '$',
    start_type => '$',
    end_type => '$',
    next => '_clone_mapping_part',
    open_pep => '_pep',   # <-- keep?
    close_pep => '_pep'   # <-- keep?
];


sub new  {
    my ($caller, %args) = @_;
    my $self = bless {
	_mappingL => ($args{mappings} || croak 'No mappings'),
	_partL => [],
	grouped => 0			# <-- remove this sometime
    }, ref $caller || $caller;
    return $self;
}

sub part { $_[0]->{_partL}->[$_[1]]; }
sub nr_parts { scalar($_[0]->{_partL}); }

# provide: start, end, startType, endType
sub add_part {
    my ($self, %args) = @_;
    my $part = _clone_mapping_part->new(start_type => '',
					end_type => '',
					%args);
    $self->{_partL}->[-1]->next($part) if(scalar @{$self->{_partL}});
    push @{$self->{_partL}}, $part;
    return $part;
}

1;

