#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

package AT::Tools::RandomGenomicRegionSelector;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;


@ISA = qw(AT::Root);


sub new  {
    my ($caller, %args) = @_;
    my $self = bless {
	_db => ($args{'db'} || croak 'no db argument'),
	length => ($args{'length'} || 1),
	include_chrN_random => ($args{'include_chrN_random'} || 0)
    }, ref $caller || $caller;
    $self->_init();
    return $self;
}


sub _init {
    my ($self) = @_;
    my @chr_names = $self->{_db}->get_chr_names;
    unless($self->include_chrN_random) {
	@chr_names = grep { !($_ =~ /random$/) } @chr_names;
    }
    my @max_starts = map { $self->{_db}->get_chr_size($_) } @chr_names;
    my $max_rnd = 0;
    for my $i (0..@max_starts-1) {
	$max_starts[$i] -= $self->{length} + 1;
	croak "Some chromosome is shorter that requested region length"
	    if($max_starts[$i] < 1);
	$max_rnd += $max_starts[$i];
    }
    $self->{_chr_names} = \@chr_names;
    $self->{_max_starts} = \@max_starts;
    $self->{_max_rnd} = $max_rnd;
}


sub select_one_region {
    my ($self) = @_;
    my $rnd = int(rand($self->{_max_rnd}))+1;
    my $max_starts = $self->{_max_starts};
    my $i;
    for($i = 0; $rnd > $max_starts->[$i]; $i++) {
	$rnd -= $max_starts->[$i];
    }
    return($self->{_chr_names}->[$i],
	   $rnd,
	   $rnd + $self->{length} - 1);
}

1;

