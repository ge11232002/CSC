package AT::HSP;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;


@ISA = qw(AT::Root);

sub new  {
    my ($caller, %args) = @_;
    my $self = bless { hsp_id => undef,
		       mapping_id  => undef,
		       qStart  => undef,
		       qEnd  => undef,
		       tStart  => undef,
		       tEnd  => undef,
		       qSeq  => undef,
		       tSeq => undef,
		       lFlank => undef,
		       rFlank => undef,
			   %args}, ref $caller || $caller;
    $self->{'blockSize'} = length($self->qSeq)
	unless defined($self->{'blockSize'});
    return $self;
}


sub qPos_from_tPos_exact {
    my ($self, $tPos, %args) = @_;
    return undef unless ($tPos >= $self->tStart and $tPos <= $self->tEnd);
    my $qPos = $tPos - $self->tStart + $self->qStart;
    #print STDERR "HSP: tPos = $tPos; qPos = $qPos\n";
    return $qPos;
}

sub tPos_from_qPos_exact {
    my ($self, $qPos, %args) = @_;
    return undef unless ($qPos >= $self->qStart and $qPos <= $self->qEnd);
    my $tPos = $qPos - $self->qStart + $self->tStart;
    #print STDERR "HSP: tPos = $tPos; qPos = $qPos\n";
    return $tPos;
}


1;
