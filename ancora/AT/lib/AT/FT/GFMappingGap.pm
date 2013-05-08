#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::FT::GFMappingGap

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::FT::GFMappingGap;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;

@ISA = qw/AT::Root/;

our %_jnc_info = ( 'GTAG' => ['canonical', 1],
		    'CTAC' => ['canonical', -1],
		    'GCAG' => ['minor', 1],
		    'CTGC' => ['minor', -1],
		    'ATAC' => ['minor', 1],
		    'GTAT' => ['minor', -1] );
# ^^ note: should make canonical, minor and other into enums
# may want to create package AT::DEF::SpliceJunction for that

sub new  {
    my ($caller, %args) = @_;
    my $jnc_str = ($args{jnc_str} || '');
    my $self = bless {
	is_broken => ($args{is_broken} || 0),
	is_intron => ($args{is_intron} || 0),  # this is set by compass
	jnc_str => uc($jnc_str),
	qInsert_length => ($args{qInsert_length} || 0),
    }, ref $caller || $caller;
    $self->_set_jnc_data();
    return $self;
}


sub clone {
    my ($self) = @_;
    return $self->new(is_broken => $self->is_broken,
		      is_intron => $self->is_intron,
		      jnc_str => $self->jnc_str,
		      qInsert_length => $self->qInsert_length);
}

sub _set_jnc_data
{
    my ($self) = @_;
    my $jnc_info = $_jnc_info{$self->{jnc_str}};
    if($jnc_info) {
	$self->{jnc_type} = $jnc_info->[0];
        $self->{jnc_strand} = $jnc_info->[1];
    }
    else {
	$self->{jnc_type} = 'other';
        $self->{jnc_strand} = 0;
    }
}


sub merge {
    my ($a, $b) = @_;
    my ($jnc_str, $qInsert_length);
    if($a->is_broken) {
        croak "Conflicting gap types" unless($b->is_broken);
	# leave jnc_str and qInsert_length undefined; they do not apply to broken gaps
    }
    else {
	$jnc_str = substr($a->jnc_str,0,2).substr($b->jnc_str,2,2);
	my $qInsert_a = $a->qInsert_length;
	my $qInsert_b = $b->qInsert_length;
	$qInsert_length = $qInsert_a < $qInsert_b ? $qInsert_a : $qInsert_b;
    }
    return $a->new(is_broken => $a->is_broken,
		   jnc_str => $jnc_str,
		   qInsert_length => $qInsert_length);
}

1;

