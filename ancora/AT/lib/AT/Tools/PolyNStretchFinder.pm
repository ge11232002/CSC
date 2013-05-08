package AT::Tools::PolyNStretchFinder;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;

@ISA = qw(AT::Root);


sub _get_args {
    my ($self, %args1) = @_;
    my %args2 = (
	win_size => ($args1{win_size} || $self->{win_size} || 14),
	min_content => ($args1{min_content} || $self->{min_content} || 10)
    );
    return %args2;
}


sub find_stretch_in_seq {
    my ($self, $seq, $nt, %args) = @_;
    %args = $self->_get_args(%args);
    $seq = $seq->seq if(ref($seq));
    return $self->_search_stretch($seq, $nt, %args);
}


sub _search_stretch
{
    my ($self, $seq, $nt, %args) = @_;

    my $win_size = $args{win_size};
    my $min_count = $args{min_content};

    $nt = lc $nt;
    $seq = lc $seq;
    my @seq = (split //, $seq);
    return 0 if(@seq < $min_count);

    # Loop variables
    my @cumul;
    $#cumul = @seq;
    $cumul[0] = 0;  # score base
    my $i = 1;

    # Step thru the first window before we start checking the counts
    for(; $i <= $win_size and $i <= @seq; $i++) {
	my $nt_score = $seq[$i-1] eq $nt ? 1 : 0;
	$cumul[$i] = $cumul[$i-1] + $nt_score;
    }
    # Return true if found at least $min_count in the first window
    return 1 if($cumul[$i-1] >= $min_count);
    # Otherwise keep stepping
    for(; $i <= @seq; $i++) {
	my $nt_score = $seq[$i-1] eq $nt ? 1 : 0;
	$cumul[$i] = $cumul[$i-1] + $nt_score;
	my $win_score = $cumul[$i] - $cumul[$i-$win_size];
	return 1 if($win_score >= $min_count);
    }

    return 0;
}
