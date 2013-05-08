# AT::Tools::SeqHandler module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::Tools:SeqHandler - various methods for sequence manipulation

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::Tools::SeqHandler;

use strict;
use vars '@ISA';
use AT::Root;
use Carp;

@ISA = qw(AT::Root);


=head2 rel2abs

 Title     : rel2abs
 Usage     : my $abs_coord = AT::Tools::SeqHandler->
		rel2abs($seq, $rel_coord);
	     my ($abs_start, $abs_end) = AT::Tools::SeqHandler->
		rel2abs($seq, $rel_start, $rel_end);
 Function  : Conversion from relative to absolute coords for
	     a Bio::LocatableSeq object.
	     Relative coords are in the interval (1,length)
	     and increase along the strand of the sequence.
	     Absolute coords are in the interval ($seq->start,
	     $seq->end) and increase along the plus strand.
	     As exemplified, the method can be called with
	     one or two coordinates. If two coords are given,
	     they are taken to represent a range. The first
	     coord must be less than or equal to the second.
	     For the returned range, the first coord
	     will always be less than or equal to the second.
 Returns   : One or two scalars
 Args      : Bio::LocatableSeq and one or two scalars

=cut

sub rel2abs
{
    my ($caller, $seq, $rel_start, $rel_end) = @_;
    if($rel_start < 1 or $rel_start > $seq->length or
       (defined($rel_end) and ($rel_end < $rel_start or
       $rel_end > $seq->length))) {
	croak('Invalid relative coord(s) '.$rel_start.
		(defined($rel_end) ? ('-'.$rel_end) : '').
		' (valid range = 1-'.$seq->length.')');
    }
    if($seq->strand == 1) {
	my $abs_start = $seq->start + $rel_start - 1;
	return $abs_start unless (defined $rel_end);
	my $abs_end = $seq->start + $rel_end - 1;
	return ($abs_start, $abs_end);
    }
    else {
	my $abs_end = $seq->end - $rel_start + 1;
	return $abs_end unless (defined $rel_end);
	my $abs_start = $seq->end - $rel_end + 1;
	return ($abs_start, $abs_end);
    }
}


=head2 abs2rel

 Title     : abs2rel
 Usage     : my $rel_coord = AT::Tools::SeqHandler->
		abs2rel($seq, $abs_coord);
	     my ($rel_start, $rel_end) = AT::Tools::SeqHandler->
		abs2rel($seq, $abs_start, $abs_end);
 Function  : Conversion from absolute to relative coords for
	     a Bio::LocatableSeq object.
	     See rel2abs for more details.
 Returns   : One or two scalars
 Args      : Bio::LocatableSeq and one or two scalars

=cut

sub abs2rel
{
    my ($caller, $seq, $abs_start, $abs_end) = @_;
    if($abs_start < $seq->start or $abs_start > $seq->end or
	(defined($abs_end) and ($abs_end < $abs_start or
	$abs_end > $seq->end))) {
	croak('Invalid absolute coord(s) '.$abs_start.
		(defined($abs_end) ? ('-'.$abs_end) : '').
		' (valid range = '.$seq->start.'-'.$seq->end.')');
    }
    if($seq->strand == 1) {
	my $rel_start = $abs_start - $seq->start + 1;
	return $rel_start unless (defined $abs_end);
	my $rel_end = $abs_end - $seq->start + 1;
	return ($rel_start, $rel_end);
    }
    else {
	my $rel_end = $seq->end - $abs_start + 1;
	return $rel_end unless (defined $abs_end);
	my $rel_start = $seq->end - $abs_end + 1;
	return ($rel_start, $rel_end);
    }
}


=head2 revcom

 Title     : revcom
 Usage     : my $str = AT::Tools::SeqHandler->revcom($str);
 Function  : Revcom a string.
 Returns   : String
 Args      : String

=cut

sub revcom {
    my ($caller, $str) = @_;

    $str =~
        tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    $str = reverse $str;
    return $str;
}


