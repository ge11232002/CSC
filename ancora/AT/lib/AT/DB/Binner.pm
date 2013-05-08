# AT::DB::Binner module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::DB::Binner - handle calculations required for use of UCSC's binning
scheme

=head1 SYNOPSIS

 # Make query to retrieve all features that overlap with
 # the range [$start, $end] on $chr

 use AT::DB::Binner;
 my ($chr, $start, $end) = ('chr11', 1001000, 1005000);
 my $binner = AT::DB::Binner->new;
 my $bin_string = $binner->bin_restriction_string($start,$end);
 my $query = qq!
    SELECT *
    FROM my_table
    WHERE
     $bin_string
     AND tName = "$chr"
     AND tStart <= $end
     AND tEnd >= $start
    !;         

=head1 DESCRIPTION

This module implements a binning scheme that makes it possible to
efficiently retrieve range-type features (e.g. mappings) that overlap
with a range of choice. The code is ported from Jim Kent's binRange.c. 
For further details about the binning scheme, see the UCSC Genome
Browser paper. 

This module is inherited by AT::DB::GenomeMapping, so if you use
AT::DB::GenomeMapping, you don't have to worry about binning for mapping
retrieval.

Note: for binning work efficiently, the bin field of the table in
question has to be indexed.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut



package AT::DB::Binner;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;

use constant BIN_LEVEL_OFFSETS => (512+64+8+1, 64+8+1, 8+1, 1, 0);

@ISA = qw(AT::Root);

sub bin_first_shift { return 17; }

sub bin_next_shift { return 3; }

sub bin_levels { return 5; }

sub bin_level_offset
{
    my ($self, $level) = @_;
    #my @BIN_LEVEL_OFFSETS => (512+64+8+1, 64+8+1, 8+1, 1, 0);

    my $offset = (&BIN_LEVEL_OFFSETS)[$level]; # <-- right way to treat const?
    croak "Bin offset undefined for level $level" unless (defined $offset);
    return $offset;
}


=head2 new

 Title     : new
 Usage     : my $binner = AT::DB::Binner->new();
 Function  : Constructor
 Returns   : AT::DB::Binner
 Args      : -
 
=cut

sub new
{
    my ($caller, %args) = @_;
    my $self = bless {}, ref $caller || $caller;
    return $self;
}


=head2 bin_restriction_string

 Title     : bin_restriction_string
 Usage     : my $bin_string = $binner->bin_restriction_string
		(40000, 50000);
 Function  : Given a coordinate range, return a string to be
             used in the WHERE section of a SQL SELECT statement
	     that is to select features overlapping a certain range.
	     * USE THIS WHEN QUERYING A DB *
 Returns   : String to be used as part of SQL statement
 Args      : 1. Start of coord. range. Required.
             2. End of coord. range. Required.
             3. Name of bin column. Optional. Default: bin.

=cut

sub bin_restriction_string
{
    my ($self, $start, $end, $field) = @_;
    $field = 'bin' unless($field);
    my @bin_ranges = $self->bin_ranges_from_coord_range($start, $end);
    my $str = '('.
		(join " or ", map
	        { ($_->[0] == $_->[1]) ?
	          $field.' = '.$_->[0] :
		  '('.$field.' >= '.$_->[0].' and '.$field.' <= '.$_->[1].')' }
		@bin_ranges)
	      .')';
    return $str;
}


=head2 bin_ranges_from_coord_range

 Title     : bin_ranges_from_coord_range
 Usage     : my @ranges = $binner->bin_ranges_from_coord_range
		(40000, 50000);
 Function  : Return the set of bin ranges that overlap a given
	     coordinate range. It is usually more convenient
	     to use bin_restriction string than to use this
	     method directly.
 Returns   : 2D array (list of start_bin, end_bin pairs)
 Args      : start of coord. range, end of coord. range

=cut

sub bin_ranges_from_coord_range
{
    my ($self, $start, $end) = @_;

    my $start_bin = ($start-1) >> $self->bin_first_shift;
    my $end_bin = ($end-1) >> $self->bin_first_shift;

    my @bin_ranges;

    for(my $i = 0; $i < $self->bin_levels; $i++) {
	my $offset = $self->bin_level_offset($i);
	$bin_ranges[$i] = [ $offset + $start_bin, $offset + $end_bin ];
	$start_bin >>= $self->bin_next_shift;
	$end_bin >>= $self->bin_next_shift;
    }

    return @bin_ranges;
}

=head2 bin_from_coord_range

 Title     : bin_from_coord_range
 Usage     : my $bin = $binner->bin_from_coord_range
		(40000, 50000);
 Function  : Return the bin number that should be assigned to
	     a feature spanning the given range.
	     * USE THIS WHEN CREATING A DB *
 Returns   : bin number
 Args      : start of coord. range, end of coord. range

=cut

sub bin_from_coord_range
{
    my ($self, $start, $end) = @_;
    my $start_bin = ($start-1) >> $self->bin_first_shift;
    my $end_bin = ($end-1) >> $self->bin_first_shift;
    for (my $i = 0; $i < $self->bin_levels; $i++) {
	if ($start_bin == $end_bin) {
	    return $self->bin_level_offset($i) + $start_bin;
	}
	$start_bin >>= $self->bin_next_shift;
	$end_bin >>= $self->bin_next_shift;
    }
    croak "Coords $start, $end out of range for Binner (largest bin is 512M)";
}
