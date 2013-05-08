############################################################################
# Composite_Alignment.pm - perl module for handeling a set of subalignments 


=head1 NAME

CompositeColumn

=head1 SYNOPSIS

    use AT::Alignment::CompositeColumn;
    my $composite_alignment = AT::Alignment::CompositeColumn -> new(subaln  => $subalingment,
								    column  => $column);

=head1 DESCRIPTION

B<AT::Alignment::CompositeColumn> is a module for handeling columns of a composite alignment.

=head1 METHODS DESCRIPTION

=cut

############################################################################

package AT::Alignment::CompositeColumn;


use strict;
#use vars '@ISA';
#use AT::Root;
use Carp;

#@ISA = qw/AT::Root/;

# Conventions:
# Subalignment and column numbering start at 1.
# If column is 0, the position is between the indicated subaln and the
# previous.
# (subaln, column) = (1,0) means that the position is before the first
# subalignment.
# To allow indication of a position after the last subalignment, subalignment
# nubering can exceed the number of subalignments by 1.

sub new  {
    my ($caller, %args) = @_;
    croak 'No subaln arg' unless defined ($args{subaln});
    croak 'No column arg' unless defined ($args{column});
    my $self = bless {
        subaln => $args{subaln},
        column => $args{column},
    }, ref $caller || $caller;
    return $self;
}

# return true if column is within a subalignment
sub is_within_subalignment { $_[0]->{'column'} }


sub DESTROY { }



1;
