# AT::GFX::SeqColorIndex
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::GFX::SeqColorIndex - module to determine color for a sequence

=head1 SYNOPSIS

 my $seq_color_index = AT::GFX::SeqColorIndex->new();
 my @color = $seq_color_index->get_seq_color("chr10");
   # Here the return value is a 3-element (R,G,B) array

=head1 DESCRIPTION

Currently the module only implements the UCSC color schema.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


package AT::GFX::SeqColorIndex;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;


@ISA = qw(AT::Root);


# Color tables

my $UCSC_COLOR_TABLE =
    [[0,0,0],          # default
     [0x99,0x66,0x00], # chr1
     [0x66,0x66,0x00], # chr2
     [0x99,0x99,0x1e], # chr3	       
     [0xcc,0x00,0x00], # ...
     [0xff,0x00,0x00],
     [0xff,0x00,0xcc], # chr6
     [0xff,0xcc,0xcc],    
     [0xff,0x99,0x00], # chr8
     [0xff,0xcc,0x00],
     [0xff,0xff,0x00], # chr10
     [0xcc,0xff,0x00],
     [0x00,0xff,0x00], # chr12
     [0x35,0x80,0x00],
     [0x00,0x00,0xcc],
     [0x66,0x99,0xff],
     [0x99,0xcc,0xff],
     [0x00,0xff,0xff],
     [0xcc,0xff,0xff],
     [0x99,0x00,0xcc],
     [0xcc,0x33,0xff],
     [0xcc,0x99,0xff],
     [0x66,0x66,0x66], # chr22
     [0x99,0x99,0x99], # chrX
     [0xcc,0xcc,0xcc], # chrY
     [0xcc,0xcc,0x99], # chrM
     [0x79,0xcc,0x3d]  # chrUn
     ];


my $DARK_COLOR_TABLE =
    [[0,0,0],          # default
     [0x99,0x66,0x00], # chr1
     [0x66,0x66,0x00], # chr2
     [0x99,0x99,0x1e], # chr3	       
     [0xcc,0x00,0x00], # chr4
     [0xff,0x00,0x00], # chr5
     [0xff,0x00,0xcc], # chr6
     [0xff,0x99,0x99], # chr7   
     [0xff,0x99,0x00], # chr8
     [0xff,0xcc,0x00], # chr9
     [0xcc,0x50,0x00], # chr10
     [0xcc,0xff,0x00], # chr11
     [0x00,0xff,0x00], # chr12
     [0x00,0x70,0x00], # chr13
     [0x00,0x00,0xcc], # chr14
     [0x66,0x99,0xff], # chr15
     [0x99,0xcc,0xff], # chr16
     [0x00,0xff,0xff], # chr17
     [0xcc,0xa0,0xa0], # chr18
     [0x99,0x00,0xcc], # chr19
     [0xcc,0x33,0xff], # chr20
     [0xcc,0x99,0xff], # chr21
     [0x66,0x66,0x66], # chr22
     [0x99,0x99,0x99], # chrX
     [0xcc,0xcc,0xcc], # chrY
     [0xcc,0xcc,0x99], # chrM
     [0x79,0xcc,0x3d]  # chrUn
     ];


## Public methods

sub new
{
    my ($caller, %args) = @_;

    # Here we can allow selection between different color schemes
    my $scheme = $args{'scheme'} || 'UCSC';
    unless($scheme eq 'UCSC' or $scheme eq 'dark') {
	carp "Invalid color scheme $scheme";
	return undef;
    }

    my $self = bless {
	scheme => $scheme
    }, ref $caller || $caller;

    return $self;   
}


sub get_seq_color
{
    my ($self, @args) = @_;

    my $scheme = $self->{scheme};
    if($scheme eq 'UCSC') {
	return $self->_get_ucsc_seq_color(@args);
    }
    elsif($scheme eq 'dark') {
	return $self->_get_dark_seq_color(@args);
    }
    else {
	croak "invalid color scheme $scheme";
    }
}

sub get_legend
{
    my $self = shift;
    my $scheme = $self->{scheme};
    if($scheme eq 'UCSC') {
	return $self->_get_ucsc_legend();
    }
    elsif($scheme eq 'dark') {
	return $self->_get_dark_legend();
    }
    else {
	croak "invalid color scheme $scheme";
    }
}


## Private methods for UCSC color scheme

sub _get_ucsc_seq_color {
    my $self = shift;
    return $self->_get_seq_color_from_ucsc_type_table($UCSC_COLOR_TABLE, @_);
}

sub _get_ucsc_legend
{
    my $self = shift;
    return $self->_get_legend_for_ucsc_type_table($UCSC_COLOR_TABLE, @_);
}


## Private methods for dark color scheme

sub _get_dark_seq_color {
    my $self = shift;
    return $self->_get_seq_color_from_ucsc_type_table($DARK_COLOR_TABLE, @_);
}

sub _get_dark_legend
{
    my $self = shift;
    return $self->_get_legend_for_ucsc_type_table($DARK_COLOR_TABLE, @_);
}

## Private methods for accessing UCSC-type color tables

sub _get_seq_color_from_ucsc_type_table
{
    my ($self, $colors, $seq_name) = @_;

    my $color;

    if($seq_name =~ /^(?:chr|group)(.+)/i) {
	$color = $self->_get_chr_color_from_ucsc_type_table($colors, $1);
    }
    elsif($seq_name =~ /^(?:scaffold|contig|supercont|super)_?(\d+)/i or $seq_name =~ /^(\d+)/) {
	$color = $self->_get_scaffold_color_from_ucsc_type_table($colors, $1);
    }
    else {
	$color = $colors->[0];
    }

    return @$color;
}

sub _get_chr_color_from_ucsc_type_table
{
    my ($self, $colors, $chr_nr) = @_;
    my $i;
    if($chr_nr =~ /^(\d+)/) { $i = $1; }
    elsif($chr_nr =~ /^U/) { $i = 26; }
    elsif($chr_nr =~ /^X/) { $i = 23; }
    elsif($chr_nr =~ /^Y/) { $i = 24; }
    elsif($chr_nr =~ /^M/) { $i = 25; }
    elsif($chr_nr =~ /^V/) { $i = 5; }
    elsif($chr_nr =~ /^IV/) { $i = 4; }
    elsif($chr_nr =~ /^III/) { $i = 3; }
    elsif($chr_nr =~ /^II/) { $i = 2; }
    elsif($chr_nr =~ /^I/) { $i = 1; }
    else { $i = 0 }
    $i = 0 if ($i >= scalar(@$colors));
    return $colors->[$i];
}

sub _get_scaffold_color_from_ucsc_type_table
{
    my ($self, $colors, $scaffold_nr) = @_;
    my $i = $scaffold_nr % (scalar(@$colors)-1); # In the UCSC schema, the last color is reserved for undefined chr
    $i = 0 if($i < 0);
    return $colors->[$i];
}


sub _get_legend_for_ucsc_type_table
{
    my ($self, $colors) = @_;

    my @names = (map({ "chr$_" } (1..22)), 'chrX, chr23','chrY, chr24','chrM, chr25','chrUn', 'Other');
    my @c = @$colors;
    push @c, shift @c;

    return { 
	colors => \@c,
	names => \@names,
	note => "Scaffolds are assigned colors from the same color table based on their number, recycling colors as needed."
    };
}

1;

