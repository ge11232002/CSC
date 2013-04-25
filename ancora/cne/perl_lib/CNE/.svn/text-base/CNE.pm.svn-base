# CNE::CNE
#
# Copyright 2007 Lenhard Group, BCCS, University of Bergen
# You may distribute this module under the same terms as perl itself
#

=head1 NAME

CNE::CNE - conserved noncoding element

=head1 SYNOPSIS

=head1 METHODS

This section details each of the object
methods. Internal methods are usually preceded with a _

=cut

package CNE::CNE;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use AT::Tools::SeqHandler;
use Data::Dumper;
use Bio::SimpleAlign;


@ISA = qw(AT::Root);

=head2 new

 Title     : new
 Usage     :
 Function  : Constructor
 Returns   : CNE::CNE
 Args      : One corresponding to each get/set method below.

=cut

sub new  {
    my ($caller, %args) = @_;
    my $self = bless { db => $args{db},
		       table_name => $args{table_name},
		       id => $args{id},
		       assembly1 => $args{assembly1},
		       assembly2 => $args{assembly2},
		       chr1 => $args{chr1},
		       chr2 => $args{chr2},
		       start1 => $args{start1},
		       start2 => $args{start2},
		       end1 => $args{end1},
		       end2 => $args{end2},
		       strand => $args{strand},
		       similarity => $args{similarity},
		       cigar => $args{cigar}
		       }, ref $caller || $caller;
    #print Dumper($self), "\n";
    croak "Can't create CNE object: required attribute(s) undefined" unless $self->validate;
    return $self;
}


=head2 Get/set methods

    db         CNE database object (CNE::DB)
    table_name Name of database table containing CNE
    id         Database id for this CNE
    assembly1  Name of first assembly (e.g. hg18)
    assembly2  Name of second assembly (e.g. danRer4)
    chr1       Sequence name, 1st assembly
    chr2       Sequence name, 2nd assembly
    start1     Start coord, 1st assembly. Coords are 1-based, inclusive.
    start2
    end1
    end2
    strand     Character representation of alignment strand, + or -
    similarity Percent identity (0-100, can be a decimal number)
    cigar      Cigar string

    If the CNE object was obtained from through a CNE::DB object,
    the cigar string is lazy-loaded (loaded when needed) to save
    memory.

=cut


sub db { my $self = shift; if(my $x = shift) { return $self->{db} = $x; } return $self->{db}; } 
sub table_name { my $self = shift; if(my $x = shift) { return $self->{table_name} = $x; } return $self->{table_name}; } 
sub id { my $self = shift; if(my $x = shift) { return $self->{id} = $x; } return $self->{id}; } 
sub assembly1 { my $self = shift; if(my $x = shift) { return $self->{assembly1} = $x; } return $self->{assembly1}; } 
sub assembly2 { my $self = shift; if(my $x = shift) { return $self->{assembly2} = $x; } return $self->{assembly2}; } 
sub chr1 { my $self = shift; if(my $x = shift) { return $self->{chr1} = $x; } return $self->{chr1}; } 
sub chr2 { my $self = shift; if(my $x = shift) { return $self->{chr2} = $x; } return $self->{chr2}; } 
sub start1 { my $self = shift; if(my $x = shift) { return $self->{start1} = $x; } return $self->{start1}; } 
sub start2 { my $self = shift; if(my $x = shift) { return $self->{start2} = $x; } return $self->{start2}; } 
sub end1 { my $self = shift; if(my $x = shift) { return $self->{end1} = $x; } return $self->{end1}; } 
sub end2 { my $self = shift; if(my $x = shift) { return $self->{end2} = $x; } return $self->{end2}; } 
sub strand { my $self = shift; if(my $x = shift) { return $self->{strand} = $x; } return $self->{strand}; } 
sub similarity { my $self = shift; if(my $x = shift) { return $self->{similarity} = $x; } return $self->{similarity}; } 

sub cigar {
    my $self = shift;
    my $cigar = shift;
    if($cigar) {
	$self->{cigar} = $cigar;
    }
    else {
	$cigar = $self->{cigar};
	unless($cigar) {
	    if($self->db) {
		$self->db->load_cigar_for_cne($self);
		$cigar = $self->{cigar};
	    }
	    else {
		croak "Cannot load cigar for CNE: database undefined";
	    }
	}
    }
    return $cigar;
}


=head2 strand_numeric

 Title     : strand_numeric
 Usage     : $cne->strand_numeric();
 Function  : Get strand of CNE as number
 Returns   : -1 or 1
 Args      : None

=cut

sub strand_numeric {
    my $strand = shift->strand;
    if($strand eq '+') {
	return 1;
    }
    elsif($strand eq '-') {
	return -1;
    }
    else {
	croak "Unrecognized $strand for CNE";
    }
}


=head2 clone

 Title     : clone
 Usage     : $new_cne = $cne->clone
 Function  : Make an identical copy of the object
 Returns   : A new CNE object
 Args      : None

=cut

sub clone
{
    my $self = shift;
    return $self->new(%$self);
}


=head2 validate

 Title     : validate
 Usage     : $cne->validate
 Function  : Check that all required attributes are defined.
 Returns   : 1 (success) or 0 (failure)
 Args      : None

=cut

sub validate
{
    my $self = shift;
    if($self->assembly1 and $self->chr1 and $self->start1 and $self->end1 and 
       $self->assembly2 and $self->chr2 and $self->start2 and $self->end2 and
       $self->strand and $self->similarity and
       (($self->db and $self->table_name) or $self->cigar))  # cigar can be lazy-loaded
    {
	return 1;
    }
    else {
	return 0;
    }
}


=head2 swap_locations

 Title     : swap_locations
 Usage     : $cne->swap_locations
 Function  : Exchange location 1 and 2. This exchanges
             $cne->assembly1 with $cne->assembly2,
             $cne->chr1 with $cne->chr2,
             etc...
             Also flips the cigar string if one has been loaded.
 Returns   : Self
 Args      : None

=cut

sub swap_locations
{
    my $s = shift;
    ($s->{assembly1}, $s->{assembly2}) = ($s->{assembly2}, $s->{assembly1});
    ($s->{chr1}, $s->{chr2}) = ($s->{chr2}, $s->{chr1});
    ($s->{start1}, $s->{start2}) = ($s->{start2}, $s->{start1});
    ($s->{end1}, $s->{end2}) = ($s->{end2}, $s->{end1});
    $s->flip_cigar() if($s->{cigar});
    return $s;
}


=head2 flip_cigar

 Title     : flip_cigar
 Usage     : $cne->flip_cigar()
 Function  : Flip cigar string.
 Returns   : -
 Args      : None

 NOTE: This method should only be called by the CNE object itself or by its
 database handler. Other modules/scripts can use the method swap_locations
 to swap assemblies for a CNE. That method will call flip_cigar() if necessary.

=cut


sub flip_cigar
{
    my $self = shift;
    my $cigar = $self->cigar;
    $cigar =~ tr/DI/ID/;
    $cigar = _reverse_cigar($cigar) if($self->strand eq '-');
    $self->cigar($cigar);
}


=head2 get_alignment

 Title     : get_alignment
 Usage     : my $aln = $cne->get_alignment($asm1, $asm2);
 Function  : Recreates the alignment for the CNE from its cigar
             string and the provided genome sequence database objects.           
 Returns   : A Bio::SimpleAlign object
 Args      : Two AT::DB::GenomeAssembly objects that should
             correspond to $cne->assembly1 and $cne->assembly2

=cut

sub get_alignment
{
    my ($self, $asm1, $asm2) = @_;
    croak "Missing assembly argument(s)" unless($asm1 and $asm2);

    # Get sequences
    my $seq1 = $asm1->get_genome_seq_str(chr => $self->chr1, start => $self->start1, end => $self->end1);
    my $seq2 = $asm2->get_genome_seq_str(chr => $self->chr2, start => $self->start2, end => $self->end2);
    croak "Failed to retrieve genome sequence(s)" unless($seq1 and $seq2);
    $seq2 = AT::Tools::SeqHandler->revcom($seq2) if($self->strand eq '-');

    # Get cigar and split it into parts
    my $cigar = $self->cigar() or croak "Cigar string undefined";
    my @cigar_parts = $cigar =~ /(\d+)([MDI])/g;

    # Produce gapped sequence strings
    my ($aln_str1, $aln_str2) = ('','');
    while (@cigar_parts) {
	my $n = shift @cigar_parts;
	my $type = shift @cigar_parts;
	if($type eq 'M') {
	    $aln_str1 .= substr($seq1, 0, $n, '');
	    $aln_str2 .= substr($seq2, 0, $n, '');
	}
	elsif($type eq 'I') {
	    $aln_str1 .= substr($seq1, 0, $n, '');
	    $aln_str2 .= ('-' x $n);
	}
	elsif($type eq 'D') {
	    $aln_str1 .= ('-' x $n);
	    $aln_str2 .= substr($seq2, 0, $n, '');
	}
    }

    # Create alignment object and return it
    my $aln = Bio::SimpleAlign->new();
    my $aln_seq1 = Bio::LocatableSeq->new(-id => $self->assembly1.'.'.$self->chr1,
					  -start => $self->start1,
					  -end => $self->end1,
					  -strand => 1,
					  -seq => $aln_str1);
    my $aln_seq2 = Bio::LocatableSeq->new(-id => $self->assembly2.'.'.$self->chr2,
					  -start => $self->start2,
					  -end => $self->end2,
					  -strand => $self->strand_numeric,
					  -seq => $aln_str2);
    $aln->add_seq($aln_seq1);
    $aln->add_seq($aln_seq2);
    return $aln;
}


sub _reverse_cigar {
    my $fwd = shift;
    my $rev = join('', reverse $fwd =~ /(\d+[MDI])/g);
    die "Illegal cigar string $fwd" if(length($rev) != length($fwd));
    return $rev;
}


=head1 AUTHOR

Par Engstrom <par.engstrom@bccs.uib.no>

Copyright 2007 Lenhard Group, BCCS, University of Bergen, Norway

This library is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut

1;

