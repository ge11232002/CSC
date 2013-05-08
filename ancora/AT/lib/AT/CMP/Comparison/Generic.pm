# AT::CMP::Comparison::Generic module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::CMP::Comparison::Generic - comparison between two or more genomic
predictions

=head1 SYNOPSIS

See AT::CMP::ComparisonFactory::GenomicAln for a usage example!

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::CMP::Comparison::Generic;

use strict;
use vars '@ISA';
use Carp;
use AT::CMP::CorrespCoords;
use AT::Tools::SeqHandler;
use AT::Root;


@ISA = qw(AT::Root);


=head2 new

 Title     : new
 Usage     : 
 Function  : Constructor
	     Do not call this method directly. Create object using
	     a ComparisonFactory.
 Returns   : AT::CMP::Comparison::Generic
 Args      : 

=cut

sub new  {
    my ($caller, %args) = @_;
    my $self = bless {
	_compared_objL => ($args{compared_obj_list}
			  or croak 'No compared_obj_list argument'),
	_genomic_seqL => ($args{genomic_seq_list}
			  or croak 'No genomic_seq_list argument'),
	_obj_gseq_map => ($args{obj_gseq_map} || []
			  or croak 'No obj_gseq_map argument'),
	genomic_aln => ($args{genomic_aln}
			  or croak 'No genomic_aln argument'),
	genomic_aln_strand => ($args{genomic_aln_strand} || 1),
    }, ref $caller || $caller;
    return $self;
}


=head2 nr_objects

 Title     : nr_objects
 Usage     : print "We have compared ", $cmp->nr_objects,
		   " objects.\n";
 Function  : Returns number of objects in comparison.
 Returns   : Scalar
 Args      : -

=cut

sub nr_objects
{
    my ($self) = @_;
    return scalar @{$self->{'_compared_objL'}};
}


=head2 get_corresp_coords

 Title     : get_corresp_coords
 Usage     : my $cc = $cmp->get_corresp_coords(ref_pos => 2,
                                               rel_start => 86,
                                               rel_end => 121);
 Function  : Finds the range in each compared object that
	     corresponds to a given range in one of the objects.
 Returns   : AT::CMP::CorrespCoords
 Args      : Supply either ref_obj or ref_pos:
	     ref_obj	Object to use as reference
	     ref_pos	Position in comparison of the object
			to use as reference (value between 1 and
			number of objects)
	    
	     And give the range using either
	     abs_start and abs_end  (chromosomal coords)
	     or
	     rel_start and rel_end  (coords starting at 1 and
		increasing in the direction given by the strand
		of the object)
		
=cut

sub get_corresp_coords
{
    my ($self, %args) = @_;

    # find reference object from args
    my $ref_obj;
    unless($ref_obj = $args{'ref_obj'}) {
	if(my $ref_nr = $args{'ref_pos'})
	{
	    if($ref_nr < 1 or $ref_nr > $self->nr_objects)
	    { croak 'argument ref_pos out of range'; }
	    $ref_obj = $self->{'_compared_objL'}->[$ref_nr-1];
	}
	else { croak 'No ref_obj or ref_pos argument'; }
    }

    # find start and end from args
    my($start, $end);
    if($args{'abs_start'} and $args{'abs_end'}) {
	$start = $args{'abs_start'};
	$end = $args{'abs_end'};
    }
    elsif($args{'rel_start'} and $args{'rel_end'}) {
	($start, $end) = $ref_obj->rel2abs($args{'rel_start'},
					   $args{'rel_end'});
    }
    else {
	croak 'No abs_start/abs_end or rel_start/rel_end argument pair';
    }
    if ($start > $end) {
	croak "Start ($start) must be less than or equal to end ($end)";
    }
    #print STDERR "In: $start $end\n";

    # find relative start & end in corresponding aligned sequence
    my $gseq_idx = $self->_genomic_seq_idx($ref_obj);
    my ($gseq_start, $gseq_end) = AT::Tools::SeqHandler->abs2rel
	($self->{'_genomic_seqL'}->[$gseq_idx], $start, $end);
    #print STDERR "Relative in aln seq: $gseq_start $gseq_end\n";

    #$self->print_genomic_aln(\*STDERR);

    # Reverse coords if seq is 2nd in aln and aln is +/-
    # (Could generalize this so that we have a strand for each gseq)
#    if($ref_nr == 2 and $self->genomic_aln_strand == -1) {
#	($start, $end) = ($ref_obj->length-$end+1, $ref_obj->length-$start+1);
#    }

    # get corresponding relative coords for all aligned seqs
    my $gseq_coords_list =
	$self->_get_corresp_coords($gseq_idx+1, $gseq_start, $gseq_end);
    return undef unless ($gseq_coords_list);
    
    # convert to absolute coords
    for(my $i = 0; $i < @$gseq_coords_list; $i++) {
	my $coords = $gseq_coords_list->[$i];
	@$coords = AT::Tools::SeqHandler->rel2abs
	    ($self->{'_genomic_seqL'}->[$i], @$coords);
    }

    # calculate relative coords for each object
    my @object_coords_list;
    foreach my $object ($self->compared_obj_list) {
	my ($corr_abs_start, $corr_abs_end) =
	    @{$gseq_coords_list->[$self->_genomic_seq_idx($object)]};
	if($corr_abs_start <= $object->end and
	   $corr_abs_end >= $object->start) {
	    # overlap with object
	    $corr_abs_start = $object->start
		if ($object->start > $corr_abs_start);
	    $corr_abs_end = $object->end
		if ($object->end < $corr_abs_end);
	}
	else {
	    # no overlap with object
	    undef $corr_abs_start;
	    undef $corr_abs_end;
	}
	push @object_coords_list,
	    [ $object, $object->abs2rel($corr_abs_start, $corr_abs_end) ];

    }

    # Reverse relative coords of 2nd aligned seq if aln is +/-
#    if($self->genomic_aln_strand == -1) {
#	my $length = $self->{'_compared_objL'}->[1]->length;
#	my $coord_pair = $coords_list->[1];
#	@$coord_pair[1,2] =
#	    ($length-$coord_pair->[2]+1, $length-$coord_pair->[1]+1);
#    }
    
    # create CorrespCoords object and return it
    my $cc = AT::CMP::CorrespCoords->new(coords => \@object_coords_list);
    return $cc;
}


#sub _genomic_seq_rel_coords
#{
#    my ($self, $obj, $start, $end);
#    my $gseq = $self->{'_genomic_seqL'}->[$self->_genomic_seq_idx($obj)];
#    return $gseq->start - 
#}


sub _genomic_seq_idx
{
    my ($self, $obj) = @_;
    my $idx = $self->{'_obj_gseq_map'}->{$obj};
    unless (defined $idx) { croak "Object $obj is not in comparison"; }
    return $idx;
}


# Note: this method only works for two objects for now
sub _get_corresp_coords
{
    my ($self, $ref_nr, $start, $end) = @_;

    # find reference seq and its number from args
    unless($ref_nr == 1 or $ref_nr == 2)
    { croak 'argument ref_nr out of range'; }
    my $ref_seq = $self->{'_genomic_seqL'}->[$ref_nr-1];
    my $other_nr = 3 - $ref_nr;

    # find start and end from args
    if ($start > $end) {
	croak "Start ($start) must be less than or equal to end ($end)";
    }

    #print STDERR "_get_corresp_coords called with ",$ref_nr, " ",$start, "-",
	#$end, "\n";

    # get columns in alignment
    my $a = $self->genomic_aln;

    # get positions in alignment for the other object
    my ($corr_start_loc, $corr_end_loc);
    if($a->isa('AT::Alignment::Composite')) {
	my $start_col =
	    $a->column_from_residue_number($ref_nr, $start);
	my $end_col =
	    $a->column_from_residue_number($ref_nr, $end);
	$corr_start_loc =
	    $a->location_from_column($other_nr, $start_col);
	$corr_end_loc =
	    $a->location_from_column($other_nr, $end_col);
    }
    elsif($a->isa('Bio::SimpleAlign')) {
	my $start_col =
	    $a->column_from_residue_number($ref_seq->id, $start);
	my $end_col =
	    $a->column_from_residue_number($ref_seq->id, $end);
	my $other_gseq = $a->get_seq_by_pos($other_nr);
	$corr_start_loc = $other_gseq->location_from_column($start_col)
	    or return undef;
	$corr_end_loc = $other_gseq->location_from_column($end_col)
	    or return undef;
    }
    else { croak 'Unsupported alignment object type'; }

    # use narrowest range
    my $corr_start_coord = $corr_start_loc->end;
    my $corr_end_coord = $corr_end_loc->start;

    # check ambiguity of returned positions (i.e. whether we landed in a gap)
    # if ambiguity: recurse with narrowest range
    #print STDERR "\tCorr start-end: ", $corr_start_coord,'-',$corr_end_coord;
    if ($corr_start_loc->location_type ne 'EXACT' or
	$corr_start_loc->length != 1 or
	$corr_end_loc->location_type ne 'EXACT' or
	$corr_end_loc->length != 1) {
        # there is some ambiguity
	if($corr_start_coord > $corr_end_coord or
	  ($corr_start_coord == $corr_end_coord and
	   $corr_start_loc->location_type ne 'EXACT' and
	   $corr_end_loc->location_type ne 'EXACT')) {
           # start and end are in the same gap
	    #print STDERR "\treturning undef\n";
	    return undef;
	}
	# recurse
	#print STDERR "\trecursing\n";
	return $self->_get_corresp_coords($other_nr,
					  $corr_start_coord,
					  $corr_end_coord);
    }

    # return the coords we found in a 2D array
    #print STDERR "\thappy!\n";
    my @coords =
	( [ $start, $end ],
          [ $corr_start_coord, $corr_end_coord ] );
    if($other_nr == 1) { @coords = reverse @coords };

    return \@coords;
}


sub print_genomic_aln
{
    my ($self, $stream) = @_;
    $stream = \*STDOUT unless (defined $stream);

    my $a = $self->genomic_aln;
    unless($a->isa('AT::Alignment::Composite')) {
	warn "Unsupported alignment class";
	return;
    }

    print $stream "Genomic alignment for ",
	(join ', ', (map {$_->id_loc_str} $self->compared_obj_list)),
	' [',$a->get_strand,' ',$a->get_score,"]\n";
    my @sa = @{$a->get_subalignments};
    @sa = reverse @sa		    # order @sa according to 1st obj
	if($self->{_compared_objL}->[0]->strand *
	   $self->{_genomic_seqL}->
	   [$self->_genomic_seq_idx($self->{_compared_objL}->[0])]->strand
	   == -1);
    foreach my $aln (@sa) {
	my $alnobj = $aln->get_alignment_obj();
	my @row;
	foreach my $obj ($self->compared_obj_list) {
	    my $i = $self->_genomic_seq_idx($obj);
	    my ($s, $e);
	    if($i == 0)		{ ($s, $e) = @{$aln->get_pos_seq1}; }
	    elsif($i == 1)	{ ($s, $e) = @{$aln->get_pos_seq2}; }
	    else		{ croak "invalid genomic seq index value"; }
	    my ($abs_s, $abs_e) = AT::Tools::SeqHandler->rel2abs
				($self->{_genomic_seqL}->[$i], $s, $e);
	    my ($rel_s, $rel_e) = $obj->abs2rel($abs_s, $abs_e);
	    push @row, "A $rel_s-$rel_e\t$abs_s-$abs_e\t(".($e-$s+1).")";
	}
	print $stream join ("\t", @row), "\n";
    }
}


#sub _aligned_seq_for_object
#{
#    my ($self, $obj) = @_;
#    return $self->{_obj_gseq_map}->{$obj};
#}

1;

