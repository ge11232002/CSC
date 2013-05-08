# AT::DB::MultizAlignment module
#
# Copyright Boris Lenhard and Par Engstrom and 
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::DB::MultizAlignment - interface to flat-file database of
whole-genome multiple alignments

=head1 SYNOPSIS

 use AT::DB::MultizAlignment;

 my $hs_alndb = AT::DB::MultizAlignment->connect(
    db_name => "/path/to/multiz/dir"
    );

 
 my @alignments = $hs_alndb->get_alignments_for_region(
    chr => 'chr22',             
    start => 14472596,
    end => 14473215,
    confine => 1);

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::DB::MultizAlignment;

use lib '/home/boris/DEVEL/AT/lib';
use strict;
use warnings;
use vars '@ISA';
use Bio::Seq;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use DBI;
use Carp;
use AT::Root;

@ISA = qw(AT::Root);

sub connect {
	my ($caller, %args) = @_;
	my $self = bless {}, ref $caller ||$caller;
	if (-e $args{db_name}) {
		$self->{db_name} = $args{db_name};
	}
	else {
		croak("Database ".$args{'db_name'}."does not exist or not a directory");
	}
	return $self;
}

sub get_alignments_for_region {
	my ($self, %args) = @_;
	my @alignments;
	my $chr =  $args{chr};
	my $start = $args{start};
	my $end = $args{end};
	$self->{fh} = $self->_get_fh_for_chr($chr);
	my $record_pointer = $self->_move_fh_to_start($start, $end);
	if (defined $record_pointer) {
		@alignments = $self->_read_alignments($start, $end);
		return @alignments;
	}
	else {
		return undef;	
	}
}

sub get_pairwise_alignment_hash_for_region {
	my ($self, %args) = @_;
	my @alignments = $self->get_alignments_for_region(%args);
	my %aln_hash;
	foreach my $aln (@alignments) {
		my @seqs = $aln->each_seq;
		my $ref_seq = shift @seqs;
		foreach my $seq (@seqs) {
			my ($assembly, $chr) = split ".", $seq->display_id;
			my $pairwise_aln = Bio::SimpleAlign->new ();
			$pairwise_aln->add_seq($ref_seq);
			$pairwise_aln->add_seq($seq);
			push @{$aln_hash{$assembly}}, $pairwise_aln;
		}
	}
	return %aln_hash;
}

sub get_alignments_of_seqs  {
	my ($self, @args) = @_;
}

sub _get_fh_for_chr {
	my ($self, $chr) = @_;
	open FH, $self->db_name. "/$chr.maf"
		or croak("Could not open file $chr.maf in database dir ".$self->db_name);
	return \*FH;
}

sub _move_fh_to_start {
	my ($self, $start, $end, $confine) = @_;
	my $filesize = -s $self->fh;
	#print STDERR "FILESIZE: $filesize\n";
	my $step = $filesize;
	my $direction = 1;
	my ($s, $e, $record_pointer);
	my %pos = (start=>$s, end=>$e, rp =>$record_pointer);
	my ($previous_start, $previous_direction, $previous_rp, $converged, 
		$closest_upstream_start, $closest_upstream_pointer, %start_count);
	$previous_rp = tell $self->fh;
	$closest_upstream_pointer = 0;
	$previous_direction = 1;
	while (!$converged) {
		($s, $e, $record_pointer) = $self->_bounds_of_next_alignment;
		return undef if !defined $s;
		#print STDERR "D,S,E,step, RP, PRP, CUP: $direction $s $e $step $record_pointer $previous_rp $closest_upstream_pointer\n";
		if ($e>$start or $s>=$start) {
			$direction = -1;
		}
		elsif ($e<$start) {
			$direction = 1;			
		}
		if ($s<$start and (!defined($closest_upstream_start) or $s>$closest_upstream_start)) {
			$closest_upstream_start = $s;
			#print STDERR "CUS, CUP: $closest_upstream_start, '$previous_rp'\n";
			$closest_upstream_pointer = ($record_pointer or 0);
		}
		$converged = ($step<=1 and $s == $previous_start and $s<=$start);
		if ($s<=$start and $e>=$end) { $converged=1; }
		$start_count{$s}++;
		if ($s>$start and $start_count{$s}>4){
			$converged = 1;
			$record_pointer = $closest_upstream_pointer;
		}
		if (!$converged) {
			if ($s==$start) { $converged = 1 ;}
			if ($previous_rp == 0 or !($direction==-1 and $s==$previous_start and $s>$start)) {
				$step = (int($step/2) or 2);
			}
			$previous_start = $s;
			if ($direction == $previous_direction) {
				seek $self->fh, $previous_rp, 0;
			}	
			$previous_direction = $direction;
			if ($previous_rp +$direction*$step <0) {
				seek $self->fh,0,0;
			}
			else {
				seek $self->fh, $direction * $step, 1;
			}
			$previous_rp = tell $self->fh;
		}
			
		
	}
	
	seek $self->fh, $closest_upstream_pointer, 0;
	return $closest_upstream_pointer;
}

sub _bounds_of_next_alignment {
	my ($self, %args) = @_;
	my $fh = $self->fh;
	<$fh>; my $first_whole_line = <$fh>;	
	unless ($first_whole_line =~/^\#/) {
		local $/ = "\n\n";
		<$fh>;
	}
	my $record_pointer = tell $self->fh;
	my $aln = 	$self->next_aln;
	if (defined $aln) {
		my $seq =  $aln->get_seq_by_pos(1);
		return ($seq->start, $seq->end, $record_pointer);
	}
	else {
		return undef;
	}
}

sub _read_alignments {
	my ($self, $start, $end) = @_;
	my @alignments;
	while (my $aln = $self->next_aln){
	 	my $seq = $aln->get_seq_by_pos(1);
	 	#print STDERR "SEQSTART ", $seq->start, " SEQEND: ", $seq->end, "\n";
	 	if ($seq->end < $start) {
	 		next
	 	}
	 	elsif ($seq->start > $end) {
	 		last
	 	}
	 	elsif ($seq->start >= $start and $seq->end <= $end) {
	 		push @alignments, $aln;
	 	}
	 	else {
	 		my ($slice_start, $slice_end);
	 		if ($seq->start >= $start) {
	 			$slice_start = 1;
	 		}
	 		else {
	 			$slice_start = $seq->column_from_residue_number($start);
	 		}
	 		if ($seq->end <= $end) {
	 			$slice_end = $aln->length;
	 		}
	 		else {
	 			$slice_end = $seq->column_from_residue_number($end);
	 		}
	 		my $aln_slice = $aln->slice($slice_start, $slice_end);
			for my $i (1..$aln->no_sequences) {
			    my $orig_seq = $aln->get_seq_by_pos($i);
			    my $orig_id = $orig_seq->id;

			    foreach my $seq_slice($aln_slice->each_seq){
				my $id_slice = $seq_slice->id;
				if($id_slice eq $orig_id){
#				    my $seq_slice = $aln_slice->get_seq_by_pos($i);
				    $seq_slice->accession_number($orig_seq->accession_number);
				    $seq_slice->strand($orig_seq->strand);
				}
			    }
			}
			push @alignments, $aln_slice;
	 	}
	}
	return @alignments;
}

sub next_aln {
	my ($self) = @_;
	my $alnstring;
	{
		local $/ = "\n\n";
		my $fh = $self->fh;
		$alnstring = <$fh>;
	}
	return unless $alnstring;
	my @cols = split "\n", $alnstring;
	while (@cols and $cols[0] !~ /^a\s+score/) { shift @cols; }
	return unless(@cols);
	my ($aline, @seqlines) = @cols;
	my ($score)= $aline =~ /a score=(.+)/;
	my $alnobj = Bio::SimpleAlign->new;
	foreach my $seqline (@seqlines) {
		next if $seqline =~/^\#/;
		my ($s, $seqname, $start, $length, $strand, $cons_pointer, $seq) = split /\s+/, $seqline;
		last if !$s;  next if $s ne "s";
		#last if !$s; # blank line = end of $alnstring
		my ($assembly, $chr) = split /\./, $seqname;
		my $seqobj = Bio::LocatableSeq->new(-id => $assembly,
						    -accession_number => $chr,
						    -seq => $seq,
						    -start=>$start+1,
						    -end =>$start+$length,
						    -strand =>($strand eq "+" ? 1 : -1)
						    );
		$alnobj->add_seq($seqobj);
	}
	return $alnobj;
}


sub DESTROY {
	my $self = shift;
	close $self->fh;
}

1;
