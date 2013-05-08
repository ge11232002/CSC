# AT::Prediction::NaiveMapPartitioner module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::FT::Factory::GFMappingMachine module

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::FT::Factory::GFMappingMachine;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use AT::FT::GFMapping;
use AT::FT::GFMappingHSP;
use AT::FT::GFMappingGap;

@ISA = qw/AT::Root/;


=head2 new

 Title     : new
 Usage     : 
 Function  : Constructor
 Returns   : 
 Args      :
 
=cut


sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
      	verbose => ( $args{verbose} || 0 ),
	trim_internal_exons => ( defined $args{trim_internal_exons} ? $args{trim_internal_exons} : 1),
	trim_external_exons => ( defined $args{trim_external_exons} ? $args{trim_external_exons} : 1),
	min_external_exon_size => ( defined $args{min_external_exon_size} ? $args{min_external_exon_size} : 10 ),
    }, ref $caller || $caller;
    
    return $self;
}


sub run
{
    my ($self, %args) = @_;
    my $in_mappings = $args{mappings} || croak "No mappings arg";
    my $target_seq = $args{target_seq} || croak "No target_seq arg";    
    my $n = @$in_mappings-1;

    my @structures;
    $#structures = $n;
    for my $i (0..$n) {
	$structures[$i] = $self->_mapping_structure($in_mappings->[$i], $target_seq);
    }
    $self->_verbose_msg("GFMappingMachine: created ",scalar(@structures)," structures\n");

    my $jnc_hash = $self->_create_jnc_hash(\@structures, $target_seq);

    if($self->{trim_internal_exons}) {
	for my $i (0..$n) {
	    $self->_trim_internal_exons($structures[$i], $in_mappings->[$i], $jnc_hash);
	}
    }

    my @gf_mappings;
    for my $i (0..$n) {
	if(@{$structures[$i]} == 0) {
	    warn "GFMappingMachine: skipping mapping of ", $in_mappings->[$i]->qName, " due to empty structure\n";
	    next;
	}
	push @gf_mappings, $self->_create_gf_mapping($in_mappings->[$i],
						     $structures[$i],
						     $jnc_hash);
    }
    
    return \@gf_mappings;
}

#1. create structures for all mappings
#2. make jnc hash, key: start_end coords, values jnc_type, jnc_strand, jnc_seq
#3. create all gfmappings

sub _create_gf_mapping
{
    my ($self, $mapping, $HSPs, $jnc) = @_;
    my @out;
    my $best_for_query = 'F';
    if($mapping->status eq 'P') {
      $best_for_query = 'T';
    }
    elsif($mapping->status eq 'A') {
      $best_for_query = $mapping->status_note > 1 ? 'A' : 'T';
    }
    my $gf_mapping = AT::FT::GFMapping->new
	(best_for_query => $best_for_query,
	 primary_mappings => [$mapping],
	 aln_strand => ($mapping->strand eq '+' ? 1 : -1));
    for my $i (0..@$HSPs-1) {
	my $HSP = AT::FT::GFMappingHSP->new
	    (start => $HSPs->[$i]->[0],
             end => $HSPs->[$i]->[1]);
	my $jnc = $i ? 
	    AT::FT::GFMappingGap->new
		(jnc_str => ($jnc->{$HSPs->[$i-1]->[1]+1}).
			    ($jnc->{$HSPs->[$i]->[0]-2}),
		 qInsert_length => $HSPs->[$i-1][2])
	    : undef;
	$gf_mapping->add_right_HSP($HSP, $jnc);
    }
    return $gf_mapping;
}


sub _mapping_structure
{
    my ($self, $mapping, $target_seq) = @_;
    my $trim_internal = $self->{trim_internal_exons};
    my $min_gap = 4;			# ignore target gaps < 4 bp; should check how this works with the merging
    my @structure;
    my @hsps = $mapping->all_HSPs;

    my $i = 0;
    while ($i < @hsps) {
	my $j = $i+1;
	while($j < @hsps and $hsps[$j]->tStart <= $hsps[$j-1]->tEnd + $min_gap) {
	    $j++;
	}
	# make one HSP from HSPs $i..$j-1
	my $qInsert_length = $j < @hsps ? $hsps[$j]->qStart - $hsps[$j-1]->qEnd - 1 : 0;
	my $qBlockSize = $hsps[$j-1]->qEnd - $hsps[$i]->qStart + 1;
	push @structure, [$hsps[$i]->tStart, $hsps[$j-1]->tEnd, $qInsert_length, $qBlockSize];
	$i = $j;
    }

    if($self->trim_external_exons) {
	$self->_trim_external_exons(\@structure, $target_seq, $mapping);
    }
   
    return \@structure;
}


sub _trim_internal_exons
# Trim internal if < 9 or < 20 and not GT-AG on each side
# unless there are exons of size >= 20 within 10 bp on both sides
{
    my ($self, $s, $mapping, $jnc) = @_;

    my $thr1 = 19;
    my $thr2 = 8;
    my $trimmed_any = 0;

    my $nr_exons = @$s;
    return if($nr_exons < 3);

    my @s2;
    for my $i (1..$nr_exons-2) {
	my $exon = $s->[$i];
	my $start = $exon->[0];
	my $end = $exon->[1];
	my $length = $end-$start+1;
	my $trim = 0;
	if($length <= $thr1 and
	   !($start - $s->[$i-1][1] <= 11 and
	     $s->[$i-1][1]-$s->[$i-1][0] > 19 and
	     $s->[$i+1][0] - $end <= 11 and
	     $s->[$i+1][1]-$s->[$i+1][0] > 19))
	{
	    if($length <= $thr2) {
		$trim = 1;
	    }
	    else {
		my $jnc_str =
		    $jnc->{ $s->[$i-1]->[1]+1 } .
		    $jnc->{ $exon->[0]-2 } .
		    $jnc->{ $exon->[1]+1 } .
		    $jnc->{ $s->[$i+1]->[0]-2 };
		$jnc_str = uc $jnc_str;
		warn "jnc: ", $mapping->qName, " ", $exon->[0], "-", $exon->[1], " $jnc_str\n";
		$trim = 1 unless($jnc_str eq 'GTAGGTAG' or $jnc_str eq 'CTACCTAC');
	    }
	}
	if($trim) {
	    @s2 = map { $s->[$_] } (0..$i-1);
	    $s2[-1][2] += $exon->[2] + $exon->[3];
	    $trimmed_any = 1;
	}
	elsif($trimmed_any) {
	    push @s2, $exon;
	}
    }

    if($trimmed_any) {
	#print STDERR "Trimmed ", $mapping->qName, "\n";
	#print STDERR "  ", join(",", map { ($_->[0]).'-'.($_->[1]) } @$s), "\n";
	#print STDERR "  ", join(",", map { ($_->[0]).'-'.($_->[1]) } (@s2,$s->[-1])), "\n";
    }

    @$s = (@s2,$s->[-1]) if ($trimmed_any);
}


sub _trim_external_exons
{
    my ($self, $s, $target_seq, $mapping) = @_;
       
    my $nr_exons = @$s;
    return if($nr_exons == 1);

    my $min_exon_size = $self->min_external_exon_size();
    my $max_stretch_frac = 79;

    my $first_exon_size = $s->[0][1] - $s->[0][0] + 1;
    #my $first_gap_size = $s->[1][0] - $s->[0][1] - 1;
    my $first_exon_stretch_frac = $self->_stretch_fraction($target_seq, $s->[0][0], $s->[0][1]);

    my $last_exon_size = $s->[-1][1] - $s->[-1][0] + 1;
    #my $last_gap_size = $s->[-1][0] - $s->[-2][1] - 1;
    my $last_exon_stretch_frac = $self->_stretch_fraction($target_seq, $s->[-1][0], $s->[-1][1]);

    #print STDERR join("\t", $mapping->qName, $nr_exons, "F", $first_exon_size, $first_exon_stretch_frac), "\n";
    #print STDERR join("\t", $mapping->qName, $nr_exons, "L", $last_exon_size, $last_exon_stretch_frac), "\n";

    # remove first "exon" ?
    if($first_exon_size < $min_exon_size or
       $first_exon_stretch_frac > $max_stretch_frac) {
	shift @$s;
	#print STDERR "<chop first> ", join(",",map{($_->[0]).'-'.($_->[1])}@$s),"\n";
    }

    # remove last "exon" ?
    if($last_exon_size < $min_exon_size or
       $last_exon_stretch_frac > $max_stretch_frac) {
	pop @$s;
	#print STDERR "<chop last> ", join(",",map{($_->[0]).'-'.($_->[1])}@$s),"\n";
    }

}


sub _stretch_fraction
{
# this criterion may be too strict: add check that we don't have consensus splice junctions?
    my ($self, $seq, $start, $end) = @_;
    my $size = $end-$start+1;
    return 0 if($size > 30);
    my $offset = $seq->start - 1;
    my $seqstr = $seq->subseq($start - $offset, $end - $offset);
    my $nr_Ts = $seqstr =~ tr/Tt//;
    my $nr_As = $seqstr =~ tr/Aa//;
    my $count = $nr_Ts > $nr_As ? $nr_Ts : $nr_As;
    return int(100*$count/$size+.5);
}


sub _create_jnc_hash
{
    my ($self, $structures, $target_seq) = @_;
    my %jnc;
    foreach my $s (@$structures) {
	for my $i (0..@$s-2) {
	    foreach my $coord ($s->[$i]->[1]+1, $s->[$i+1]->[0]-2) {
		unless($jnc{$coord}) {
			$jnc{$coord} =
			    $target_seq->subseq($coord -$target_seq->start +1,
    		    		                $coord -$target_seq->start +2);
		}
	    }
   	}
    }
    return \%jnc;
}

1;
