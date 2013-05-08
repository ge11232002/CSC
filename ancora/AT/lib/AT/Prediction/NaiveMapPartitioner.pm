# AT::Prediction::NaiveMapPartitioner module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::Prediction::NaiveMapPartitioner module

=head1 SYNOPSIS

This is class is meant to provide basic map partitioning without
questioning the supporting mappings much.

It is still being developed and thus subject to change.


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::Prediction::NaiveMapPartitioner;

use strict;
use vars '@ISA';
use Carp;
use AT::Prediction::MapPartitionerI;
use AT::Prediction::Gene;
use AT::Prediction::SpliceForm;

use constant MIN_INTRON_LEN_FOR_SS_CHECK => scalar 12;


@ISA = qw(AT::Prediction::MapPartitionerI);


=head2 new

 Title     : new
 Usage     : my $partitioner = AT::Prediction::NaiveMapPartitioner->new();
 Function  : Constructor
 Returns   : AT::Prediction::NaiveMapPartitioner
 Args      : exclute_on_revcom_ss  If true, excludes mappings with
	                           more revcom'd canonical than
				   canonical splice sites.

=cut


sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
	_exclude_on_revcom_ss => ($args{'exclude_on_revcom_ss'} || 0)
	}, ref $caller || $caller;
    
    return $self;
}


=head2 partition

 Title     : partition
 Usage     : my @genes = $partitioner->partition(mappings => \@mappings,
			             	         target_seq => $seq);
 Function  : 
 Returns   : An array of AT::Prediction::Gene objects
 Args      : mappings    Ref to array of AT::Mapping objects
             target_seq	 A sequence spanning all the mappings
			 (Bio::LocatableSeq-compliant)

=cut

sub partition
{
    my ($self, %args) = @_;
    my $mappings = $args{mappings} || croak "No mappings arg";
    my $target_seq = $args{target_seq} || croak "No target_seq arg";
    my $excluded_mappings = [];

    if($self->{'_exclude_on_revcom_ss'}) {
	($mappings, $excluded_mappings) =
	    $self->_part_on_revcom_ss($mappings, $target_seq);
    }

    my $parts = $self->_partition_into_genes($mappings);

    if(0) {  # Debug output
	print "\nPartitions:\n";
	foreach my $gene (@$parts) {
	    foreach my $m (@$gene) {
		print $m->qName;
		foreach my $hsp ($m->all_HSPs) {
		    print "\t", $hsp->tStart, "-", $hsp->tEnd;
		}
		print "\n";
	    }
	    print "***\n";
	}
	print "\n";
    }

    my @genes;
    foreach my $part (@$parts) {
	# find limits of part and get subseq that spans part
	my ($start, $end);
	foreach my $m (@$part) {
	    $start = $m->tStart if(!defined($start) || $start > $m->tStart);
	    $end = $m->tEnd if(!defined($end) || $end < $m->tEnd);
	}
	my $seqstr = $target_seq->subseq($start - $target_seq->start + 1,
					 $end - $target_seq->start + 1);
	my $seq = Bio::LocatableSeq->new(-id => $target_seq->id,
					 -start => $start,
					 -end => $end,
					 -strand => $target_seq->strand,
					 -seq => $seqstr);
	# create hsp_clusters
	my $hsp_clusters = $self->_create_hsp_clusters($part);
	# create objects
	my $sf = AT::Prediction::SpliceForm->new(mapping_list => $part);
       	my $gene = AT::Prediction::Gene->new
	    (spliceform_list => [$sf],
	    hsp_cluster_list => $hsp_clusters,
	    excluded_mapping_list => $excluded_mappings,
	    seq => $seq);
	push @genes, $gene;
    }

    return @genes;
}


# This method is overly complicated for doing what it does. The reason is that
# it was made by removing stuff from a more advanced version...
sub _partition_into_genes
{
    my ($self, $mappings) = @_;

    my @mtr_queue;              # = matched target regions
    my $min_gap = 12;
    #my @se_list;                # start, end list; should be replaced

    # Make mtr queue
    for(my $i = 0; $i < @$mappings; $i++) {
	my @hsps = $mappings->[$i]->all_HSPs;
	for(my ($j, $k) = (0, 0); $k < @hsps; $k++) {
	    if($k == (@hsps-1) or
	       $hsps[$k+1]->tStart >= $hsps[$k]->tEnd + $min_gap) {
		push @mtr_queue, { mapping => $i,
				   start => $hsps[$j]->tStart,
				   end => $hsps[$k]->tEnd };
		$j = $k+1;
	    }
	}
    }

    # Sort mtr queue
    @mtr_queue = sort {$a->{start} <=> $b->{start}} @mtr_queue;
    
    # Make graph
    my @graph;
    for my $i (0 .. @$mappings-1) { $graph[$i] = [ (0) x @$mappings ] };
    my @mtr_cluster = (shift @mtr_queue);
    my $min_end = $mtr_cluster[0]->{end};
    while(my $mtr = shift @mtr_queue) {
	if($mtr->{start} <= $min_end) {
	    # join each mapping in the cluster with the current mapping
	    foreach my $cl_mtr (@mtr_cluster) {
		$graph[$cl_mtr->{mapping}][$mtr->{mapping}] = 1;
		$graph[$mtr->{mapping}][$cl_mtr->{mapping}] = 1;
	    }
	    $min_end = $mtr->{end} if ($min_end > $mtr->{end});
	    push @mtr_cluster, $mtr;
	}
	else {
	    # Pick out each mtr in the cluster that overlaps with the current
	    my @new_mtr_cluster;
	    $min_end = $mtr->{end};
	    foreach my $cl_mtr (@mtr_cluster)
	    {
		if ($cl_mtr->{end} >= $mtr->{start}) {
		    $min_end = $cl_mtr->{end} if ($min_end > $cl_mtr->{end});
		    push @new_mtr_cluster, $cl_mtr;
		    my $a = $cl_mtr->{mapping};

		    # join the mapping to the current
		    $graph[$cl_mtr->{mapping}][$mtr->{mapping}] = 1;
		    $graph[$mtr->{mapping}][$cl_mtr->{mapping}] = 1;
		}
	    }
	#    unless (@new_mtr_cluster) {
	#	push @se_list, $self->_new_se(@mtr_cluster);
	    #}
      	    @mtr_cluster = (@new_mtr_cluster, $mtr);
	}
    }
    #push @se_list, $self->_new_se(@mtr_cluster);

    # Print graph (DEBUG)
    if(0) {
	print "\t\t", join (" ", (0..@graph-1)) , "\n";
	foreach my $u (0..@graph-1) {
	    print $u, " ", $mappings->[$u]->qName, " [", ,
	    "]\t";
	    foreach my $v (0..@graph-1) {
		print (($u==$v) ? "- " : ($graph[$u][$v]." "));
	    }
	    print "\n";
	}
    }

    # Disentangle graph into genes using a BFS
    my @genes;
    my @visited;
    foreach my $vertex (0..@graph-1) {
	unless($visited[$vertex]) {
	    $visited[$vertex] = 1;
	    my @gene_mappings = ($mappings->[$vertex]);
	    my @Q = ($vertex);
	    while (@Q) {
		my $u = shift @Q;               # get vertex from queue
		foreach my $v (0..@graph-1) {   # visit each neighbour
		    if($graph[$u][$v] and not $visited[$v]) {
			push @gene_mappings, $mappings->[$v];
			$visited[$v] = 1;
			push @Q, $v;
		    }
		}
	    }
	    push @genes, [ @gene_mappings ];
	}
    }

    return \@genes;
}


sub _create_hsp_clusters
{
    my ($self, $mappings) = @_;

    my @clusters;

    my @hsps = sort {$a->tStart <=> $b->tStart}
    (map {$_->all_HSPs} @$mappings);

    my $start = $hsps[0]->tStart;
    my $end = $hsps[0]->tEnd;
    my @cluster_hsps = ($hsps[0]);
    for(my $i = 1; $i < @hsps; $i++) {
	if($hsps[$i]->tStart > $end+1) {
	    push @clusters, {
		hsps => [@cluster_hsps],
		start => $start,
		end => $end };
	    @cluster_hsps = ($hsps[$i]);
	    $start = $hsps[$i]->tStart;
	    $end = $hsps[$i]->tEnd;
	}
	else {
	    push @cluster_hsps, $hsps[$i];
	    $end = $hsps[$i]->tEnd if ($end < $hsps[$i]->tEnd);
	}
    }
    push @clusters, {
	hsps => [@cluster_hsps],
	start => $start,
	end => $end };

    return \@clusters;
}


sub _part_on_revcom_ss
{
    my ($self, $mappings, $seq) = @_;
    my (@ok, @rev);
    foreach my $mapping (@$mappings) {
	my @ssCount = (0) x 7;
	my @intron_coords = $mapping->tInsert_coords
			(min_tInsert_length => MIN_INTRON_LEN_FOR_SS_CHECK);
	# note: check that hsps are in genomic order!
	foreach my $intron_pos (@intron_coords) {
	    my ($start, $end) = @$intron_pos;
	    #print STDERR "$start\t$end\n";
	    $start -= $seq->start - 1;
	    $end -= $seq->start - 1;
	    #print STDERR "$start\t$end\n";
	    my $ss1 = $seq->subseq($start, $start + 1);
	    my $ss2 = $seq->subseq($end - 1, $end);
	    #print STDERR "$ss1 $ss2\n";
	    if($mapping->strand eq '-') {
		($ss1, $ss2) = ($self->_revcom($ss2), $self->_revcom($ss1));
	    }	
	    $ssCount[$self->_ss_type($ss1, $ss2)]++;
	}
	my $ok_count = $ssCount[1]+$ssCount[2]+$ssCount[3];
	my $rev_count = $ssCount[4]+$ssCount[5]+$ssCount[6];
	if($rev_count > $ok_count) {
	    #print STDERR "Mapping ",$mapping->qName,"/",$mapping->mapping_id,
	    #" (", (join ' ', @ssCount), ") $ok_count $rev_count"; 
	    #print STDERR "\t=> rev\n";
	    push @rev, $mapping;
	}
	else {
	    push @ok, $mapping;
	}
    }
    return (\@ok, \@rev);
}


sub _ss_type
{
    my ($self, $ss1, $ss2) = @_;
    my $both = lc($ss1.$ss2);
    return 1 if ($both eq "gtag");
    return 2 if ($both eq "gcag");
    return 3 if ($both eq "atac");
    return 4 if ($both eq "ctac");
    return 5 if ($both eq "ctgc");
    return 6 if ($both eq "gtat");
    return 0;
}


# this should be put in at::root or some lib package
sub _revcom {
    my ($self, $str) = @_;
    
    $str =~ 
	tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    $str = reverse $str;
    return $str;
}

1;
