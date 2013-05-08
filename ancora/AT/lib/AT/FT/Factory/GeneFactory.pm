# AT::FT::Factory::GeneFactory module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::FT::Factory::GeneFactory - factory object that produces gene predictions
from cDNA-genome mappings

=head1 SYNOPSIS

 # Connect to databases

 my $mapping_db = AT::DB::GenomeMapping->connect
    ( -dbname => "AT_HS_JUN02",
      -dbhost => "localhost",
      -dbuser => "",
      -dbpass => "");
 my $assembly_db = AT::DB::GenomeAssembly->connect
    ( -dbname => "HS_JUN02",
      -dbhost => "localhost",
      -dbuser => "",
      -dbpass => "");

 # Create GeneFactory

 my $gf = AT::FT::Factory::GeneFactory->new
    ( mapping_db => $mapping_db,
      assembly_db => $assembly_db);

 # Create genes for a 400 kb region on chr2

 my @genes = $gf->create_genes_for_region
    (chr => "chr2",
     start => 150000000,
     end =>   150400000,
     strand => -1);

 # Get some mappings and create genes genes for them

 my @mappings = $mapping_db->get_mappings_for_acc("NM_024796");
 push @genes, $gp->create_genes_for_mappings
    (mappings => \@mappings,
     expand => 1);

 # Output locations of predicted genes

 foreach my $gene (@genes) {
     print "Gene at ",$gene->loc_str, "\n";
 }

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::FT::Factory::GeneFactory;

use strict;
use vars '@ISA';
use Carp;
use AT::Mapping;
use AT::FT::Factory::M2GProcess;
use AT::FT::Factory::GeneStream;

@ISA = qw/AT::Root/;

=head2 new

 Title     : new
 Usage     : my $gpredictor = AT::Prediction::GenePredictor->new
		(mapping_db => $mapping_db,
		 assembly_db => $assembly_db,
		 map_partitioner => $map_partitioner,
		 gene_annotators => [$gene_annotator]);
 Function  : Constructor
 Returns   : An AT::Prediction::GenePredictor object
 Args      : mapping_db	 cDNA-genome mapping db to use. Required.
			 (AT::DB::GenomeMapping object)
	     assembly_db  Assembly db to use. Required.
			  (AT::DB::GenomeAssembly object)
	     map_partitioner  Strategy to partition mappings into
			      genes and spliceforms. Required.
			      (AT::Prediction::MapPartitionerI-compliant)
	     gene_annotators  Strategies for annotating
			      features, e.g. exons and introns. Optional.
			      (Ref to array of
			      AT::Prediction::GeneAnnotatorI-compliants)

=cut

sub new
{
    my ($caller, %args) = @_;

    my $self = bless { }, ref $caller || $caller;
    
    $self->{_query_db} = $args{mapping_db} || croak "No mapping_db";
    $self->{_target_db} = $args{assembly_db} || croak "No assembly_db";
    $self->{_seq_db} = $args{transcript_db};
    unless($self->{_seq_db}) { warn "GeneFactory: warning: no seq_db; will not attach query sequences\n"; }
    #$self->{_map_partitioner} = $args{map_partitioner} || 
	#croak "No map_partitioner";
    #$self->{_gene_annotators} = $args{gene_annotators} || [];
    $self->{m2g_process} = AT::FT::Factory::M2GProcess->new();
    $self->{_qTypeL} = $args{qType_list} || ['mRNA'];

    return $self;
}


=head2 create_genes_for_mappings

 Title     : create_genes_for_mappings
 Usage     : my @genes = $gp->create_genes_for_mappings
			    (mappings => \@mappings,
			     expand => 1);
 Function  : Create gene predictions using a given set of mappings.
	     The mappings need not overlap.
 Returns   : An array of AT::Prediction::Gene objects
 Args      : mappings   Reference to an array of AT::Mapping objects
	     expand     If true, mappings that overlap with the supplied
			mappings (on the same strand) will be retrieved
			from the query db and used in creating the gene
			predictions. Still, only gene predictions
			based on the supplied mappings are returned.
	     restrict	If true, only the supplied mappings will be
			used. (The opposite of expand)
	     Either 'expand' or 'restrict' must be true.

=cut


sub create_genes_for_mappings
{
    my ($self, %args) = @_;

    croak "No mappings arg" unless ($args{'mappings'});
    my @mappings = @{$args{'mappings'}};
    my @genes;

    if($args{'restrict'}) {
    # Use only the supplied mappings to create genes
    
	# Sort mappings by 1.assembly, 2.chromosome, 2.strand, 3.start
	@mappings = sort {$a->target_db cmp $b->target_db ||
			$a->tName cmp $b->tName ||
			#$a->strand <=> $b->strand ||
			$a->tStart <=> $b->tStart} @mappings;

	# The loop below does the following until the list @mappings is empty.
	# 1. Take the first mapping off the list
	#    Set the cluster to span that mapping.
	# 2. If there is a next mapping in the list and it overlaps with the cluster,
	#    extend the cluster with that mapping and take it off the list. Repeat
	#    until the next mapping in the list does not overlap with the cluster.	
	#    (Same strand is a requirement for overlap. Adjacency is not counted as
	#    overlap.)
	# 3. Create genes for the cluster of mappings
    
	my $i = 0;
	while ($i < @mappings) {
	    my $cluster_end = $mappings[$i]->tEnd;
	    my $j = $i+1;
	    while($j < @mappings &&
		$mappings[$i]->target_db eq $mappings[$j]->target_db &&
		$mappings[$i]->tName eq $mappings[$j]->tName &&
		#$mappings[$i]->strand == $mappings[$j]->strand &&
		$mappings[$j]->tStart <= $cluster_end)
	    {
		if ($cluster_end < $mappings[$j]->tEnd) {
		    $cluster_end = $mappings[$j]->tEnd;
		}
		$j++;
	    }
	    my $gene_list_ref = $self->_create_genes
		    (#$mappings[$i]->tStart, $cluster_end,
		     undef,
		     $mappings[$i]->tStart, $cluster_end,
		     [@mappings[$i..$j-1]]);
	    push @genes, @{$gene_list_ref};
	    $i = $j;
	}
    }
    elsif($args{'expand'}) {
    # Create genes containing the supplied mappings, but include more if possible

	my @all_genes = $self->create_genes_for_mapped_regions(mappings => \@mappings);

	# the part below should be moved to filter_by_mappings
	foreach my $gene (@all_genes) {
	    my $gene_contains_supplied_mapping;
	    foreach my $m ($gene->mapping_list) {
		for(my $j = 0; $j < @mappings; $j++) {
		    if($mappings[$j]->mapping_id == $m->mapping_id) {
		        # this won't work if we allow several query db's
		        $gene_contains_supplied_mapping = 1;
		    }
		}
	    }
	    push @genes, $gene if ($gene_contains_supplied_mapping);
        }	
    }
    else {
	croak "Need arg expand or restrict";
    }

    return @genes;
}



sub create_genes_for_mapped_regions
{
    my ($self, %args) = @_;

    croak "No mappings arg" unless ($args{'mappings'});
    my @mappings = @{$args{'mappings'}};

    my @genes;
    my %mappings_used;

    for(my $i = 0; $i < @mappings; $i++) {
        next if ($mappings_used{$mappings[$i]->mapping_id});
        my ($gene_list_ref, undef, $mappings_used) = 
    	$self->_create_genes_for_mapping_cluster
	    ($mappings[$i]->tName,
	     undef,
	     $mappings[$i]->tStart,
	     $mappings[$i]->tEnd);
        foreach my $m (@$mappings_used) {
	    $mappings_used{$m->mapping_id} = 1;
	}
        push @genes, @{$gene_list_ref} if($gene_list_ref);
    }

    return @genes;
}



#
# use hash %mapping_usage;
#  keys: mapping_ids of mappings that were used; values: mapping objects
#  the mapping creation can add to a supplied hash or create a new one
# create_genes_for_mappings goes through all mappings until they all have
# been used.
# filtering by region must not be done; we filter by mappings
#
#


=head2 create_genes_for_region

 Title     : create_genes_for_region
 Usage	   : my @genes = $gp->create_genes_for_region
		    (chr => "chr15",
		     start => 1145000,
		     end => 1148000,
		     strand => -1);
 Function  : Create gene predictions for a given genomic region.
	     If a strand argument is given, only that strand is considered.
 Returns   : An array of AT::Prediction::Gene objects
 Args      : chr, start, end - required
	     strand - optional

=cut


sub create_genes_for_region
{
    my ($self, %args) = @_;
    my ($chr, $strand, $start, $end) = $self->_parse_region_args(%args);
    my @gene_list;
    my $gene_list_ref;
    my $cluster_end;
    ($gene_list_ref, $cluster_end) =
	$self->create_genes_for_first_mapping_cluster(chr => $chr,
						      strand => $strand,
						      start => $start,
						      end => $end);
    push @gene_list, @{$gene_list_ref} if ($gene_list_ref);
    while($cluster_end and $cluster_end < $end) {
	($gene_list_ref, $cluster_end) = 
	    $self->create_genes_for_next_mapping_cluster(chr => $chr,
							strand => $strand,
							start => $cluster_end+1,
							end => $end);
	push @gene_list, @{$gene_list_ref} if ($gene_list_ref);
    }
    $gene_list_ref = $self->_filter_genes_by_region(\@gene_list, $strand, $start, $end);
    return wantarray ? @$gene_list_ref : $gene_list_ref;
}


sub create_genes_for_chr_in_stream
{
    my ($self, $chr) = @_;

    # get chr size
    my $tSize = $self->_target_db->get_chr_size($chr);

    # select all unique tStarts
    my $tStarts = $self->_query_db->dbh->selectcol_arrayref(
	"select distinct tStart from MAPPING where tName = '$chr' order by tStart");
    
    return AT::FT::Factory::GeneStream->new(gene_factory => $self,
					    tName => $chr,
					    tStarts => [@$tStarts], # best to copy the data; don't know what dbh could do to the orignal array
					   );#tEnd => $tSize);
}


=head2 create_genes_for_region_in_stream

 Title     : create_genes_for_region_in_stream
 Usage	   : my @genes = $gp->create_genes_for_region_in_stream
		    (chr => "chr15",
		     start => 1,
		     end => 10000000,
		     strand => 1);
 Function  : Successively create gene predictions for a given
	     genomic region. Only one strand is considered.
	     Preferable for processing large regions.
	     * THIS METHOD IS NOT YET FUNCTIONAL. *
 Returns   : An AT::Prediction::GeneStream object
 Args      : chr, start, end - required
	     strand - optional

=cut


sub create_genes_for_region_in_stream
{
    my ($self, %args) = @_;
    my ($chr, $strand, $start, $end) = $self->_parse_region_args(%args);
    my $tStarts = $self->_query_db->dbh->selectcol_arrayref(
	"select distinct tStart
	 from MAPPING
	 where tName = '$chr'
	   and tEnd >= $start
	   and tStart <= $end
	 order by tStart");
    return AT::FT::Factory::GeneStream->new(gene_factory => $self,
					    tName => $chr,
					    tStarts => [@$tStarts], # best to copy the data; don't know what dbh could do to the orignal array
					   );#tEnd => $end);
}


=head2 create_genes_for_first_mapping_cluster

 Title     : create_genes_for_first_mapping_cluster
 Usage	   : my @genes = $gp->create_first_mapping_cluster
		    (chr => "chr15",
		     start => 1,
		     end => 10000000,
		     strand => 1);
 Function  : This method is primarily meant to be used internally and by
	     AT::Prediction::GeneStream.
	     Create gene predictions for the first cluster of
	     overlapping mappings in a given genomic region.
	     If a strand argument is given, only that strand is considered.
 Returns   : If successful:
	     1. A reference to a set of AT::Prediction::Gene objects
	     2. The position just after the end of the mapping cluster
	        (for use as start in the next call) or undef if that
		position is past the end of the given region.
	     If there are no mappings in the region, undef is returned.
 Args      : chr, start, end  - all required
	     strand - optional

=cut


sub create_genes_for_first_mapping_cluster
{
    my ($self, %args) = @_;
    my ($chr, $strand, $start, $end) = $self->_parse_region_args(%args);

    # change this so that we get any mapping overlapping the first base
    # if none, call create_genes_for_next_mapping_cluster
    $self->_verbose_msg("GeneFactory: getting the leftmost mapping over $chr:$start-$end\n");
    my ($m1) = $self->_query_db->get_mappings_in_region
	(chr => $chr, #strand => $strand,
	 start => $start, end => $end,
	 first => 1, types => $self->{_qTypeL});
    return undef unless ($m1);

    return $self->_create_genes_for_mapping_cluster($chr, $strand,
						    $m1->tStart, $m1->tEnd);
}


sub create_genes_for_next_mapping_cluster
{
    my ($self, %args) = @_;
    my ($chr, $strand, $start, $end) = $self->_parse_region_args(%args);

    $self->_verbose_msg("GeneFactory: getting the leftmost mapping starting within $chr:$start-$end\n");
    # this is not as quick as it should be; can we be more rough
    # and just get a chunk of mappings?
    my ($m1) = $self->_query_db->get_first_mapping_starting_in_region
	(chr => $chr, #strand => $strand,
	 start => $start, end => $end,
	 types => $self->{_qTypeL});
    return undef unless ($m1);

    return $self->_create_genes_for_mapping_cluster($chr, $strand,
						    $m1->tStart, $m1->tEnd);
}


sub _create_genes_for_mapping_cluster
{
    my ($self, $chr, $strand, $cluster_start, $cluster_end) = @_;

    # Step 2: let cluster = all mappings that touch the region of $m1 
    #my ($cluster_start, $cluster_end) =	($m1->tStart, $m1->tEnd);
    my ($prev_cluster_start, $prev_cluster_end) =
	($cluster_start, $cluster_end);
    $self->_verbose_msg("GeneFactory: getting all mappings overlapping $chr:$cluster_start-$cluster_end\n");
    my @mappings = $self->_query_db->get_mappings_in_region
	(chr => $chr, #strand => $strand,
	 start => $cluster_start,  end => $cluster_end,
	 types => $self->{_qTypeL},
	 attach_query_seq => $self->{_seq_db});
    foreach my $m (@mappings) {
	$cluster_start = $m->tStart if($m->tStart < $cluster_start);
	$cluster_end = $m->tEnd if($m->tEnd > $cluster_end);
    }

    # Step 3: expand cluster to the right
    while($cluster_end > $prev_cluster_end) {
        $self->_verbose_msg("GeneFactory: getting all mappings starting within $chr:",$prev_cluster_end+1,"-$cluster_end\n");
	my @some_mappings = $self->_query_db->get_mappings_in_region
	    (chr => $chr, #strand => $strand,
	     start => $prev_cluster_end+1, end => $cluster_end,
	     confine_start => 1, types => $self->{_qTypeL},
	     attach_query_seq => $self->{_seq_db});
	push @mappings, @some_mappings;
	$prev_cluster_end = $cluster_end;
	foreach my $m (@some_mappings) {
	    $cluster_end = $m->tEnd if($m->tEnd > $cluster_end); 
	}
    }

    # Step 4: expand cluster to the left
    while($cluster_start < $prev_cluster_start) {
        $self->_verbose_msg("GeneFactory: getting all mappings ending within $chr:$cluster_start-",$prev_cluster_start-1,"\n");
	my @some_mappings = $self->_query_db->get_mappings_in_region
	    (chr => $chr, #strand => $strand,
	     start => $cluster_start, end => $prev_cluster_start-1,
	     confine_end => 1, types => $self->{_qTypeL},
	     attach_query_seq => $self->{_seq_db});
	push @mappings, @some_mappings;
	$prev_cluster_start = $cluster_start;
	foreach my $m (@some_mappings) {
	    $cluster_start = $m->tStart if($m->tStart < $cluster_start); 
	}
    }
    $self->_verbose_msg("GeneFactory: ...got all mappings\n");

    # Create gene using the mappings
    my $genes = $self->_create_genes($strand,
				     $cluster_start, $cluster_end,
				     \@mappings);

    # Return genes and the new start
    return ($genes, $cluster_end, \@mappings);
}



=head2 _create_genes

 Title     : _create_genes
 Usage     : 
 Function  : Private method.
 Returns   : 
 Args      : 

=cut


sub _create_genes
{
    my ($self, $strand_req, $min_start, $max_end, $mappings) = @_;

    # For debug
    if(0) {
	print "Creator: ";
	print $mappings->[0]->tName, "\t";
	print $mappings->[0]->strand, "\t";
	print $min_start, "-", $max_end ,"\t";
	print scalar(@$mappings), "\n";
    }

    # Attach transcript seqs to mappings
    if($self->{_seq_db}) {
	$self->_verbose_msg("GeneFactory: attaching query sequences\n");
	$mappings = [
	    grep {
		if($_->query_seq) {
		    1;
		}
		else {
		    my $seq = $self->{_seq_db}->get_seq($_->qName);
		    if($seq) {
		        $_->attach_query_seq($seq);
		        1;
		    }
		    else {
		        warn "No query seq with acc ".$_->qName;
		        0;
		    }
		}
	    } @$mappings
	];
    }
    return [] unless(@$mappings);

    # Get genomic sequence for the region
    $self->_verbose_msg("GeneFactory: getting genome sequence\n");
    my $tSeq_start = $min_start - 1000;  # usually 100 nt context is enough; instead of getting this much we should make genome sequence expandable (see ideas)
    $tSeq_start = 1 if($tSeq_start < 1);
    my $tSeq_end = $max_end + 1000;
    $tSeq_end = $mappings->[0]->tSize if($tSeq_end > $mappings->[0]->tSize);
    my $target_seq = $self->_target_db->get_genome_seq
	(chr => $mappings->[0]->tName,
	 #strand => $mappings->[0]->strand,
	 start => $tSeq_start,
	 end => $tSeq_end);

    # Create genes
    $self->_verbose_msg("GeneFactory: Running M2Gprocess with ", scalar(@$mappings), " mappings\n");
    my @genes = $self->m2g_process->run(mappings => $mappings,
				        target_seq => $target_seq,
				        strand => $strand_req || undef
					     );

    if(0) {
	foreach my $gene (@genes) {
	    print "gene ",$gene->start,"-",$gene->end, "\t";
	    print join " ", (map { $_->qName } $gene->mapping_list);
	    print "\n";
	    print "\t", join (" ", (map { $_->{start}."-".$_->{end} }
				    $gene->hsp_cluster_list)), "\n";
	}
    }

    return \@genes;
}


sub create_genes_for_gfmappings
{
    my ($self, %args) = @_;

    croak "No gfmappings arg" unless ($args{'gfmappings'});
    my @mappings = @{$args{'gfmappings'}};
    return unless(@mappings);

    my $min_start = $mappings[0]->start;
    my $max_end = $mappings[0]->end;
    foreach my $m (@mappings[1..@mappings-1]) {
	my $start = $m->start;
	my $end = $m->end;
	$min_start = $start if($min_start > $start);
	$max_end = $end if($max_end < $end);
    }
  

    $self->_verbose_msg("GeneFactory: getting genome sequence\n");
    my $chr = $mappings[0]->chr;
    my $tSeq_start = $min_start - 1000;
    $tSeq_start = 1 if($tSeq_start < 1);
    my $tSeq_end = $max_end + 1000;
    my $chr_size = $self->_target_db->get_chr_size($chr);
    $tSeq_end = $chr_size if($tSeq_end > $chr_size);
    my $target_seq = $self->_target_db->get_genome_seq
	(chr => $chr,
	 start => $tSeq_start,
	 end => $tSeq_end);

    # Create genes
    $self->_verbose_msg("GeneFactory: Running M2Gprocess with ", scalar(@mappings), " GFMappings\n");
    my @genes = $self->m2g_process->run_from_gfmappings(gfmappings => \@mappings,
							target_seq => $target_seq);

    return @genes;
}


sub _filter_genes_by_region
{
    my ($self, $genes, $strand, $start_bound, $end_bound) = @_;

    print STDERR "filtering: $start_bound-$end_bound\n";

    # Keep only those genes that overlap with the region
    my @filtered_genes =
	$strand ?
	    (grep { $_->start <= $end_bound and $_->end >= $start_bound and
	            $_->strand == $strand } @$genes)
	    :
	    (grep { $_->start <= $end_bound and $_->end >= $start_bound } @$genes);
    $self->_verbose_msg("GeneFactory: Keeping ",scalar(@filtered_genes),
			" of ",scalar(@$genes)," genes\n");
    return \@filtered_genes;
}


sub _filter_genes_by_mappings
{
    my ($self, $genes, $mappings) = @_;

    # Keep only those genes that overlap with the region
    my @filtered_genes = grep { $_ == 1 } @$genes;
    $self->_verbose_msg("GeneFactory: Keeping ",scalar(@filtered_genes),
			" of ",scalar(@$genes)," genes\n");
    return \@filtered_genes;
}


=head2 _parse_region_args

 Title     : _parse_region_args
 Usage     : 
 Function  : Private method.
 Returns   : 
 Args      : 

=cut


sub _parse_region_args
{
    my ($self, %args) = @_;
    my $chr = $args{chr} || croak "No chr";
    my $strand = $args{strand} || undef;
    my $start = $args{start} || croak "No start";
    my $end = $args{end} || croak "No end";
    return ($chr, $strand, $start, $end);
}


#sub create_genes_for_first_mapping_cluster
#{
#    my ($self, %args) = @_;
#    my ($chr, $strand, $start, $end) = $self->_parse_region_args(%args);
#
#    # Get the first cluster of overlapping mappings in the region
#
#    # Step 1: let $m1 = the leftmost mapping touching the region
#    $self->_verbose_msg("GeneFactory: getting the leftmost mapping over $chr:$start-$end\n");
#    my ($m1) = $self->_query_db->get_mappings_in_region
#	(chr => $chr, #strand => $strand,
#	 start => $start, end => $end,
#	 first => 1, types => $self->{_qTypeL});
#    return undef unless ($m1);
#
#    # Step 2: let cluster = all mappings that touch the region of $m1 
#    my ($cluster_start, $cluster_end) =	($m1->tStart, $m1->tEnd);
#    my ($prev_cluster_start, $prev_cluster_end) =
#	($cluster_start, $cluster_end);
#    $self->_verbose_msg("GeneFactory: getting all mappings overlapping $chr:$cluster_start-$cluster_end\n");
#    my @mappings = $self->_query_db->get_mappings_in_region
#	(chr => $chr, #strand => $strand,
#	 start => $cluster_start,  end => $cluster_end,
#	 types => $self->{_qTypeL});
#    foreach my $m (@mappings) {
#	$cluster_start = $m->tStart if($m->tStart < $cluster_start);
#	$cluster_end = $m->tEnd if($m->tEnd > $cluster_end);
#    }
#
#    # Step 3: expand cluster to the right
#    while($cluster_end > $prev_cluster_end) {
#        $self->_verbose_msg("GeneFactory: getting all mappings starting within $chr:",$prev_cluster_end+1,"-$cluster_end\n");
#	my @some_mappings = $self->_query_db->get_mappings_in_region
#	    (chr => $chr, #strand => $strand,
#	     start => $prev_cluster_end+1, end => $cluster_end,
#	     confine_start => 1, types => $self->{_qTypeL});
#	push @mappings, @some_mappings;
#	$prev_cluster_end = $cluster_end;
#	foreach my $m (@some_mappings) {
#	    $cluster_end = $m->tEnd if($m->tEnd > $cluster_end); 
#	}
#    }
#
#    # Step 4: expand cluster to the left
#    while($cluster_start < $prev_cluster_start) {
#        $self->_verbose_msg("GeneFactory: getting all mappings ending within $chr:$cluster_start-",$prev_cluster_start-1,"\n");
#	my @some_mappings = $self->_query_db->get_mappings_in_region
#	    (chr => $chr, #strand => $strand,
#	     start => $cluster_start, end => $prev_cluster_start-1,
#	     confine_end => 1, types => $self->{_qTypeL});
#	push @mappings, @some_mappings;
#	$prev_cluster_start = $cluster_start;
#	foreach my $m (@some_mappings) {
#	    $cluster_start = $m->tStart if($m->tStart < $cluster_start); 
#	}
#    }
#    $self->_verbose_msg("GeneFactory: ...got all mappings\n");
#
#    # Create gene using the mappings
#    my $genes = $self->_create_genes($start, $end, $strand,
#				     $cluster_start, $cluster_end, \@mappings);
#
#    # Determine start of remaining subregion
#    my $new_start = ($cluster_end < $end) ? $cluster_end+1 : undef;
#
#    # Return genes and the new start
#    return ($genes, $new_start);
#}


1;
