package AT::CMP::CrossMatchFactory;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use AT::AlignmentFactory;
use AT::CMP::CrossMatch;

use constant DEFAULT_REL_START => -2000;
use constant DEFAULT_REL_END => 0;

@ISA = qw(AT::Root);


sub new  {
    my ($caller, %args) = @_;
    my $self = bless { assembly1_db => $args{'assembly1_db'},
		       mapping1_db  => $args{'mapping1_db'},
		       assembly2_db => $args{'assembly2_db'},
		       mapping2_db  => $args{'mapping2_db'},
		       min_coverage => ($args{'min_coverage'} || 0),
		       min_identity => ($args{'min_identity'} || 0),
		       rel_start => ($args{'rel_start'} || DEFAULT_REL_START),
		       rel_end   => ($args{'rel_end'  } || DEFAULT_REL_END)#,
		       #_mappingL    => []
		      
		   }, ref $caller || $caller;
    


    return $self;
}



sub run  {
    my ($self, %args) = @_;
    my ($mapping1, $mapping2);
    if ($args{'acc1'})  {
	my @mappings = $self->mapping1_db->get_mappings_for_acc($args{'acc1'});
	foreach my $m (@mappings) {
	    if ($self->_is_mapping_good_enough($m))  {
		$mapping1 = $m;
		last;
	    }
	}
    }
    elsif ($args{'mapping1'})  {
	$mapping1 = $args{'mapping1'};
    }
    else  {
	croak "No acc1 or mapping1 argument";
    }
    return undef unless $mapping1; # not found
    $mapping1 = $self->extend_5prime($mapping1, $self->mapping1_db) 
	if $args{'extend'};
    
    if ($args{'acc2'})  {
	my @mappings = $self->mapping2_db->get_mappings_for_acc($args{'acc2'});
	foreach my $m (@mappings) {
	    if ($self->_is_mapping_good_enough($m))  {
		$mapping2 = $m;
		last;
	    }
	}
    }
    elsif ($args{'mapping2'})  {
	$mapping2 = $args{'mapping2'};
    }
    else  {
	croak "No acc2 or mapping2 argument";
    }
    unless ($mapping2) { #temporary
	return undef;}
    $mapping2 = $self->extend_5prime($mapping2, $self->mapping2_db) 
	if $args{'extend'};;

    unless ($mapping2) {
	return undef;
	#$mapping2 = $self->find_best_crossmatch2($mapping1);
    }

    return undef unless $mapping2;
    #############################################################
    # get alignment, exons and introns for both mappings

    # 1
    
    my @alignments;

    my $af = AT::AlignmentFactory->new(query_db  => $self->mapping1_db,
				       target_db => $self->assembly1_db);
    $af->run(mapping=>$mapping1);

    @alignments = $af->get_alignments();
    return undef unless @alignments;
    my $aln1 = shift @alignments;


    my ($start1, $end1)  = $aln1->rel_to_abs_target
	($self->rel_start, $self->rel_end); 

    my $genomic1 = 
	$self->assembly1_db->get_genome_seq(start=>$start1, 
					    end=>$end1,
					    strand=>$mapping1->strand,
					    chr=>$mapping1->tName);

    # 2

    $af = AT::AlignmentFactory->new(query_db  => $self->mapping2_db,
				    target_db => $self->assembly2_db);


    $af->run( mapping => $mapping2 );

    @alignments = $af->get_alignments();
    return undef unless @alignments;
    my $aln2 = shift @alignments;

    my ($start2, $end2)  = $aln2->rel_to_abs_target
	($self->rel_start, $self->rel_end); 

    my $genomic2 = 
	$self->assembly1_db->get_genome_seq(start=>$start2, 
					    end=>$end2,
					    strand=>$mapping2->strand,
					    chr=>$mapping2->tName);




    ##############################################################
    # generate crossmatch object



    my $crossmatch = 
	AT::CMP::CrossMatch->new(mapping1   => $mapping1,
				 mapping2   => $mapping2,
				 alignment1 => $aln1,
				 alignment2 => $aln2,
				 genomic1   => $genomic1,
				 genomic2   => $genomic2);

    return $self->{'crossmatch'} = $crossmatch;
    
}


sub extend_5prime  {
    my ($self, $orig_mapping, $mapping_db) = @_;
    my $the_mapping = $orig_mapping;
    # printf("ORIGINAL: %s: %d - %d [%s] %s %d\n", $orig_mapping->tName,
#	    $orig_mapping->tStart,
#	    $orig_mapping->tEnd,
#	    $orig_mapping->strand,
#	   $orig_mapping->qName,
#	   $orig_mapping->matches);
    my @mappings = 
	$mapping_db->get_mappings_in_region
	  (chr   => $the_mapping->tName,
	   start => $the_mapping->tStart,
	   end   => $the_mapping->tEnd,
	   confine => 0,
	   target_db => $orig_mapping->target_db);
    foreach my $m (@mappings) {
	if (
	    $self->_is_mapping_good_enough($m)
	    and
	    $self->_test_equivalence($m, $the_mapping)
	    and
	    (($m->strand eq "+" and $m->tStart < $the_mapping->tStart)
	     or
	    ($m->strand eq "-" and $m->tEnd > $the_mapping->tEnd)
	    ))
	{ $the_mapping = $m; }
    }
    if  ($the_mapping eq $orig_mapping) {
	return $the_mapping;
    }
    else  {
	$self->extend_5prime($the_mapping, $mapping_db);
    }
}


sub find_best_crossmatch2  {
    my ($self, $mapping1) = @_;
    my @cross_matches = 
	$self->get_mappings_for_acc($mapping1->qName, 
				    $self->assembly2_db->dbname);
    return undef unless @cross_matches;
    my $mapping1to2 = shift @cross_matches;
    
    my @mappings2 = 
	$self->mapping2_db->get_mappings_in_region
	  (start  =>$mapping1to2->start,
	   end=>$mapping1to2->end,
	   target_db => $self->assembly2_db->dbname);
    if (my $best_mapping = shift @mappings2) {
	return $self->extend_5prime($best_mapping);
    }
    else  {
	return undef;
    }

}


sub _test_equivalence  {
    my ($self, $map1, $map2) = @_;
    if ($map1->strand eq $map2->strand)  {return 1;}
    else {return 0;}
    # supposed to check whether the two mappings have a comparable 
    # intron/exon structure

    # later
   
}

sub _is_mapping_good_enough {
    my ($self, $mapping) = @_;
    if (($mapping->matches+$mapping->misMatches) / $mapping->qSize 
	< $self->min_coverage) { return 0; }
    elsif ($mapping->matches / ($mapping->matches+$mapping->misMatches)
	   < $self->min_coverage)  { return 0; }
    else { return 1; }
}



1;
