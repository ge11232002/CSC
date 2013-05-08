# AT::FT::Factory::GeneMachine module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::FT::Factory::GeneMachine module

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::FT::Factory::GeneMachine;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use AT::Prediction::Gene;
use AT::Prediction::Exon;
use AT::Prediction::Intron;
use AT::Prediction::SpliceForm;
use AT::FT::TSS;
use AT::FT::PAS;
use AT::FT::Gene::Dynamic;


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
	min_pep_score => ($args{min_pep_score} || 25),
	min_string_score => ($args{min_string_score} || 25),
	gene_stringency => ($args{gene_stringency} || 0),
	min_gfmapping_score_sum => ($args{min_gfmapping_score_sum} || 0)
    }, ref $caller || $caller;
    
    return $self;
}


sub run
{
    my ($self, %args) = @_;
    my $pep_graphs = $args{pep_graphs} || croak "No pep_graphs arg";
    my $target_seq = $args{target_seq} || croak "No target_seq arg";
    my $strand = $args{strand} || croak "No strand arg";
    my $excluded_mappings = $args{excluded_mappings} || croak "No excluded_mappings arg";

    my $genes = $self->_pepgraphs2genes($pep_graphs, $target_seq, $strand,
					$excluded_mappings);
    
}


sub _pepgraphs2genes
{
    my ($self, $pep_graphs, $target_seq, $strand, $excl) = @_;
    my $min_pep_score = $self->min_pep_score;
    my $min_string_score = $self->min_string_score;
    my $min_gfmap_score_sum = $self->min_gfmapping_score_sum;
    my @genes;
    foreach my $pep_graph (@$pep_graphs) {
        my $pep_sets = $pep_graph->visit(min_pep_score => $min_pep_score,
					 min_string_score => $min_string_score);
	my $gfmappings = $self->_get_gfmappings($pep_sets);
	#my ($unique_mappings, $nonunique_mappings) =
	#    $self->_find_unique_and_nonunique_mappings($gfmappings);
	for my $i (0..@$pep_sets-1) {
	    next if($min_gfmap_score_sum and
		    $self->_gfmap_score_sum($gfmappings->[$i]) < $min_gfmap_score_sum);
	    push @genes, $self->_peps2gene($pep_sets->[$i],
					   $target_seq,
					   $strand,
					   undef, #$unique_mappings->[$i],
					   undef, #$nonunique_mappings->[$i],
					   $gfmappings->[$i],
					   $excl);
    	}
    }
    return \@genes;
}


sub _get_gfmappings
# could put this in gene::dynamic
{
    my ($self, $pep_sets) = @_;
    my @ps_gfmappings;
    foreach my $peps (@$pep_sets) {
	my %m;
	foreach my $m 
	    (map { $_->mapping_list }
	     @$peps) {
	    $m{$m} = $m;
        }
	push @ps_gfmappings, [values(%m)];
	#print STDERR "PS ",scalar(@$peps)," : ", join(' ', map {$_->qName} values(%m)), "\n";
    }
    return \@ps_gfmappings;
}


sub _gfmap_score_sum
{
    my ($self, $gfmappings) = @_;
    my $sum = 0;
    foreach my $m (@$gfmappings) {
	$sum += $m->score1; 
    }
    return $sum;
}


sub _find_unique_and_nonunique_mappings
{
    my ($self, $ps_gfmappings) = @_;

    # get non-redundant set of mappings for each pep set
    my @ps_mappings;
    foreach my $gfmapping_set (@$ps_gfmappings) {
	push @ps_mappings, [ map { $_->primary_mapping_list } @$gfmapping_set ];
    }
#    foreach my $peps (@$pep_sets) {
#	my %m;
#	foreach my $m 
#	    (map { $_->primary_mapping_list } # use primary mappings for now since impl in AT::Predicton::Gene
#	     map { $_->mapping_list }
#	     #map { $_->mapping }
#	     #map { ($_->open_support_list, $_->close_support_list) }
#	     @$peps) {
#	    $m{$m} = $m;
#        }
#	push @ps_mappings, [values(%m)];
#	#print STDERR "PS ",scalar(@$peps)," : ", join(' ', map {$_->qName} values(%m)), "\n";
#    }

    # count the occurence of each mapping
    my %count;
    foreach my $m (map { @{$_} } @ps_mappings) {
	$count{$m}++;
    }

    # for each pep set, partition mappings into unique and non-unique
    my (@ps_u, @ps_nu);
    for my $i (0..@ps_mappings-1) {
	my (@u, @nu);
	foreach my $m (@{$ps_mappings[$i]}) {
	    if($count{$m} == 1) { push @u, $m; }
	    else { push @nu, $m; }
	}
	$ps_u[$i] = \@u;
	$ps_nu[$i] = \@nu;
    }

#    for my $i (0..@ps_mappings-1) {
#        print STDERR "PSU ",scalar(@{$pep_sets->[$i]}), " : ",
#	    join(' ', map {$_->qName} @{$ps_u[$i]}), "\n";
#    }
#
#    for my $i (0..@ps_mappings-1) {
#        print STDERR "PSNU ",scalar(@{$pep_sets->[$i]}), " : ",
#	    join(' ', map {$_->qName} @{$ps_nu[$i]}), "\n";
#    }

    return (\@ps_u, \@ps_nu);
}


sub _peps2gene {
    my ($self, $peps_ref, $target_seq, $strand, $unique_mappings,
	$nonunique_mappings, $gfmappings, $excl) = @_;

    my @peps = sort {$a->start <=> $b->start} @$peps_ref;
    return unless (@peps);

    my $is_spliced = 0;
    foreach my $m (@$gfmappings) {
	if($m->nr_introns) { $is_spliced = 1; last; }
    }

    my $gene_start = $peps[0]->start;
    my $gene_end = $peps[-1]->end;
    my $seqstr = $target_seq->subseq($gene_start - $target_seq->start + 1,
				    $gene_end - $target_seq->start + 1);
    my $seq = Bio::LocatableSeq->new(-id => $target_seq->id,
				    -start => $gene_start,
				    -end => $gene_end,
				    -strand => $strand,
				    -seq => $seqstr);

    my $gene = AT::FT::Gene::Dynamic->new
	(seq => $seq,
	 peps => \@peps,
	 strand => $strand,
	 stringency => $self->gene_stringency,
	 gfmapping_list => $gfmappings,
	 is_spliced => $is_spliced
	 #unique_mapping_list => $unique_mappings,
	 #nonunique_mapping_list => $nonunique_mappings
    );

return $gene;
			
    #my $sf = AT::Prediction::SpliceForm->new(mapping_list => $unique_mappings);
    my $gene = AT::Prediction::Gene->new
	(unique_mapping_list => $unique_mappings,
	 nonunique_mapping_list => $nonunique_mappings,
	 seq => $seq);
    # the non-unique mappings should be annotated somehow; put them in excl for now
    # let's forget about the actual excluded mappings; can't keep track of them anyway

    # Create exons and introns
    my (@exons, @introns);
    my $i = 0;
    while ($i < @peps) {
	my $j = $i+1;
	while($j < @peps and $peps[$j]->start == $peps[$j-1]->end + 1) {
	    $j++;
	}
	# make exon from peps $i..$j-1
	#print STDERR "E: ",$peps[$i]->start,"-",$peps[$j-1]->end,"\n";
	push @exons, AT::Prediction::Exon->new
	    ( -start => $peps[$i]->start - $gene_start + 1,
	      -end => $peps[$j-1]->end - $gene_start + 1,
	      -strand => $strand,
	      -abs_start => $peps[$i]->start,
	      -abs_end => $peps[$j-1]->end
	      );
	if($j < @peps) {
	    #print STDERR "I: ",$peps[$j-1]->end+1,"-",$peps[$j]->start-1,"\n";
	    push @introns, AT::Prediction::Intron->new
		( -start => $peps[$j-1]->end - $gene_start + 2,
		-end => $peps[$j]->start - $gene_start,
		-strand => $strand,
		-abs_start => $peps[$j-1]->end + 1,
		-abs_end => $peps[$j]->start - 1
		);
	}
	$i = $j;
    }
    $gene->set_exons_introns(\@exons, \@introns);
    $self->_dbg_print_open_close(\@peps, 2);

    my @open_pos =
	map { [$_->start, $_->open_score] }
	grep { $_->open_score >= 20 } @peps;
    my @close_pos =
	map { [$_->end, $_->close_score] }
	grep { $_->close_score >= 20 } @peps;
    my ($putative_TSS_pos, $putative_PAS_pos) =
	($strand == 1) ?
	(\@open_pos, \@close_pos) :
	(\@close_pos, \@open_pos);
    my $TSS = $self->_create_borders($putative_TSS_pos,
				     $gene_start,
				     $strand,
				     'AT::FT::TSS');
    my $PAS = $self->_create_borders($putative_PAS_pos,
				     $gene_start,
				     $strand,
				     'AT::FT::PAS');
    $gene->TSS_list($TSS);
    $gene->PAS_list($PAS);

    if($self->debug) {
        print STDERR "Created gene: ", $gene->loc_str, "\n";
        $gene->print_exons_introns(\*STDERR);
        $gene->print_borders(\*STDERR);
        print STDERR "\n";
    }

    return $gene;
}


sub _create_borders
{
    my ($self, $a, $gene_start, $strand, $border_class) = @_;
    my $max_d = 10;
    my @borders;
    for(my $i = 0; $i < @$a;) {
	my ($start, $score) = @{$a->[$i]};
	if ($score >= 3 or ($i+1 < @$a and $a->[$i+1][0] == $start+1))
	{
	    while($i+1 < @$a and $a->[$i+1][0] <= $a->[$i][0] + $max_d) {
		$i++;
		$score += $a->[$i][1];
	    }
	    my $end = $a->[$i][0];
	    my ($rel_start, $rel_end) = ($start - $gene_start + 1,
					 $end - $gene_start + 1);
	    push @borders, $border_class->new
		(-start => $rel_start,
		 -end => $rel_end,
		 -strand => $strand,
		 -abs_start => $start,
		 -abs_end => $end,
		 -score => $score);
	}
	$i++;
    }
    return \@borders;
}


#sub _create_borders
#{
#    my ($self, $a, $dir, $gene_start, $strand, $border_class) = @_;
#    my $max_d = 9;
#    my @borders;
#    for(my $i = 0; $i < @$a;) {
#	my ($pos, $score) = @{$a->[$i]};
#	if ($score == 3 or ($i+1 < @$a and $a->[$i+1][0] == $pos+$dir))
#	{
#	    while($i+1 < @$a and $a->[$i+1][0] <= $a->[$i][0] + $max_d
#		             and $a->[$i+1][0] >= $a->[$i][0] - $max_d) {
#		$i++;
#		$score += $a->[$i][1];
#	    }
#	    my $rel_pos = $pos - $gene_start + 1;
#	    push @borders, $border_class->new
#		(-start => $rel_pos,
#		 -end => $rel_pos,
#		 -strand => $strand,
#		 -abs_start => $pos,
#		 -abs_end => $pos,
#		 -score => $score);
#	}
#	$i++;
#    }
#    return \@borders;
#}
#

sub _dbg_print_open_close
{
    my ($self, $peps, $min_score) = @_;
    return unless ($self->debug);
    foreach my $pep (@$peps) {
	my $open_score = $pep->open_score;
	my $close_score = $pep->close_score;
	print STDERR "OPEN $open_score\t",$pep->start,"\t",join(',',@{$pep->{_open_qNames}}),"\n" if($open_score >= $min_score);	
	print STDERR "CLOSE $close_score\t",$pep->end,"\t",join(',',@{$pep->{_close_qNames}}),"\n" if($close_score >= $min_score);
    }
}

1;
