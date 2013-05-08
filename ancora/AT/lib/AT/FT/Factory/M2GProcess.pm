# AT::FT::Factory::M2GProcess module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD


package AT::FT::Factory::M2GProcess;

use strict;
use vars '@ISA';
use Carp;
use Bio::Graphics::Panel;
use AT::Root;
use AT::FT::Factory::MappingEndMover;
use AT::FT::Factory::MappingEndMover_old;
use AT::FT::Factory::MappingMerger;
use AT::FT::Factory::GFMappingMachine;
use AT::FT::Factory::MappingCompass;
use AT::FT::Factory::MappingEndValidator;
use AT::FT::Factory::MappingScorer;
use AT::FT::Factory::MappingGrouper;
use AT::FT::Factory::PEPGraphMachine;
use AT::FT::Factory::GeneMachine;
use AT::GFX::PEPGraph::ScoreCurve;
use File::Temp 'tempfile';
use Data::Dumper;

use constant PROBLEM_TAG_METHOD => scalar 'main';

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
	show_chr_panel => ($args{'show_chr_panel'} || 0),
	write_chr_panel => ($args{'show_chr_panel'} || 0),
	show_pep_graph => ($args{'show_pep_graph'} || 0),
	write_pep_graph => ($args{'write_pep_graph'} || 0),
	verbose => ($args{'verbose'} || 0),
	gf_mapping_machine => AT::FT::Factory::GFMappingMachine->new(),
	mapping_merger => AT::FT::Factory::MappingMerger->new(),
	mapping_end_mover => AT::FT::Factory::MappingEndMover->new(),
	mapping_compass => AT::FT::Factory::MappingCompass->new(),
	mapping_end_validator => AT::FT::Factory::MappingEndValidator->new(),
	mapping_scorer => AT::FT::Factory::MappingScorer->new(),
	mapping_grouper => AT::FT::Factory::MappingGrouper->new(),
	pep_graph_machine => AT::FT::Factory::PEPGraphMachine->new(),
	gene_machine => AT::FT::Factory::GeneMachine->new(),
	do_endmoving => (defined($args{'do_endmoving'}) ? $args{'do_endmoving'} : 1),
	do_grouping => (defined($args{'do_grouping'}) ? $args{'do_grouping'} : 1),
	do_scoring => (defined($args{'do_scoring'}) ? $args{'do_scoring'} : 1),
	do_endvalidation => (defined($args{'do_endvalidation'}) ? $args{'do_endvalidation'} : 1),
	finish_at_gfmappings => ($args{'finish_at_gfmappings'} || 0),
	print_excluded_gfmappings => ($args{'print_excluded_gfmappings'} || 0),
	return_excluded_gfmappings => ($args{'return_excluded_gfmappings'} || 0),

    }, ref $caller || $caller;
    
    return $self;
}

# instead of the do_x methods, we should create a switch for each machine
# e.g. $gf->m2g_process->mapping_end_mover->is_active(0);

sub run
{
    my ($self, %args) = @_;
    my $mappings = $args{mappings} || croak "No mappings arg";
    my $target_seq = $args{target_seq} || croak "No target_seq arg";
    my $strand = $args{strand} || 0;

    # init chr panel
    $self->_init_chr_panel($mappings, $target_seq);

    # create gfmappings
    my ($plus_mappings, $minus_mappings, $excluded_mappings) =
	$self->get_gfmappings(%args);
    if($self->print_excluded_gfmappings) {
        foreach my $m (@$excluded_mappings) {
	    print STDERR "MAPPING EXCLUDED: mRNA=", $m->mRNA_acc_string, " EST=", $m->EST_acc_string, "; ", $m->loc_str, 
	      "; REASONS: ", join(',', $m->problem_tag_list), "\n";
	}
    }
    if($self->finish_at_gfmappings) {
        return $self->return_excluded_gfmappings ? (@$plus_mappings, @$minus_mappings, @$excluded_mappings)
						 : (@$plus_mappings, @$minus_mappings);
    }
    

    # create genes
    my @genes = $self->create_genes_from_gfmappings($plus_mappings, $minus_mappings, $excluded_mappings, $target_seq);

    $self->_show_chr_panel();

    return @genes;
}


sub run_from_gfmappings
{
    my ($self, %args) = @_;
    my $mappings = $args{gfmappings} || croak "No gfmappings arg";
    my $target_seq = $args{target_seq} || croak "No target_seq arg";

    my (@plus_mappings, @minus_mappings);
    foreach my $m (@$mappings) {
	if($m->strand == 1) {
	    push @plus_mappings, $m;
	}
	else {
	    push @minus_mappings, $m;
	}
    }
    $self->_init_chr_panel($mappings, $target_seq);
    my @genes = $self->create_genes_from_gfmappings(\@plus_mappings, \@minus_mappings, [], $target_seq);
    $self->_show_chr_panel();
    return @genes;
}


sub get_gfmappings
{
    my ($self, %args) = @_;
    my $mappings = $args{mappings} || croak "No mappings arg";
    my $target_seq = $args{target_seq} || croak "No target_seq arg";
    my $strand = $args{strand} || 0;

    # create and orient mappings
    $self->_verbose_msg("Creating GFMappings\n");
    my $gf_mappings = $self->gf_mapping_machine->run(mappings => $mappings,
	      				             target_seq => $target_seq);
    $self->_verbose_msg("Moving mapping ends\n");
    $self->mapping_end_mover->run(mappings => $gf_mappings, target_seq => $target_seq)
	if($self->do_endmoving);
	
    #my $excluded_mappings;
    #($gf_mappings, $excluded_mappings) = $self->_exclude_low_quality_mappings($gf_mappings);   # exclude mappings with <75% of query mapped
    $self->_verbose_msg( "Merging mappings\n");
    $gf_mappings = $self->mapping_merger->run(mappings => $gf_mappings);
    $self->_verbose_msg( "Orienting mappings\n");
    my ($plus_mappings, $minus_mappings, $excluded_mappings) =
	$self->mapping_compass->run(mappings => $gf_mappings,
		       target_seq => $target_seq);
    if($strand == 1) { $minus_mappings = []; }
    elsif($strand == -1) { $plus_mappings = []; }
    #push @$excluded_mappings, @$more_excluded_mappings;

    if($self->do_endvalidation) {
        $self->_verbose_msg( "Validating mapping ends\n");
        $self->mapping_end_validator->run(mappings => $plus_mappings);
        $self->mapping_end_validator->run(mappings => $minus_mappings);
    }
    if($self->do_scoring) {
        $self->_verbose_msg( "Scoring mappings\n");
        $self->mapping_scorer->run(mappings => $plus_mappings);
        $self->mapping_scorer->run(mappings => $minus_mappings);  
    }

    return ($plus_mappings, $minus_mappings, $excluded_mappings);
}


sub create_genes_from_gfmappings
{
    my ($self, $plus_mappings, $minus_mappings, $excluded_mappings, $target_seq) = @_;

    # show mappings on panel
    $self->_verbose_msg( "Drawing mappings\n");
    $self->_draw_gf_mappings_on_panel
	('Excluded mappings', @$excluded_mappings);
    unless($self->do_grouping) {
	$self->_draw_gf_mappings_on_panel
	    ('Minus strand mappings', @$minus_mappings);
	$self->_draw_gf_mappings_on_panel
	    ('Plus strand mappings', @$plus_mappings);
    }

    # create genes on plus strand
    $self->_verbose_msg("Creating for + strand\n");
    my $plus_genes = $self->_mappings2genes
	($plus_mappings, 1, $target_seq, $excluded_mappings);
    # create genes on minus strand
    $self->_verbose_msg("Creating for - strand\n");
    my $minus_genes = $self->_mappings2genes
	($minus_mappings, -1, $target_seq, $excluded_mappings);
    # remove likely FP genes
    $self->_verbose_msg("Removing reversed genes");
    ($plus_genes, $minus_genes) = $self->_remove_reversed_genes($plus_genes, $minus_genes);
    # show genes on panel
    $self->_verbose_msg( "Drawing genes\n");
    $self->_draw_genes_on_panel("- genes", $minus_genes, "+ genes", $plus_genes);
    #$self->_draw_extra_genes($minus_genes, $plus_genes);

    return (@$plus_genes, @$minus_genes);
}


sub _exclude_low_quality_mappings
# exclude mappings with <150nt and <75% of query sequence mapped
# only works on mappings with one primary for now
{
    my ($self, $mappings1) = @_;
    my (@incl, @excl);
    foreach my $m (@$mappings1) {
	my $mapped_nt = 0;
	foreach my $hsp ($m->HSP_list) {
	    $mapped_nt += $hsp->length;
	}
	if($mapped_nt >= 150) {
	    push @incl, $m;
	}
	else {
	    my $seqstr = ($m->primary_mapping_list)[0]->query_seq->seq;
	    $seqstr =~ s/^T*//i;  $seqstr =~ s/A*$//i;
	    if($mapped_nt >= 0.75*length($seqstr)) {
		push @incl, $m;
	    }
	    else {
		push @excl, $m;
		$m->add_problem_tag(PROBLEM_TAG_METHOD, 'fraction_mapped');
	    }
	}
    }
    return (\@incl, \@excl);
}


sub _draw_extra_genes
{
    my ($self, $minus_genes, $plus_genes) = @_;
    my (@histr_minus_genes, @histr_plus_genes);
    my $orig_stringency;
    foreach my $gene (@$minus_genes) {
	$orig_stringency = $gene->stringency;
	$gene->stringency(25);
	push @histr_minus_genes, $gene->primary_gene if($gene->primary_gene);
    }
    foreach my $gene (@$plus_genes) {
	$orig_stringency = $gene->stringency;
	$gene->stringency(25);
	push @histr_plus_genes, $gene->primary_gene if($gene->primary_gene);
    }
    $self->_draw_genes_on_panel("- genes lostr", \@histr_minus_genes, "+ genes lostr", \@histr_plus_genes);
    my (@nostretch_minus_genes, @nostretch_plus_genes);
    foreach my $gene (@$minus_genes) {
	$gene->use_score(2);
	push @nostretch_minus_genes, $gene->primary_gene if($gene->primary_gene);
	$gene->use_score(1);
	$gene->stringency($orig_stringency);
    }
    foreach my $gene (@$plus_genes) {
	$gene->use_score(2);
	push @nostretch_plus_genes, $gene->primary_gene if($gene->primary_gene);
	$gene->use_score(1);
	$gene->stringency($orig_stringency);
    }
    $self->_draw_genes_on_panel("- genes lostr no_stretch", \@nostretch_minus_genes, "+ genes lostr no_stretch", \@nostretch_plus_genes);

}


sub _mappings2genes
{
    my ($self, $mappings, $strand, $target_seq, $excluded_mappings) = @_;

    # Group mappings
    my $groups;
    if($self->do_grouping) {
        $self->_verbose_msg( "Grouping mappings\n");
        $groups = $self->mapping_grouper->run(mappings => $mappings,
						strand => $strand);
        $self->_verbose_msg( "Drawing groups\n");
        $self->_draw_map_groups_on_panel($groups, $strand);
    }
    else {
	$groups = [$mappings];
    }

    # Make pep graphs from groups
    $self->_verbose_msg( "Making PEP graph\n");
    my @pep_graphs =
	map { $self->pep_graph_machine->run(mappings => $_) } @$groups;
    for my $i (0..@pep_graphs-1) {
	$pep_graphs[$i]->prune(min_pep_score => 3,
			       min_string_score => 3);
	$self->_display_pep_graph($pep_graphs[$i], $i+1, $strand);
	#$self->_display_pep_graph_curve($pep_graphs[$i]) ;#if ($i == 0);
    }

    # Make genes from pep graphs
    $self->_verbose_msg( "Making genes\n");
    my $genes = $self->gene_machine->run(pep_graphs => \@pep_graphs,
					 strand => $strand,
					 target_seq => $target_seq,
					 excluded_mappings => $excluded_mappings);

    return $genes;
}



# * the reversed-gene-removal methods should be modularized * 

sub _remove_reversed_genes
{
    my ($self, $genes1, $genes2) = @_;
    my $kept_genes1 = $self->_remove_reversed_genes_on_one_strand($genes1, $genes2);
    my $kept_genes2 = $self->_remove_reversed_genes_on_one_strand($genes2, $genes1);
    return ($kept_genes1, $kept_genes2);
}

sub _remove_reversed_genes_on_one_strand
{
    my ($self, $genes1, $genes2) = @_;
    my @kept_genes;

    foreach my $g (@$genes1) {

	# Keep gene by default
	my $keep_gene = 1;

	# Is gene supported by unspliced ESTs only?
	my $only_unspliced_est_support = 1;
	foreach my $m ($g->gfmapping_list) {
	    if($m->mRNA_acc_list or $m->is_spliced) {
		$only_unspliced_est_support = 0;
		last;
	    }
	}

	# If yes, go on...
	if($only_unspliced_est_support) {

	    # Get number of ESTs supporting gene
	    my $nr_ests = $g->gfmapping_list;

	    # Get number of ESTs in opposite, overlapping genes
	    my $nr_opposite_ests = 0;
	    foreach my $opposite_gene (@$genes2) {
		if($self->_genes_have_exon_overlap($g,$opposite_gene)) {
		    foreach my $m ($opposite_gene->gfmapping_list) {
			$nr_opposite_ests++ if($m->EST_acc_list);
		    }
		}
		print STDERR "  ", $opposite_gene->loc_str, " ", $nr_opposite_ests, "\n";
	    }
	
	    # Check balance
	    $keep_gene = 0 unless($self->_est_balance($nr_ests, $nr_opposite_ests)); 

	    print STDERR "M2GProcess: gene ", $g->loc_str, "  :  $nr_ests -- $nr_opposite_ests  :  ", ($keep_gene ? "kept\n" : "tossed\n")
		if(1);#if($self->debug);
	}

	push @kept_genes, $g if($keep_gene);
    }
    return \@kept_genes;
}


sub _est_balance
{
    my ($self,$est_count1,$est_count2) = @_;
    return 1 if($est_count1 >= $est_count2);
    my $thresholds = $self->{_est_balance_thresholds};
    unless($thresholds) {
	$thresholds = $self->{_est_balance_thresholds} =
	    [ 4,    72,   215,   408,   635,   888,  1159,  1447,  1746,  2057,  2377,  2704,
	      3039,  3380,  3726,  4077,  4433,  4793,  5157,  5525,  5895,  6269,  6645,  7024,
	      7406,  7790,  8176,  8564,  8954,  9346,  9740, 10135, 10532, 10931, 11331, 11732,
	      12135, 12539, 12944, 13351, 13758, 14167, 14577, 14988, 15400, 15812, 16226, 16641,
	      17056, 17473, 17890, 18308, 18727, 19146, 19566, 19987, 20409, 20831, 21254, 21678,
	      22102, 22527, 22953, 23379, 23805, 24233, 24660, 25089, 25517, 25947, 26376, 26807,
	      27237, 27669, 28100, 28532, 28965, 29398, 29831, 30265, 30699, 31134, 31569, 32004,
	      32440, 32876, 33313, 33750, 34187, 34624, 35062, 35501, 35939, 36378, 36817, 37257,
	      37697, 38137, 38577, 39018
	     ];
    }
    my $t = $thresholds->[$est_count1-1];
    return (!$t or $est_count2 < $t) ? 1 : 0;
}

# In the above method, the thresholds are hardcoded and they should of course not be.
# Can be computed upon instantiation or first method call.
# Parameters: EST orienataion error rate; conf level
# (hardcoded table computed for error rate = 0.002 and conf level = 0.95)
# Array index =  nr ESTs in gene1-1
# Array value =  max nr ESTs in gene2 allowed for gene1 to be kept
#
# R code used to compute array:
# revProb <- 0.002;
# confLevel <- 0.99;
# xMax <- 100;
# result <- vector(length=xMax);
# y <- 0;
# for (x in 1:xMax) {
#  while(pbinom(x-1, x+y+1, revProb) > confLevel) y <- y+1;
#  result[x] <- y;
# }

# Thresholds for confLevel 0.95:
#	    [   24,   176,   406,   679,   981,  1301,  1637,  1983,  2340,  2704,  3075,  3452,
#              3833,  4220,  4610,  5004,  5401,  5801,  6204,  6610,  7017,  7427,  7839,  8253,
#              8669,  9086,  9505,  9925, 10347, 10770, 11195, 11620, 12047, 12475, 12903, 13333,
#	     13764, 14196, 14628, 15062, 15496, 15931, 16367, 16803, 17241, 17679, 18117, 18557,
#	     18996, 19437, 19878, 20320, 20762, 21204, 21648, 22091, 22536, 22980, 23426, 23871,
#	     24317, 24764, 25211, 25658, 26106, 26554, 27003, 27451, 27901, 28350, 28800, 29251,
#	     29701, 30152, 30604, 31055, 31507, 31960, 32412, 32865, 33318, 33772, 34225, 34679,
#	     35134, 35588, 36043, 36498, 36953, 37409, 37864, 38320, 38777, 39233, 39690, 40147,
#	     40604, 41061, 41519, 41977 ];	


sub _genes_have_exon_overlap
{
    my ($self,$g1,$g2) = @_;

    return 0 unless($g1->abs_start <= $g2->abs_end and $g1->abs_end >= $g2->abs_start);
    # Note: we don't check that genes are on same chr

    my @e1 = $g1->exon_list;
    my @e2 = $g2->exon_list;

    # sort in ascending absolute coord order (i.e. reverse array for - strand)
    @e1 = reverse @e1 if($g1->strand == -1);
    @e2 = reverse @e2 if($g2->strand == -1);

    my($i, $j) = (0, 0);

    while($i < @e1 and $j < @e2) {
	my $ol_start = ($e1[$i]->abs_start > $e2[$j]->abs_start) ?
	    $e1[$i]->abs_start : $e2[$j]->abs_start;
	my $ol_end = ($e1[$i]->abs_end < $e2[$j]->abs_end) ?
	    $e1[$i]->abs_end : $e2[$j]->abs_end;
	return 1 if($ol_start <= $ol_end);
	if ($e1[$i]->abs_end < $e2[$j]->abs_end) {
	    $i++;
	}
	else {
	    $j++;
	}
    }

    return 0;
}


# ******************
# ******************
# ** DIAGNOSTICS ***
# ******************
# ******************


sub _count_spliced {
    my ($self, $mappings) = @_;
    my $count = 0;
    foreach my $m (@$mappings) {
	my $is_spliced;
	for my $i (2..$m->nr_HSPs) {
	    $is_spliced = 1 if($m->HSP($i)->left_gap->jnc_type ne 'other');	    
	}
	$count++ if($is_spliced);
    }
    my $tot = scalar(@$mappings);
    my $pc = int(100*$count/$tot+.5);
    return($count, $tot, $pc);
}


# ****************
# ****************
# ** GFX STUFF ***
# ****************
# ****************

sub _init_chr_panel
{
    my ($self, $mappings, $seq) = @_;

    unless($self->{show_chr_panel} or $self->{write_chr_panel}) {
	$self->{_chr_panel} = 0;
	return;
    }

#    if(@$mappings > 1000) {
#	# Drawing many features is very slow, so skip it
#	warn "_init_chr_panel: recieved ".scalar(@$mappings)." mappings; will not create panel\n";
#	$self->{_chr_panel} = 0;
#	return;
#    }

    $self->{_chr_panel} = Bio::Graphics::Panel->new
         ( -length => $seq->length + 40,
	   -offset => $seq->start-21,
	   -width => 1200,
	   -grid => 1,
	   -pad_left => 50,
	   -pad_right => 50,
	   -pad_top => 10,
	   -pad_bottom => 10,
	   -key_style => 'between'
    );

    $self->{_chr_panel_id} = ($seq->id).'_'.($seq->start).'-'.($seq->end);
}


sub _show_chr_panel
{
    my ($self) = @_;
    my $panel = $self->{_chr_panel} or return;
    
    # Add scale track
    my $full_length = Bio::SeqFeature::Generic->new(-start=> $panel->start,
						    -end=> $panel->end,
						    -primary=>'ruler');
    $panel->unshift_track('arrow',
		      $full_length,
		      -tick => 2,
		      -fgcolor => 'black',
		      -double => 1);
    
    my $png = $panel->png;
    
    $self->_display_png($png) if($self->show_chr_panel);

    if($self->write_chr_panel) {
	open PNG_OUT, '>'.($self->{_chr_panel_id}).'.png';
	print PNG_OUT $png;
	close PNG_OUT;
    }
}


sub _draw_mappings_on_panel
{
    my ($self, @data) = @_;
    my $panel = $self->{_chr_panel} or return;

    my (@tracks, @track_args);

    for(my $i = 0; $i < @data; $i+=2) {
	push @track_args, { key => $data[$i] };
	push @tracks, [ map { $self->_mapping_to_SeqFeature($_) }
		        @{$data[$i+1]} ];
    }

    my %color = ( 'EST' => 'purple', 'mRNA' => 'khaki',
		  'refSeq' => 'red', 'refSeqRev' => 'springgreen');

    foreach my $i (0..@tracks-1) {
	$panel->unshift_track
	    ('segments' => $tracks[$i],
	     -bgcolor => sub { $color{$_[0]->type}; },
	     -label => 1,
	     %{$track_args[$i]}
	    );
    }

}

sub _draw_map_groups_on_panel
{
    my ($self, $groups, $strand) = @_;
    $self->{_chr_panel} or return;
    foreach my $i (0..@$groups-1) {
	$self->_draw_gf_mappings_on_panel
	    ("strand ".($strand==1?'+':'-')." group ".($i+1),
	     @{$groups->[$i]});
    }
}


sub _draw_gf_mappings_on_panel
{
    my ($self, $key, @m) = @_;
    my $panel = $self->{_chr_panel} or return;

    my @feats = map { $self->_gf_mapping_to_SeqFeature($_) } @m;

    $panel->unshift_track
	    ('segments' => \@feats,
	     -label => 1,
	     -glyph => 'segments',
	     -connector => sub { (shift->each_tag_value('connector'))[0] },
	     -bgcolor => sub { shift->type->[0] },
	     -fgcolor => sub { shift->type->[1] },
	     -key => $key
	    );
}


sub _gf_mapping_to_SeqFeature
{
    my ($self, $m) = @_;
    my @unbroken; my $i = 0;
    my $bgcolor = $m->mRNA_acc_list ? 'springgreen' : 'khaki';
    my $fgcolor = $m->internally_primed ? 'lightblue' : 'black';
    foreach my $part ($m->HSP_list) {
	push @{$unbroken[$i]}, $part;
	$i++ if($part->right_gap and
		$part->right_gap->is_broken);
    }
    my $group = Bio::Graphics::Feature->new
	(-type => [$bgcolor, $fgcolor],  # encode color here since attr are forgotten
	 -name => $m->acc_string,
	 -attributes => {connector => 'dashed'});  # attr are forgotten. why?
    foreach my $unbroken_set (@unbroken) {
	$group->add_segment(Bio::Graphics::Feature->new
	    (-type => [$bgcolor, $fgcolor],
	     -segments => [ map { [$_->start, $_->end ] } @$unbroken_set ],
	     -attributes => {connector => 'solid'}));

    }
    return $group;
}


sub _draw_genes_on_panel
{
    my ($self, @data) = @_;
    my $panel = $self->{_chr_panel} or return;

    my (@tracks, @track_args);

    my %glyph = ('gene' => 'segments', 'tss_pas' => 'segments', 'tss' => 'arrow', 'pas' => 'box');
    my %color = ('tss' => 'green', 'pas' => 'red', 'undef_part' => 'lightblue');

    for(my $i = 0; $i < @data; $i+=2) {
	push @track_args, { key => $data[$i] };
	push @tracks, [ map { $self->_gene_to_exon_intron_SeqFeature($_),
			      $self->_gene_to_TSS_PAS_SeqFeature($_)}
		        @{$data[$i+1]} ];
    }

    foreach my $i (0..@tracks-1) {
	$panel->unshift_track
	    ($tracks[$i],
	     -label => 1,
     	     -glyph => sub { $glyph{shift->primary_tag}; },
	     -northeast => sub { shift->strand == 1 ? 1: 0 },
	     -bgcolor => sub { $color{shift->primary_tag} || 'blue' },
	     -fgcolor => sub { $color{shift->primary_tag} || 'black' },
	     -base => 1,
	     %{$track_args[$i]}
	    );
    }
}



sub _mapping_to_SeqFeature
{
    my ($self, $m) = @_;
    my $f = Bio::Graphics::Feature->new
	(-segments => [ map { [$_->tStart, $_->tEnd ] } $m->all_HSPs ],
	 -strand => ($m->strand eq '+') ? 1 : -1,
	 -type => $self->_mapping_category($m),
	 -source_tag => 'BLAT',
	 -name => $m->qName);
    return $f;
}


sub _gene_to_exon_intron_SeqFeature
{
    my ($self, $g) = @_;
    my $gene_ft = Bio::SeqFeature::Generic->new
	(-strand => $g->strand,
	 -primary => 'gene',
	 #-source_tag => 'AT',
	 #-seq_id => $m->tName,
	 -name => $g->loc_str);
    foreach my $exon ($g->exon_list) {
	my $exon_ft = Bio::SeqFeature::Generic->new(-start => $exon->abs_start,
					   -end => $exon->abs_end,
					   -strand => $exon->strand,
					   -primary => 'exon');
	$gene_ft->add_SeqFeature($exon_ft, 'EXPAND');
    }
    foreach my $undef ($g->undef_part_list) {
	my $undef_ft = Bio::SeqFeature::Generic->new(-start => $undef->abs_start,
					   -end => $undef->abs_end,
					   -strand => $undef->strand,
					   -primary => 'undef_part');
	$gene_ft->add_SeqFeature($undef_ft, 'EXPAND');
    }
    return $gene_ft;
}


#sub _gene_to_exon_intron_SeqFeature
#{
#    my ($self, $g) = @_;
#    my $f = Bio::Graphics::Feature->new
#	(-segments => [ map { [ $_->abs_start, $_->abs_end ] } $g->exon_list ],
#	 -strand => $g->strand,
#	 -type => 'segments',
#	 -source_tag => 'AT',
#	 #-seq_id => $m->tName,
#	 -name => $g->loc_str);
#    return $f;
#}


sub _gene_to_TSS_PAS_SeqFeature
{
    my ($self, $g) = @_;
    return () unless($g->TSS_list or $g->PAS_list);
    my $group = Bio::SeqFeature::Generic->new
	(-primary => 'tss_pas',
	 -name => '');
    foreach my $tss ($g->TSS_list) {
	$group->add_SeqFeature(Bio::SeqFeature::Generic->new
	    (-primary => 'tss',
	     -start => $tss->abs_start,
	     -end => $tss->abs_end,
	     -strand => $g->strand), 'EXPAND');
    }
    foreach my $pas ($g->PAS_list) {
	$group->add_SeqFeature(Bio::SeqFeature::Generic->new
	    (-primary => 'pas',
	     -start => $pas->abs_start,
	     -end => $pas->abs_end), 'EXPAND');
    }
    return $group;
}


#sub _gene_to_TSS_PAS_SeqFeature
#{
#    my ($self, $g) = @_;
#    my $group = Bio::Graphics::Feature->new
#	(-type => 'segments',
#	 -name => '');
#    foreach my $tss ($g->TSS_list) {
#	$group->add_segment(Bio::Graphics::Feature->new
#	    (-type => 'triangle',
#	     -start => $tss->abs_start,
#	     -end => $tss->abs_end,
#	     -score => $tss->score,
#	     -strand => $g->strand));
#	print "TSS\n";
#    }
#    foreach my $pas ($g->PAS_list) {
#	$group->add_segment(Bio::Graphics::Feature->new
#	    (-type => 'dot',
#	     -start => $pas->abs_start,
#	     -end => $pas->abs_end,
#	     -score => $pas->score));
#	print "PAS\n";
#    }
#    return $group;
#}


sub _display_pep_graph_curve
{
    my ($self, $pepg) = @_;
    my $img = AT::GFX::PEPGraph::ScoreCurve->new(orientation => 'landscape');
    $img->draw($pepg);
    $img->display;
}

sub _display_pep_graph
{
    my ($self, $pepg, $nr, $strand) = @_;
    return unless ($self->{show_pep_graph} or $self->{write_pep_graph});

    require GraphViz;
    require Tk;
    require Tk::GraphViz;

    my $g = GraphViz->new(); # Create GraphViz object
    my $direction = $strand == 1 ? 'forward' : 'back';

    # Assign a name to each pep and add nodes
    my %names;
    my $major = ord('A');
    my $minor = 1;
    for(my $pep = $pepg->head->next; $pep; $pep = $pep->next) {
	my $name = chr($major).$minor;
	$names{$pep} = $name;
	$g->add_node($name,
		     label => ($pep->end-$pep->start+1)." ".$pep->score,
		     cluster => chr($major),
		     fontsize => 10);
	if($pep->next and $pep->next->start > $pep->end + 1) {
	    $major++;
	    $minor = 1;
	}
	else {
	    $minor++;
	}
    }

    # Add edges
    for(my $pep = $pepg->head->next; $pep; $pep = $pep->next) {
	my $name = $names{$pep};
	if($pep->next and $pep->next->start <= $pep->end + 1) {
	    $g->add_edge($name => $names{$pep->next},
			 dir => $direction);
	}
	foreach my $string (values %{$pep->right_strings}) {
	    $g->add_edge($name => $names{$string->right_pep},
			 label => ($string->right_pep->start - $pep->end - 1).' '.
				  $string->score,
			 dir => $direction);
	}
    }

    # Display graph
    if($self->{show_pep_graph})
	{ $self->_display_png($g->as_png); }
    if($self->{write_pep_graph}) {
	#my ($fh) = tempfile("graph_XXXX");
	#print $fh $g->as_png;
	open GRAPH_OUT, ">pep_graph_${nr}_$strand.png";
	print GRAPH_OUT $g->as_png;
	close GRAPH_OUT;
    }
#    my $mw = MainWindow->new;
#    my $gv = $mw->Scrolled ( 'GraphViz',
#                              -background => 'white',
#                              -scrollbars => 'sw' )
#             ->pack ( -expand => '1', -fill => 'both' );
#    $gv->show($g);
#    MainLoop;
}


sub _display_png
{
    my ($self, $png) = @_;
    my ($fh, $fn) = tempfile(SUFFIX => '.png');
    print $fh $png;
    close $fh;
    my $cmd = "perl -e \"system(\'kview $fn\');unlink('$fn');\" &";
    system($cmd);
}


#sub _display_png
#{
#    my ($self, $png) = @_;
#    my $pid;
#    if($pid = open (CHILD, "|-")) {
#	CHILD->autoflush(1);
#	print CHILD $png;
#    }
#    else {
#	die "cannot fork: $!" unless defined($pid);
#	system("display -");
#	exit;
#    }
#}


1;
