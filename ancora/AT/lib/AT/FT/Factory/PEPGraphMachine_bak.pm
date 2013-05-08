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

package AT::FT::Factory::PEPGraphMachine;

use strict;
use vars '@ISA';
use Carp;
use Bio::Graphics::Panel;
use AT::Root;
use AT::FT::Factory::MappingEndMover;
use AT::FT::Factory::MappingMerger;
use AT::FT::Factory::GFMappingMachine;
use AT::FT::Factory::MappingCompass;
use AT::Prediction::Gene;
use AT::Prediction::Exon;
use AT::Prediction::Intron;
use AT::Prediction::SpliceForm;
use GraphViz;
use Tk;
use Tk::GraphViz;
use File::Temp 'tempfile';


use Class::Struct '_map_element' => [
    start => '$',
    end => '$',
    startType => '$',
    endType => '$'
];

use Class::Struct '_map_exon' => [
    start => '$',
    end => '$',
    next => '_map_exon',
    set => '_map_exon_set',
    open_pep => '_pep',
    close_pep => '_pep'
];

use Class::Struct '_map_exon_set' => [
    mapping => 'AT::Mapping',
    grouped => '$',
    exons => '@',
    pep_score => '$'
];

use Class::Struct '_pep_string' => [
    left_pep => '_pep',
    right_pep => '_pep',
    score => '$'
];

use Class::Struct '_pep' => [
    start => '$',
    end => '$',
    score => '$',
    prev => '_pep',
    next => '_pep',
    staple_score => '$',      # currently not used?
    left_strings => '%',
    right_strings => '%',
    open_support => '@',
    close_support => '@',
    visited => '$'
];

use Class::Struct '_pos' => [
    pos => '$',
    open => '$',
    score => '$',
    exon => '_map_exon'
];


use constant MIN_INTRON_LEN_FOR_SS_CHECK => scalar 12;


@ISA = qw/AT::Root/;


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
	_exclude_on_revcom_ss => ($args{'exclude_on_revcom_ss'} || 0),
	show_chr_panel => ($args{'show_chr_panel'} || 0),
	show_pep_graph => ($args{'show_pep_graph'} || 0),
	write_pep_graph => ($args{'write_pep_graph'} || 0),
	gf_mapping_machine => AT::FT::Factory::GFMappingMachine->new(),
	mapping_merger => AT::FT::Factory::MappingMerger->new(),
	mapping_end_mover => AT::FT::Factory::MappingEndMover->new()
    }, ref $caller || $caller;
    
    return $self;
}


sub _init_chr_panel
{
    my ($self, $seq) = @_;

    unless($self->{show_chr_panel}) {
	$self->{_chr_panel} = 0;
	return;
    }

    $self->{_chr_panel} = Bio::Graphics::Panel->new
         ( -length => $seq->end - $seq->start + 1,
	   -offset => $seq->start-1,
	   -width => 800,
	   -grid => 1,
	   -pad_left => 10,
	   -pad_right => 10,
	   -pad_top => 10,
	   -pad_bottom => 10,
	   -key_style => 'between'
    );
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
    
    $self->_display_png($panel->png);
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
	$self->_draw_mappings_on_panel
	    ("strand ".($strand==1?'+':'-')." group ".(@$groups-$i),
	     [ map { $_->mapping } @{$groups->[$i]} ]);
    }
}


sub _draw_gf_mappings_on_panel
{
    my ($self, $key, @m) = @_;
    my $panel = $self->{_chr_panel} or return;

    my @feats = map { $self->_gf_mapping_to_SeqFeature($_) } @m;

    $panel->unshift_track
	    (\@feats,
	     -label => 1,
	     -glyph => 'segments', #sub { shift->type; },
	     -connector => sub { (shift->each_tag_value('connector'))[0] },
	     -key => $key
	    );
}


sub _gf_mapping_to_SeqFeature
{
    my ($self, $m) = @_;
    my @unbroken; my $i = 0;
    foreach my $part ($m->HSP_list) {
	push @{$unbroken[$i]}, $part;
	$i++ if($part->right_gap and
		$part->right_gap->is_broken);
    }
    my $group = Bio::Graphics::Feature->new
	(-type => 'segments',
	 -name => join(",", map {$_->qName} $m->primary_mapping_list),
	 -attributes => {connector => 'dashed'});
    foreach my $unbroken_set (@unbroken) {
	$group->add_segment(Bio::Graphics::Feature->new
	    (-type => 'segments',
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

    for(my $i = 0; $i < @data; $i+=2) {
	push @track_args, { key => $data[$i] };
	push @tracks, [ map { $self->_gene_to_SeqFeature($_) }
		        @{$data[$i+1]} ];
    }

    foreach my $i (0..@tracks-1) {
	$panel->unshift_track
	    ('segments' => $tracks[$i],
	     -label => 1,
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


sub _gene_to_SeqFeature
{
    my ($self, $g) = @_;
    my $f = Bio::Graphics::Feature->new
	(-segments => [ map { [ $_->abs_start, $_->abs_end ] } $g->exon_list ],
	 -strand => $g->strand,
	 -type => ref($g),
	 -source_tag => 'AT',
	 #-seq_id => $m->tName,
	 -name => $g->loc_str);
    return $f;
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
    my $strand = $args{strand} || 0;

    # init chr panel
    $self->_init_chr_panel($target_seq);

    # create and orient mappings
    my $gf_mappings = $self->gf_mapping_machine->run(mappings => $mappings,
	      				             target_seq => $target_seq);
    $self->mapping_end_mover->run(mappings => $gf_mappings, target_seq => $target_seq);
    $gf_mappings = $self->mapping_merger->run(mappings => $gf_mappings);
    my $orienter = AT::FT::Factory::MappingCompass->new;
    my ($plus_mappings, $minus_mappings, $excluded_mappings) =
	$orienter->run(mappings => $gf_mappings,
		       target_seq => $target_seq);

    # show mappings on panel
    $self->_draw_gf_mappings_on_panel
	('Excluded mappings', @$excluded_mappings);
    $self->_draw_gf_mappings_on_panel
	('Minus strand mappings', @$minus_mappings);
    $self->_draw_gf_mappings_on_panel
	('Plus strand mappings', @$plus_mappings);
    $self->_show_chr_panel();

    return @$gf_mappings;

    # create genes on plus strand
    my $plus_genes = $self->_mappings2genes
	($plus_mappings, 1, $target_seq, $excluded_mappings);
    # create genes on minus strand
    my $minus_genes = $self->_mappings2genes
	($minus_mappings, -1, $target_seq, $excluded_mappings);
    # show genes on panel
    $self->_draw_genes_on_panel("+ genes", $plus_genes, "- genes", $minus_genes);

    $self->_show_chr_panel();

    return (@$plus_genes, @$minus_genes);
}


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


sub _mappings2genes
{
    my ($self, $mappings, $strand, $target_seq, $excluded_mappings) = @_;

    my ($singles, $pep_graphs) = $self->_mappings2pepgraphs($mappings, $strand);  
    my $genes = $self->_pepgraphs2genes($pep_graphs, $target_seq, $strand,
					$excluded_mappings);
    
    free_peps($singles, $pep_graphs);

    return $genes;
}


# In the end, this method shouldn't be in this module
# one module should group the mappings
# another module should make the pep_graphs (or combine these first 2)
# another module should then infer the gene structures
# yet another module should then infer splice variants etc
# summarize all this in some module, e.g. M2GProcess or Factory::Process::M2G
# (we could then put the machines in Factory::Machine::...)

sub _pepgraphs2genes
{
    my ($self, $pep_graphs, $target_seq, $strand, $excl) = @_;
    my @genes;
    foreach my $pep_graph (@$pep_graphs) {
	for(my $pep = $pep_graph->next; $pep; $pep = $pep->next) {
	    next if($pep->visited);
	    my $peps = $self->_visit_connected_peps($pep);
	    push @genes, $self->_peps2gene($peps, $target_seq, $strand, $excl);
    	}
    }
    return \@genes;
}


sub _visit_connected_peps
{
    my ($self, $pep) = @_;

    my $pep_score_thr = 2;
    my $string_score_thr = 2;

    $pep->visited(1);
    return [] if($pep->score < $pep_score_thr);

    print STDERR "---\n";
    my @q = ($pep);
    for(my $i = 0; $i < @q; $i++) {
	my $pep = $q[$i];
	print STDERR "pep: ",$pep->start,"-",$pep->end," ",$pep->score,"\n";
	my $next = $pep->next;
	if($next and !($next->visited) and $next->start == $pep->end+1) {
	    print STDERR " next: ", $next->start, " ", ($next->visited || 0), "\n";
	    $next->visited(1);
	    if($next->score >= $pep_score_thr) {
		push @q, $next;
	    }
	}
	my $prev = $pep->prev;
	if($prev and !($prev->visited) and $prev->end == $pep->start-1) {
	    print STDERR " prev: ", $prev->end, " ", ($prev->visited || 0), "\n";
	    $prev->visited(1);
	    if($prev->score >= $pep_score_thr) {
		push @q, $prev;
	    }
	}
	foreach my $string (values %{$pep->left_strings}) {
	    print STDERR " left string: ",
		$string->right_pep->start, "->", $string->left_pep->end, " ",
		$string->score," ",
		($string->left_pep->visited || 0),"\n";
	    if(!($string->left_pep->visited) and
	       $string->score >= $string_score_thr) {
		$string->left_pep->visited(1);
		if($string->left_pep->score >= $pep_score_thr) {
		    push @q, $string->left_pep;
		    print STDERR "    followed string\n";
		}
	    }
	}
	foreach my $string (values %{$pep->right_strings}) {
	    print STDERR " right string: ",
		$string->left_pep->end, "->", $string->right_pep->start, " ",
		$string->score, " ",
		($string->right_pep->visited || 0),"\n";
	    if(!($string->right_pep->visited) and
	       $string->score >= $string_score_thr) {
		$string->right_pep->visited(1);
		if($string->right_pep->score >= $pep_score_thr) {
		    push @q, $string->right_pep;
		    print STDERR "    followed string\n";
		}
	    }
	}

    }

    print STDERR "---\n";

    return \@q;
}


sub _peps2gene {
    my ($self, $peps_ref, $target_seq, $strand, $excl) = @_;

    my @peps = sort {$a->start <=> $b->start} @$peps_ref;
    return unless (@peps);

    my $gene_start = $peps[0]->start;
    my $gene_end = $peps[-1]->end;
    my $seqstr = $target_seq->subseq($gene_start - $target_seq->start + 1,
				    $gene_end - $target_seq->start + 1);
    my $seq = Bio::LocatableSeq->new(-id => $target_seq->id,
				    -start => $gene_start,
				    -end => $gene_end,
				    -strand => $strand,
				    -seq => $seqstr);
    my $gene = AT::Prediction::Gene->new
	(excluded_mapping_list => $excl,
	 seq => $seq);

    my (@exons, @introns);
    my $i = 0;
    while ($i < @peps) {
	my $j = $i+1;
	while($j < @peps and $peps[$j]->start <= $peps[$j-1]->end + 1) {
	    $j++;
	}
	# make exon from peps $i..$j-1
	push @exons, AT::Prediction::Exon->new
	    ( -start => $peps[$i]->start - $gene_start + 1,
	      -end => $peps[$j-1]->end - $gene_start + 1,
	      -strand => $strand,
	      -abs_start => $peps[$i]->start,
	      -abs_end => $peps[$j-1]->end
	      );
	if($j < @peps) {
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

    print STDERR "Created gene: ", $gene->loc_str, "\n";
    $gene->print_exons_introns(\*STDERR);
    print "\n";

    return $gene;
}


sub _mappings2pepgraphs
{
    my ($self, $mappings, $strand) = @_;

    # Create repr of mappings as '_map_exon_sets' (or 'singles')
    my $min_gap = 12;
    #my (@spliced, @unspliced);
    my @singles;
    foreach my $mapping (@$mappings) {
	my @hsps = $mapping->all_HSPs;
	my $set = _map_exon_set->new(mapping => $mapping, grouped => 0);
	my $exon_head = _map_exon->new();
	my $prev_exon = $exon_head;
	my @exons;
	for(my ($j, $k) = (0, 0); $k < @hsps; $k++) {
	    if($k == (@hsps-1) or
	       $hsps[$k+1]->tStart >= $hsps[$k]->tEnd + $min_gap) {
		my $exon = _map_exon->new(start => $hsps[$j]->tStart,
					  end => $hsps[$k]->tEnd,
					  set => $set);
		$prev_exon->next($exon);
		$prev_exon = $exon;
		push @exons, $exon;
		$j = $k+1;
	    }
	}
	$set->exons(\@exons);
	#if(@exons > 1) { push @spliced, $set; }
	#else { push @unspliced, $set; }
	push @singles, $set;
    }

    my $groups = $self->_group_singles(\@singles, $strand);
    $self->_draw_map_groups_on_panel($groups, $strand);
    #$self->_draw_map_groups_on_panel([\@unspliced], $strand);

    my @pep_graphs = map { $self->_make_pep_graph($_) } @$groups;
    foreach my $pep_graph (@pep_graphs) {
	$self->_display_pep_graph($pep_graph, $strand);
    }

    return (\@singles, \@pep_graphs);
}


sub free_peps
{
    my ($exon_sets, $pep_graphs) = @_;

    my $dummy_pep = _pep->new;
    my @dummy_exons;
    my %dummy_hash;
    foreach my $es (@$exon_sets) {
	foreach my $e (@{$es->exons}) {
	    $e->open_pep($dummy_pep);
	    $e->close_pep($dummy_pep);
	}
	$es->exons(\@dummy_exons);
    }
    foreach my $pep_graph (@$pep_graphs) {
	for(my $pep = $pep_graph->next; $pep; $pep = $pep->next) {
	    $pep->prev($dummy_pep);
	    $pep->left_strings(\%dummy_hash);
	    $pep->right_strings(\%dummy_hash);
    	}
    }
}


sub _make_pep_graph
{
    my ($self, $group) = @_;

    # Create '_pos' objects that will help us create the graph
    # Each object is a start or end position of an '_exon' object
    # (the open attribute is true for starts, false for ends)
    # We sort the positions in increasing order
    my @pos =
	sort {$a->pos <=> $b->pos or $b->open <=> $a->open}
	(map { _pos->new(pos => $_->start, exon => $_, open => 1,
			 score => 1), #$_->set->pep_score),
	       _pos->new(pos => $_->end, exon => $_, open => 0,
			 score => -1)} #-$_->set->pep_score) }
	(map { @{$_->exons} } @$group));

    # Create linked list of peps by looking at each pair of adjacent positions
    # Identical positions are treated together.
    # i..j-1 is a range of identical positions in @pos to be used for a pep start
    # j..k-1 is a range of identical positions in @pos to be used for a pep end
    my $pep_head = _pep->new(start => 0, end => 0, score => 0);
	# ^^ dummy head for linked list
    my $prev_pep = $pep_head;
    my @exon_pep_assoc;
    my $score_change;
    my $score = $pos[0]->score;
    my ($i, $j) = (0,1);
    for(; $j < @pos and $pos[$i]->pos == $pos[$j]->pos; $j++) {
	    $score += $pos[$j]->score;
    }
    while($j < @pos) {

	# find next range; calc new score change
	$score_change = $pos[$j]->score;
	my $k = $j+1;
	for(; $k < @pos and $pos[$j]->pos == $pos[$k]->pos; $k++) {
	    $score_change += $pos[$k]->score;
	}

	if($score) {

	    # create pep for this region
	    my $pep = _pep->new(start => $pos[$i]->pos + ($pos[$i]->open ? 0 : 1),
				end => $pos[$j]->pos - ($pos[$j]->open ? 1 : 0),
				score => $score);

	    # here: need to tell peps which exons support them

	    # tell peps which exons they close
	    #$pep->open_support([map {$_->exon if($_->open)} @pos[$i..$j-1]]);
	    my @close_support = map {$_->open ? () : $_->exon} @pos[$j..$k-1];
	    $pep->close_support(\@close_support);
	    #  ^^^ beware of circularities

	    # tell exons which pep opens them
	    foreach my $pos (@pos[$i..$j-1]) {
		$pos->exon->open_pep($pep) if ($pos->open);
	    }
	    #foreach my $pos (@pos[$j..$k-1]) {
		#$pos->exon->close_pep($pep) unless($pos->open);
	    #}

	    # add pep to doubly linked list
	    $pep->prev($prev_pep);
	    $prev_pep->next($pep);
	    $prev_pep = $pep;
	}

	# move range fwd and update score
	($i, $j) = ($j, $k);
	$score += $score_change;
    }

    # Now add 'strings' (intronic connections) between peps
    #  for each pep: collect all supporting left exons
    #  they will give 0-more right exons
    #  each right exon will give one right pep
    #  for each left-right connection: increase the score
    for(my $left_pep = $pep_head; $left_pep; $left_pep = $left_pep->next) {
	foreach my $left_exon (@{$left_pep->close_support}) {
	    my $right_exon = $left_exon->next;
	    next unless ($right_exon);
	    my $score = $self->_string_score($left_exon, $right_exon);
	    my $right_pep = $right_exon->open_pep;
	    my $string = $left_pep->right_strings->{$right_pep};
	    unless($string) {
		$left_pep->right_strings->{$right_pep} =
		$right_pep->left_strings->{$left_pep} =
		    _pep_string->new(left_pep => $left_pep,
				     right_pep => $right_pep,
				     score => $score); 
	    }
	    else {
		$string->score($string->score + $score);
	    }
	}
    }

    return $pep_head;
}


sub _string_score { 1; }


sub _display_pep_graph
{
    my ($self, $pepg, $strand) = @_;
    return unless ($self->{show_pep_graph} or $self->{write_pep_graph});

    my $g = GraphViz->new(); # Create GraphViz object
    my $direction = $strand == 1 ? 'forward' : 'back';

    # Assign a name to each pep and add nodes
    my %names;
    my $major = ord('A');
    my $minor = 1;
    for(my $pep = $pepg->next; $pep; $pep = $pep->next) {
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
    for(my $pep = $pepg->next; $pep; $pep = $pep->next) {
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
	my ($fh) = tempfile("graph_XXXX");
	print $fh $g->as_png;
    }
#    my $mw = MainWindow->new;
#    my $gv = $mw->Scrolled ( 'GraphViz',
#                              -background => 'white',
#                              -scrollbars => 'sw' )
#             ->pack ( -expand => '1', -fill => 'both' );
#    $gv->show($g);
#    MainLoop;
}


sub _group_singles {
    my ($self, $singles, $strand) = @_;

# Cluster compatible singles into 'groups'

# make a graph connecting singles that share at least 1/2 of their spl jnc
# must compare all ag all to do this
# then to bfs to get groups of connected singles

# seed group with a single, try to include others
# include if jnc shared w at leas 1 group member
# iterate until no more can be included

    # Sort in incr exon order to make sure that seed groups
    # with unspliced before seeding groups with spliced.
    # It is probably also an optimization to seed groups with
    # mappings with a large number of exons.
    $singles = [ sort { @{$a->exons} <=> @{$b->exons} } @$singles ];

    if(1) {
	print STDERR "\nSingles: ";
	foreach my $i (0..@$singles-1) {
	    print STDERR $i.'.'.$singles->[$i]->mapping->qName.', ';
	}
	print STDERR "\n\n";
    }

    my @group_lists;
    my $ungrouped = @$singles-1; # index of some single yet ungrouped
    while($ungrouped != -1) {	# loop while there are ungrouped singles
	#print STDERR "group seed: $ungrouped\n";
	my @group_flg = (0) x scalar(@$singles);
	my @group_lst = ($singles->[$ungrouped]);
	$group_flg[$ungrouped] = 1; # start group with the ungrouped
	$singles->[$ungrouped]->grouped(1);
	my $grouped_one = 1;
	while($grouped_one and $ungrouped != -1) {
	    $grouped_one = 0;
	    $ungrouped = -1;    # assume we are out of singles to group
	    for(my $i = 0; $i < @$singles; $i++) {
		if($group_flg[$i]) {
		    # already in group - do nothing
		    #print STDERR "$i: present\n";
		}
		elsif(($singles->[$i]->grouped == 0
			#or @{$singles->[$i]->exons} == 1
			) and
		      $self->_fits_in_group($singles->[$i], \@group_lst, $strand)) {
		    # single fits in group
		    $group_flg[$i] = 1;
		    push @group_lst, $singles->[$i];
		    $singles->[$i]->grouped(1);
		    $grouped_one = 1;
		    #print STDERR "$i: added\n";
		}
		elsif($singles->[$i]->grouped == 0) {
		    # single does not fit and is ungrouped
		    $ungrouped = $i;
		    #print STDERR "$i: ungrouped\n";
		}
	    }
	}
	#print STDERR (join ' ', @group_flg), " ",
	#    (join ",", map {$_->mapping->qName} @group_lst),"\n";
	push @group_lists, \@group_lst;
    }

    return \@group_lists;
}


sub _fits_in_group {
    my ($self, $candidate, $group, $strand) = @_;
    foreach my $member (@$group) {
	if($self->_singles_compatible($candidate, $member, $strand)) {
	    return 1;
	}
    }
    return 0;
}


#PARAMETERS:
#sing_sing_minOl_bp
#sing_sing_minOl_pc
#sing_mul_maxJncGlitch_bp
#sing_mulExt_minOl_bp
#sing_mulExt_minOl_pc
#mul_mul_maxJncGlitch_bp
#mul_mul_minJncMatches_abs
#mul_mul_minJncMatches_pc

sub _singles_compatible {
    my ($self, $a, $b, $strand) = @_;
    my $spljnc_glitch = 0;
    my $conf_glitch = 10;
    my $exa = $a->exons;
    my $exb = $b->exons;

    # note: a and b are treated differently
    # a should be the candidate, b the member

    # return false if they don't overlap at all
    if($exa->[0]->start > $exb->[-1]->end or $exa->[-1]->end < $exb->[0]->start)
    { return 0 };

    if(@$exa == 1 and @$exb == 1) {
	# both single-exon
	# compatible if overlap is at least 50 bp or 50% of the smaller
	return $self->_overlap_size_check($exa->[0]->start, $exa->[0]->end,
					  $exb->[0]->start, $exb->[0]->end);
    }
    elsif(@$exa == 1) {
	# a is single-exon, b is multi-exon
	my($s1, $e1) = ($exa->[0]->start, $exa->[0]->end);
	if($e1 <= $exb->[0]->end+$conf_glitch) {
	    return $self->_overlap_size_check($s1, $e1,
					      $exb->[0]->start, $exb->[0]->end);
	}
	elsif($s1 >= $exb->[-1]->start-$conf_glitch) {
	    return $self->_overlap_size_check($s1, $e1,
					      $exb->[-1]->start, $exb->[-1]->end);
	}
	else {
	    foreach my $i (1 .. @$exb-1) {
		if($s1 >= $exb->[$i]->start-$conf_glitch and
		    $e1 <= $exb->[$i]->end+$conf_glitch) {
		    return 1;
		}
	    }
	}
    }
    elsif(@$exb == 1) {
	# a is multi-exon, b is single-exon
	return 0;
    }
    else {
	# both multi-exon
	#print STDERR "comp-check: ",$a->{mapping}->qName," ", $b->{mapping}->qName, "\n";
	my $shared_spljnc = 0;
	my($i, $j) = (0,0);
	while($i < @$exa and $j < @$exb) {
	    my($s1, $e1) = ($exa->[$i]->start, $exa->[$i]->end);
	    my($s2, $e2) = ($exb->[$j]->start, $exb->[$j]->end);
	    #print STDERR $i."[$s1-$e1]:".$j."[$s2-$e2]  ";
	    if ($i and $j and
		$s1 >= $s2 - $spljnc_glitch and $s1 <= $s2 + $spljnc_glitch) {
		$shared_spljnc++;
	    }
	    if($i < @$exa-1 and $j < @$exb-1 and
		$e1 >= $e2 - $spljnc_glitch and $e1 <= $e2 + $spljnc_glitch) {
		$shared_spljnc++;
	    }
	    if ($e1 < $e2) {
		$i++;
	    }
	    else {
		$j++;
	    }
	}
	my $min_exons = (@$exa < @$exb) ? @$exa : @$exb;
	my $min_spljnc = $min_exons * 2 - 2;
	#print STDERR "\n$shared_spljnc $min_spljnc\n";
	return (2*$shared_spljnc >= $min_spljnc) ? 1 : 0;
    }
}


sub _overlap_size_check
{
    my ($self, $s1, $e1, $s2, $e2) = @_;
    my $ol_start = ($s1 > $s2) ? $s1 : $s2;
    my $ol_end = ($e1 < $e2) ? $e1 : $e2;
    my $ol = $ol_end - $ol_start;
    return 1 if ($ol >= 50);
    my $size1 = $e1 - $s1 + 1;
    my $size2 = $e2 - $s2 + 1;
    my $min_size = ($size1 < $size2) ? $size1 : $size2;
    return (2*$ol >= $min_size) ? 1 : 0;
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


sub _mapping_category
{
    my ($self, $m) = @_;
    if($m->qType eq 'mRNA' and $m->qName =~ /^NM_/) {
	if($m->refSeqStatus eq 'Reviewed') {
	    return "refSeqRev";
	}
	else {
	    return "refSeq";
	}
    }
    else {
	return $m->qType;
    }
}


sub _display_png
{
    my ($self, $png) = @_;
    my $pid;
    if($pid = open (CHILD, "|-")) {
	CHILD->autoflush(1);
	print CHILD $png;
    }
    else {
	die "cannot fork: $!" unless defined($pid);
	system("display -");
	exit;
    }
}


1;
