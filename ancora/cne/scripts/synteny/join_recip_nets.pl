#!/usr/bin/perl
use warnings;
use strict;
use DBI;
use Sys::Hostname;
use Class::Struct;
use Getopt::Std;
use AT::Tools::RangeHandler;

struct MyChain => [start1 => '$',
		   end1 => '$',
		   chr2 => '$',
		   start2 => '$',
		   end2 => '$',
		   strand => '$',
		   ali => '$',
		   blocks1 => '@',
		   blocks2 => '@'];

struct MyChainCluster => [start1 => '$',
			  end1 => '$',
			  chr2 => '$',
			  start2 => '$',
			  end2 => '$',
			  ali => '$',
			  strand => '$',
			  chains => '@',
			  retained => '$'];

struct vertex => [ id => '$', visited => '$', edges => '%' ];
struct edge => [ labels => '%', weight => '$', visited => '$' ];

my $DB_HOST = (hostname() =~ /^angband\./) ? 'localhost' : 'asp';
my $DB_PORT = (hostname() =~ /^angband\./) ? 3306 : 3307;
my $FILTER = 0;

my %args;
getopts('o:', \%args);
my $BASE_FN = $args{'o'};
my ($ASM1, $ASM2, $CHAIN_FN, $ASM1_SEQID_FN, $ASM2_SEQID_FN, $MAX_GAP1, $MAX_GAP2) = @ARGV;
$MAX_GAP2 = $MAX_GAP1 unless(defined $MAX_GAP2);
unless(defined $MAX_GAP1) {
    die "usage: perl $0 [-o <outfile-basename>] <asm1> <asm2> <chain-file> <asm1.chr.list> <asm2.chr.list> <max-gap1> [max-gap2]\n";
}
my $TRACK_NAME1 = "$ASM2 joinedRecipNets";
my $TRACK_NAME2 = "$ASM1 joinedRecipNets";
my $NET_TABLE = 'net'.ucfirst($ASM2);
my $max_gap1_str = int($MAX_GAP1/1000)."k";
my $max_gap2_str = int($MAX_GAP2/1000)."k";
$BASE_FN = "joinedRecipNets_${ASM1}_${ASM2}_gap${max_gap1_str}${max_gap2_str}" unless($BASE_FN);
my $track_header1 = "track name=\"${TRACK_NAME1}_${max_gap1_str}_${max_gap2_str}\" ".
    "description=\"Reciprocal synteny from $CHAIN_FN, maxgap=($max_gap1_str,$max_gap2_str)\" ".
    "useScore=1 visibility=3\n";
my $track_header2 = "track name=\"${TRACK_NAME2}_${max_gap1_str}_${max_gap2_str}\" ".
    "description=\"Reciprocal synteny from $CHAIN_FN, maxgap=($max_gap1_str,$max_gap2_str)\" ".
    "useScore=1 visibility=3\n";

my $chr_list1 = read_seq_id_list($ASM1_SEQID_FN);
my $chr_list2 = read_seq_id_list($ASM2_SEQID_FN);

open BED1, ">$BASE_FN.$ASM1.bed" or die "open failed";
open BED2, ">$BASE_FN.$ASM2.bed" or die "open failed";
open DETAIL1, ">$BASE_FN.${ASM1}_detail.bed" or die "open failed";
open DETAIL2, ">$BASE_FN.${ASM2}_detail.bed" or die "open failed";
open SUMMARY, ">$BASE_FN.txt" or die "open failed";
print BED1 $track_header1;
print BED2 $track_header2;
print DETAIL1 $track_header1;
print DETAIL2 $track_header2;

foreach my $chr1 (@$chr_list1) {
    my @clusters;
    foreach my $chr2 (@$chr_list2) {
	print STDERR "joining nets for $chr1, $chr2...\n";
	my @chains;
	my $g = graph__new();
	my $i = 0;
	my $chain_filter_cmd = "chainFilter -t=$chr1 -q=$chr2 $CHAIN_FN |";
	open CHAIN, $chain_filter_cmd or die "could not execute: $chain_filter_cmd";
	while(my $chain = read_next_chain(\*CHAIN)) {
	    $chains[$i] = $chain;
	    graph__add_vertex($g, $i);
	    for my $j (0..$i-1) {
		graph__add_edge($g, $i, $j) if chains_are_close($chains[$j], $chains[$i]);
	    }
	    $i++;
	}
	my $components = graph__connected_components($g);
	foreach my $component (@$components) {
	    my $vertices = $component->[0];
	    next unless(@$vertices);
	    my @indexes = map {$_->id} @$vertices;
	    push @clusters, merge_chains(@chains[@indexes]);
	}
    }
    if($FILTER) {
	print STDERR "filtering joined nets for $chr1...\n";
	filter_chain_clusters(\@clusters);
    }
    print STDERR "printing result for $chr1...\n";
    print_chain_clusters(\*SUMMARY, \*BED1, \*BED2, \*DETAIL1, \*DETAIL2, $chr1, \@clusters);
}

close BED1;
close BED2;
close SUMMARY;


sub read_next_chain
{
    my $fh = shift;

    # Read chain header
    my $header = <$fh>;
    return unless($header);
    chomp $header;
    my ($type, undef, undef, undef, undef, $tStart, $tEnd, $qName, $qSize, $qStrand, $qStart, $qEnd) = split /\s+/, $header;

    # Read chain blocks
    my (@blocks1, @blocks2, $ali, $spacer);
    if($qStrand eq '+') {
	my $pos1 = $tStart+1;
	my $pos2 = $qStart+1;
	while(my $line = <$fh>) {
	    chomp $line;
	    my ($size, $dt, $dq) = split /\t/, $line;
	    $ali += $size;
	    push @blocks1, [$pos1, $pos1+$size-1];
	    push @blocks2, [$pos2, $pos2+$size-1];
	    if(defined $dt) {
		$pos1 += $size + $dt;
		$pos2 += $size + $dq;
	    }
	    else {
		$spacer = <$fh>;
		last;
	    }
	}
    }
    elsif($qStrand eq '-') {
	($qStart, $qEnd) = ($qSize-$qEnd, $qSize-$qStart);
	my $pos1 = $tStart+1;
	my $pos2 = $qEnd;
	while(my $line = <$fh>) {
	    chomp $line;
	    my ($size, $dt, $dq) = split /\t/, $line;
	    $ali += $size;
	    push @blocks1, [$pos1, $pos1+$size-1];
	    push @blocks2, [$pos2-$size+1, $pos2];
	    if(defined $dt) {
		$pos1 += $size + $dt;
		$pos2 -= $size + $dq;
	    }
	    else {
		$spacer = <$fh>;
		last;
	    }
	}
	@blocks2 = reverse @blocks2;
    }
    else {
	die "unknown strand $qStrand in chain file";
    }

    # Check that chain record ended with a blank line
    die "early EOF in chain file" unless(defined $spacer and $spacer =~ /^\s+$/);

    # Consistency check
    unless($blocks1[0][0] eq $tStart+1 and
	   $blocks2[0][0] eq $qStart+1 and
	   $blocks1[-1][1] eq $tEnd and
	   $blocks2[-1][1] eq $qEnd) {
	die "computed block start and/or end does not match header in chain file";
    }

    # Create chain object
    my $chain = MyChain->new(chr2 => $qName,
			     start1 => $blocks1[0][0],
			     start2 => $blocks2[0][0],
			     end1 => $blocks1[-1][1],
			     end2 => $blocks2[-1][1],
			     strand => $qStrand,
			     ali => $ali,
			     blocks1 => \@blocks1,
			     blocks2 => \@blocks2
			     );

    return $chain;
}


sub print_chain_clusters
{
    my ($summary_fh, $bed_fh1, $bed_fh2, $detail_fh1, $detail_fh2, $chr1, $clusters) = @_;
    foreach my $cluster (@$clusters) {
	my $chr2 = $cluster->chr2;
	my $start1 = $cluster->start1;
	my $start2 = $cluster->start2;
	my $name1 = $chr2 . '_' . int($start2/1000) . 'k';
	my $name2 = $chr1 . '_' . int($start1/1000) . 'k';
	my $score = $cluster->retained ? 1000 : 300;
	my $nets = $cluster->chains();
	my $strand = $cluster->strand;
	my (@blockSizes1, @blockStarts1, @blockSizes2, @blockStarts2);
	my ($block_starts1, $block_sizes1, $block_detail_starts1, $block_detail_sizes1) = collapse_nets($nets,$start1,1);
	my ($block_starts2, $block_sizes2, $block_detail_starts2, $block_detail_sizes2) = collapse_nets($nets,$start2,2);
	print $summary_fh join("\t",
			       $chr1, $start1, $cluster->end1,
			       $chr2, $start2, $cluster->end2,
			       scalar(@$nets),
			       $cluster->retained ? 'retained' : 'discarded',
			       $cluster->ali), "\n";
	print $bed_fh1 join("\t",
			    $chr1, $start1-1, $cluster->end1,
			    $name1, $score, $strand,
			    $start1-1, $cluster->end1, '0,0,0',
			    scalar(@$block_sizes1), join(',',@$block_sizes1), join(',',@$block_starts1)), "\n";
	print $detail_fh1 join("\t",
			       $chr1, $start1-1, $cluster->end1,
			       $name1, $score, $strand,
			       $start1-1, $cluster->end1, '0,0,0',
			       scalar(@$block_detail_sizes1), join(',',@$block_detail_sizes1), join(',',@$block_detail_starts1)), "\n";
	print $bed_fh2 join("\t",
			    $chr2, $cluster->start2-1, $cluster->end2,
			    $name2, $score, $strand,
			    $start2-1, $cluster->end2, '0,0,0',
			    scalar(@$block_sizes2), join(',',@$block_sizes2), join(',',@$block_starts2)), "\n";
	print $detail_fh2 join("\t",
			       $chr2, $cluster->start2-1, $cluster->end2,
			       $name2, $score, $strand,
			       $start2-1, $cluster->end2, '0,0,0',
			       scalar(@$block_detail_sizes2), join(',',@$block_detail_sizes2), join(',',@$block_detail_starts2)), "\n";
    }	
}


sub collapse_nets
{
    my ($nets,$offset,$index) = @_;

    my ($blocks, $blocks_detail) = ([],[]);
    
    if($index == 1) {
	foreach my $net (@$nets) {
	    push @$blocks, [$net->start1, $net->end1];
	    push @$blocks_detail, @{$net->blocks1};
	}
    }
    else {
	foreach my $net (@$nets) {
	    push @$blocks, [$net->start2, $net->end2];
	    push @$blocks_detail, @{$net->blocks2};
	}
    }

    $blocks = collapse_regions($blocks);
    $blocks_detail = collapse_regions($blocks_detail);

    my (@starts, @sizes, @starts_detail, @sizes_detail);
    foreach my $block (@$blocks) {
	push @starts, $block->[0]-$offset;
	push @sizes, $block->[1]-$block->[0]+1;
    }
    foreach my $block (@$blocks_detail) {
	push @starts_detail, $block->[0]-$offset;
	push @sizes_detail, $block->[1]-$block->[0]+1;
    }

    return (\@starts, \@sizes, \@starts_detail, \@sizes_detail);
}


sub collapse_regions
{
    my $blocks_ref = shift;

    my @blocks = sort { $a->[0] <=> $b->[0] } @$blocks_ref;
	
    my @iv;
    my ($start, $end) = @{$blocks[0]};
    for my $i (1..@blocks-1) {
        my ($my_start, $my_end) = @{$blocks[$i]};
        if($my_start > $end) {
            push @iv, [$start,$end];
            ($start, $end) = ($my_start, $my_end);
        }
        else {
            $end = $my_end if ($end < $my_end);
        }
    }
    push @iv, [$start,$end];
    return \@iv;
}


sub filter_chain_clusters {
    my ($clusters) = @_;
    my @sorted_clusters = sort { $b->ali <=> $a->ali } @$clusters;
    for my $i (0..@sorted_clusters-1) {
	my $c1 = $sorted_clusters[$i];
	my $overlaps;
	for my $j (0..$i-1) {
	    my $c2 = $sorted_clusters[$j];
	    if($c1->start1 <= $c2->end1 and $c1->end1 >= $c2->start1) {
		$overlaps = 1;
		last;
	    }
	}
	if($overlaps) {
	    $c1->retained(0);
	    # warn maybe
	}
	else {
	    $c1->retained(1);
	}
    }
}


sub merge_chains
{
    my (@net_list) = @_;
    my $nr_nets = scalar @net_list;
    my $first = $net_list[0];
    my $cluster = MyChainCluster->new(chr2 => $first->chr2,
				      start1 => $first->start1,
				      start2 => $first->start2,
				      end1 => $first->end1,
				      end2 => $first->end2,
				      ali => $first->ali,
				      strand => $first->strand,
				      chains => \@net_list,
				      retained => 1
				      );
    foreach my $net (@net_list[1..@net_list-1]) {
	$cluster->start1($net->start1) if($cluster->start1 > $net->start1);
	$cluster->start2($net->start2) if($cluster->start2 > $net->start2);
	$cluster->end1($net->end1) if($cluster->end1 < $net->end1);
	$cluster->end2($net->end2) if($cluster->end2 < $net->end2);
	$cluster->ali($cluster->ali + $net->ali);
	$cluster->strand('.') if($cluster->strand ne '.' and $cluster->strand ne $net->strand);
    }
    return $cluster;
}


sub chains_are_close {
    my ($net1, $net2) = @_;
    return 0 if(distance($net1->start1, $net1->end1, $net2->start1, $net2->end1) > $MAX_GAP1);
    return 0 if(distance($net1->start2, $net1->end2, $net2->start2, $net2->end2) > $MAX_GAP2);
    return 1;
}


sub distance {
    my ($s1, $e1, $s2, $e2) = @_;
    if($e1 < $s2) {
	return $s2 - $e1;
    }
    elsif($e2 < $s1) {
	return $s1 - $e2;
    }
    else {
	return 0;
    }
}


sub read_seq_id_list
{
    my $fn = shift;
    open IN, $fn or die "could not open $fn";
    my @list;
    while(my $line = <IN>) {
	chomp $line;
	next if($line =~ /^\#/ or $line =~ /^\s*$/);
	my ($id) = $line =~ /^\s*(\S+)/;
	push @list, $id;
    }
    return \@list;
}




##
## GRAPH subs
##


sub graph__new
{
    return {}; # just an emty hash of vertices
}


sub graph__add_vertices
{
    my $g = shift;
    foreach my $id (@_) {
        unless($g->{$id}) {
            $g->{$id} = vertex->new(id => $id, visited => 0);
        }
    }
    return $g;
}


sub graph__add_vertex
{
    my ($g, $id) = @_;
    my $v = $g->{$id};
    unless($v) {
        $g->{$id} = $v = vertex->new(id => $id, visited => 0);
    }
    return $v;
}


sub graph__add_edge
{
    my ($g, $id1, $id2, $weight, $label) = @_;
    my $v1 = $g->{$id1} || graph__add_vertex($g, $id1);
    my $v2 = $g->{$id2} || graph__add_vertex($g, $id2);
    my $e = $v1->edges->{$id2};
    unless($e) {
	$e = edge->new(visited => 0, weight => 0);
    }
    $e->labels->{$label} = 1 if($label);
    $e->weight($e->weight+$weight) if($weight);
    $v1->edges->{$id2} = $e;
    $v2->edges->{$id1} = $e;
}


sub graph__connected_components
{
    my $g = shift;
    my @components;
    # here: should zero all visited flags
    foreach my $v (values %$g) {
        next if($v->visited);
        push @components, [graph__dfs_visit($g,$v)];
    }
    return \@components;
}


sub graph__dfs_visit
{
    my ($g,$v1) = @_;
    my @vertices = ($v1);
    my @edges;
    $v1->visited(1);
    my $edges = $v1->edges;
    foreach my $v2_id (keys %$edges) {
	my $v2 = $g->{$v2_id};
	my $e = $edges->{$v2_id};
	unless($e->visited) {
	    $e->visited(1);
	    push @edges, $e;
	    unless($v2->visited) {
		my ($more_vertices, $more_edges) = graph__dfs_visit($g,$v2);
		push @vertices, @$more_vertices;
		push @edges, @$more_edges;
	    }
	}
    }

    return (\@vertices, \@edges);
}
