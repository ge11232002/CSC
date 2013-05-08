#! /usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use FileHandle;
use File::Basename;
use Class::Struct;
use AT::Tools::RangeHandler;
use AT::GFX::SeqColorIndex;

struct vertex => [ id => '$', visited => '$', edges => '%' ];
struct edge => [ labels => '%', weight => '$', visited => '$' ];
struct cne => [ 
    start1 => '$',
    end1 => '$',
    start2 => '$',
    end2 => '$',
    strand => '$',
    score => '$',
    cigar => '$'];

my $COLOR_INDEX = AT::GFX::SeqColorIndex->new();
my %args;
getopts('bs',\%args);
my $MAKE_BED = $args{'b'};
my $SILENT = $args{'s'};

main();
exit;

sub usage
{
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd 

Usage: perl $cmd [options] file1 file2

Options:
 -b Create bed files
 -s Silent (less progress information on stderr)

EndOfUsage
exit;
}


sub main
{
    # Process commandline args
    my ($file1, $file2) = @ARGV;
    usage() unless($file2);

    # Parse species names and percentage identity cutoff from file name
    my $win_suffix;
    my ($basename1) = fileparse($file1);
    my ($sp1, $sp2, $cutoff, $win_size) = $basename1 =~ /cne_(\w+)_(\w+)_(\d+)_(\d+)$/;
    die "Could not parse species and cutoff from filename $file1\n" unless($sp1 and $sp2 and $cutoff and $win_size);
    my $cutoff_pc = int(100*($cutoff/$win_size));
    my $cutoff_frac = $cutoff_pc / 100;

    # Read CNEs
    my ($cnes_by_chr, $chr_names1, $chr_names2) = load_cnes($file1, $file2);

    # Open output files
    my $cne_fn =  "cne2w_${sp1}_${sp2}_${cutoff}_$win_size";
    my $cne_out =  FileHandle->new(">$cne_fn") or die "could not open $cne_fn";
    my ($bed_out1, $bed_out2);
    if($MAKE_BED) {
	my $bed_fn1 = "cne2w_${sp1}_${sp2}_${cutoff_frac}_$win_size.bed";
	my $bed_fn2 = "cne2w_${sp2}_${sp1}_${cutoff_frac}_$win_size.bed";
	$bed_out1 = FileHandle->new(">$bed_fn1") or die "could not open $bed_fn1";
	$bed_out2 = FileHandle->new(">$bed_fn2") or die "could not open $bed_fn2";
	print $bed_out1 "track name=\"cne2w_${sp1}_${sp2}_${cutoff_frac}_$win_size\" description=\"CNEs 2-Way ${sp1}/${sp2} $cutoff_pc\% / $win_size cols in ${sp1} coordinates\" itemRgb=On\n";
	print $bed_out2 "track name=\"cne2w_${sp2}_${sp1}_${cutoff_frac}_$win_size\" description=\"CNEs 2-Way ${sp1}/${sp2} $cutoff_pc\% / $win_size cols in ${sp2} coordinates\" itemRgb=On\n";
    }

    # Merge CNEs and output
    print STDERR "Merging...\n" unless($SILENT);
    foreach my $chr1 (sort @$chr_names1) {
	foreach my $chr2 (sort @$chr_names2) {
	    
	    my $cnes = $cnes_by_chr->{$chr1}{$chr2};
	    next unless($cnes);
	    
	    $cnes = merge_cnes($cnes);
	    
	    # Output merged CNEs
	    foreach my $cne (@$cnes) {
		
		# Output CNE format
		print $cne_out join("\t",
				    $chr1, $cne->start1-1, $cne->end1, 
				    $chr2, $cne->start2-1, $cne->end2,
				    $cne->strand, $cne->score, $cne->cigar), "\n";
		
		# Output bed format
		if($MAKE_BED) {
		    print_bed($bed_out1,
			      $chr1, $cne->start1, $cne->end1,
			      $chr2, $cne->start2, $cne->end2,
			      $cne->strand, $cne->score);
		    print_bed($bed_out2, 
			      $chr2, $cne->start2, $cne->end2,
			      $chr1, $cne->start1, $cne->end1,
			      $cne->strand, $cne->score);
		}
	    }
	}
    }
}


# Print CNE in bed format
sub print_bed
{
    my ($fh, $chr1, $start1, $end1, $chr2, $start2, $end2, $strand, $score) = @_;
    $score = int($score * 10 + .5);
    my $color = join(',',$COLOR_INDEX->get_seq_color($chr2));
    print $fh join("\t",
		   $chr1, $start1-1, $end1,        # chr,start,end
		   $chr2.':'.$start2.'-'.$end2,    # name
		   $score, $strand,                # score, strand
		   $start1-1, $end1,               # thickStart,thickEnd
		   $color                          # color
	), "\n";
}


# Read two CNE files into memory
sub load_cnes
{
    my ($file1, $file2) = @_;

    my %cnes_by_chr;
    my (%chr_names1, %chr_names2);
    my %seen_cnes;
    my ($n1, $n2, $n2_unique) = (0,0,0);

    # Load File 1
    open(FILE1, $file1) or die "Could not open $file1\n";
    while (<FILE1>) {
	chomp;
	my ($chr1,$start1,$end1,$chr2,$start2,$end2,$strand,$score,$cigar) = split "\t";
	$n1++;
	$seen_cnes{"$chr1:$start1-$end1,$chr2:$start2-$end2"} = 1;
	my $cne = cne->new(start1 => $start1+1,
			   start2 => $start2+1,
			   end1 => $end1,
			   end2 => $end2,
			   strand => $strand,
			   score => $score,
			   cigar => $cigar);
	push @{$cnes_by_chr{$chr1}{$chr2}}, $cne;
	$chr_names1{$chr1} = 1;
	$chr_names2{$chr2} = 1;
    } 
    close FILE1;
    print STDERR "Loaded $n1 CNEs from file 1.\n" unless($SILENT);

    # Load File 1 - switch coordinates around: sp2 coordinates first!
    open FILE2, $file2 or die "could not open $file2";
    while(<FILE2>) {
	chomp;
	my ($chr2,$start2,$end2,$chr1,$start1,$end1,$strand,$score,$cigar) = split "\t";
	$n2++;
	next if($seen_cnes{"$chr1:$start1-$end1,$chr2:$start2-$end2"});
	$n2_unique++;
	#print STDERR "$chr1:$start1-$end1 $strand\n$cigar\n";
	$cigar =~ tr/DI/ID/;
	$cigar = reverse_cigar($cigar) if($strand eq '-');
	#print STDERR "$cigar\n\n";
	my $cne = cne->new(start1 => $start1+1,
			   start2 => $start2+1,
			   end1 => $end1,
			   end2 => $end2,
			   strand => $strand,
			   score => $score,
			   cigar => $cigar);
	push @{$cnes_by_chr{$chr1}{$chr2}}, $cne;
	$chr_names1{$chr1} = 1;
	$chr_names2{$chr2} = 1;
    }
    close FILE2;  
    print STDERR "Read $n2 CNEs from file 2. Loaded $n2_unique unique to file 2.\n" unless($SILENT);

    return (\%cnes_by_chr, [keys %chr_names1], [keys %chr_names2]);
}


# Reverse cigar string
sub reverse_cigar {
    my $fwd = shift;
    my $rev = join('', reverse $fwd =~ /(\d+[MDI])/g);
    die "Illegal cigar string $fwd" if(length($rev) != length($fwd));
    return $rev;
}


# Take a list of CNEs and return a nonredundant set of merged CNEs.
# This is the core function of this script.
sub merge_cnes {
    my ($cnes_ref) = @_;

    # Sort CNEs by start in assembly 1
    my @cnes = sort { $a->start1 <=> $b->start1 } @$cnes_ref;

    # Count CNEs
    my $nr_cnes = @cnes;
    return [] unless($nr_cnes); # return if nothing to do

    # Create a graph where CNEs that overlap in both genomes are linked
    my $g = graph__new();
    for my $i (0..$nr_cnes-1) {
	graph__add_vertex($g, $i);
	for my $j ($i+1..$nr_cnes-1) {
	    last if($cnes[$i]->end1+1 < $cnes[$j]->start1);
	    graph__add_edge($g, $i, $j)
		if($cnes[$i]->start2 <= $cnes[$j]->end2+1 and $cnes[$i]->end2+1 >= $cnes[$j]->start2);	    
	}
    }

    # Get the connected components of the graph
    my $components = graph__connected_components($g);

    # Process each connected component (CNE cluster)
    my @merged_cnes;
    foreach my $component (@$components) {

	# Get CNEs in component
	my $vertices = $component->[0];
	next unless(@$vertices);
	my @cne_cluster = @cnes[map {$_->id} @$vertices];
	my $cluster_size = scalar @cne_cluster;

	# Reduce the cluster to a set of nonredundant CNEs
	# This means removing CNEs that fall completely wihtin other CNEs
	for my $i (1..$cluster_size) {    
	    my $cne1 = shift @cne_cluster;
	    foreach my $cne2 (@cne_cluster) {	
		if($cne1->start1 >= $cne2->start1 and $cne1->end1 <= $cne2->end1 and 
		   $cne1->start2 >= $cne2->start2 and $cne1->end2 <= $cne2->end2) {
		    $cne1 = 0;
		    last;
		}
	    }
	    push @cne_cluster, $cne1 if($cne1);
	}

	# How many CNEs left in cluster after removing redundance?
	if(@cne_cluster == 1) {    
	    # One left: just keep that one - easy
	    push @merged_cnes, @cne_cluster;
	}
	else {
	    # Several left: merge them into one (this case is very rare)
	    push @merged_cnes, create_merged_cne(\@cne_cluster);
	}
    }

    return \@merged_cnes;
}


# Merge several CNEs into one
sub create_merged_cne
{
    my $cnes = shift;
    my ($start1, $end1, $start2, $end2, $strand);
    my $score_sum = 0;
    my @cigar_parts;
    foreach my $c (@$cnes) {
	$start1 = $c->start1 if(!defined($start1) or $start1 > $c->start1);
	$start2 = $c->start2 if(!defined($start2) or $start2 > $c->start2);
	$end1 = $c->end1 if(!defined($end1) or $end1 < $c->end1);
	$end2 = $c->end2 if(!defined($end2) or $end2 < $c->end2);
	if(!defined($strand)) {
	    $strand = $c->strand;
	}
	elsif($strand ne $c->strand) {
	    $strand = '.';
	}
	$score_sum += $c->score;
	push @cigar_parts, join(',', $c->start1.'-'.$c->end1, $c->start2.'-'.$c->end2, $c->strand, $c->cigar);
    }
    my $score = $score_sum/scalar(@$cnes);  # score mean (approximate solution)
    my $cigar = join(';',@cigar_parts); # concateate the cigar info into a "meta-cigar string" - not so pretty
    return cne->new(start1 => $start1,
		    start2 => $start2,
		    end1 => $end1,
		    end2 => $end2,
		    strand => $strand,
		    score => $score,
		    cigar => $cigar);
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

