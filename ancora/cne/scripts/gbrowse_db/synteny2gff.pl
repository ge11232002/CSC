#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

my %args;
getopts('a:g:',\%args);

my $MIN_ALI = $args{a} || 2000;
my $MIN_GAP = $args{g} || 20;

sub usage {
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;

$cmd - convert synteny files from bed to gff 
    
To convert joinedNet files:
Usage: $cmd jn <target-id> <summary-file> <bed-file>

To convert chainNetSynteny file:
Usage: $cmd cns <target-id> <bed-file>

Options:
-a <n>   Only include blocks with at least n aliged bases (for joinedNet only)
-g <n>   Close gaps of smaller than n
    
EndOfUsage
    exit;
}

my $mode = shift @ARGV;
usage() unless($mode);

if($mode eq 'jn') {
    my ($target_asm_id, $summary_file, $bed_file) = @ARGV;
    usage() unless($bed_file);
    convert_joined_net($target_asm_id, $summary_file, $bed_file);
}
elsif($mode eq 'cns') {
    my ($target_asm_id, $bed_file) = @ARGV;
    usage() unless($bed_file);
    convert_chain_net_synteny($target_asm_id, $bed_file);
}


sub convert_joined_net
{
    my ($target_asm_id, $summary_file, $bed_file) = @_;

    open SUMMARY, $summary_file or die "could not open $summary_file";
    open BED, $bed_file or die "could not open $bed_file";
    
    my $i = 1;
    
    while(my $bed_line = <BED>) {
	
	# Skip non-data lines in bed file
	next if($bed_line =~ /^track/ or $bed_line =~ /^browser/);
	
	# Read line from summary file
	my $summary_line = <SUMMARY>;
	
	# Parse input
	chomp $bed_line;
	chomp $summary_line;
	my ($chr1, $start1, $end1, $chr2, $start2, $end2, $nr_blocks, $status, $ali) = split "\t", $summary_line; 
	next unless($status eq 'retained');
	my ($chr, $start, $end, undef, undef, $strand, undef, undef, undef, undef, $blockSizes_str, $blockStarts_str) = split "\t", $bed_line;
	$start++;
	
	# Check if target is first or second in the summary file
	my ($target_chr, $target_start, $target_end);
	if($chr eq $chr1 and $start == $start1 and $end == $end1) {
	    ($target_chr, $target_start, $target_end) = ($chr2, $start2, $end2);
	}
	elsif($chr eq $chr2 and $start == $start2 and $end == $end2) {
	    ($target_chr, $target_start, $target_end) = ($chr1, $start1, $end1);
	}
	else {
	    die "input files are not ordered in the same way";
	}
	
	# Ignore if ali is below threshold
	next if($ali < $MIN_ALI);
	
	# Generate some strings
	my $source = "join_net_$target_asm_id";
	my $id = "$target_asm_id.$i";
	my $short_name = $target_chr . ':' . int($target_start/1000) . 'k';
	my $long_name = "$target_chr:$target_start..$target_end";
	
	# Output GFF
	write_gff($chr, $source, $start, $end, $ali, $strand, $id, $short_name, $long_name,
		  $blockStarts_str, $blockSizes_str);
	
	# Increment block id
	$i++;
    }
}


sub convert_chain_net_synteny
{
    my ($target_asm_id, $bed_file) = @_;

    open BED, $bed_file or die "could not open $bed_file";
    
    my $i = 1;
    
    while(my $bed_line = <BED>) {
	
	# Skip non-data lines in bed file
	next if($bed_line =~ /^track/ or $bed_line =~ /^browser/);
	
	# Parse input
	chomp $bed_line;
	my ($chr, $start, $end, $long_name, $score, $strand, undef, undef, undef, undef, $blockSizes_str, $blockStarts_str) = split "\t", $bed_line;
	$start++;
	my ($target_chr1, $target_start1) = $long_name =~ /^([^:]+):(\d+)/;
	
	# Generate some strings
	my $source = "chainNetSynteny_$target_asm_id";
	my $id = "$target_asm_id.$i";
	my $short_name = $target_chr1 . ':' . int($target_start1/1000) . 'k';
	$long_name =~ s/-/../g;

	# Output GFF
	write_gff($chr, $source, $start, $end, $score, $strand, $id, $short_name, $long_name,
		  $blockStarts_str, $blockSizes_str);

	# Increment block id
	$i++;
    }
}


sub write_gff
{
    my($chr, $source, $start, $end, $score, $strand, $id, $short_name, $long_name,
       $blockStarts_str, $blockSizes_str) = @_;

    # Write GFF line for entire block
    print join("\t",
	       $chr, $source, "synteny_block",
	       $start, $end,
	       $score, $strand, '.',
	       "SyntenyBlock $id ; ShortName $short_name ; Match $long_name"
	       ), "\n";
	
    # Write GFF line for block parts
    my ($block_starts, $block_sizes) = get_block_starts_and_sizes($blockStarts_str, $blockSizes_str, $MIN_GAP);
    while(@$block_starts) {
	my $rel_block_start = shift @$block_starts;
	my $block_size = shift @$block_sizes;
	my $block_start = $start + $rel_block_start;
	my $block_end = $block_start + $block_size - 1;
	print join("\t",
		   $chr, $source, "synteny_part",
		   $block_start, $block_end, '.',  $strand, '.',
		   "SyntenyBlock $id" 
		   ), "\n";
    }
    
}


sub get_block_starts_and_sizes
{
    my ($blockStarts_str, $blockSizes_str, $min_gap) = @_;

    my @blockStarts = split ',', $blockStarts_str;
    my @blockSizes = split ',', $blockSizes_str;
    my @newBlockStarts = (shift @blockStarts);
    my @newBlockSizes = (shift @blockSizes);
    my $lastEnd = $newBlockStarts[0] + $newBlockSizes[0];

    while(@blockStarts) {
	my $blockStart = shift @blockStarts;
	my $blockSize = shift @blockSizes;
	my $gapSize = $blockStart - $lastEnd;
	if($gapSize >= $min_gap)
	{
	    push @newBlockStarts, $blockStart;
	    push @newBlockSizes, $blockSize;
	}
	else
	{
	    $newBlockSizes[-1] += $blockSize + $gapSize;
	}
	$lastEnd = $blockStart + $blockSize;
    }

    return (\@newBlockStarts, \@newBlockSizes);
}
