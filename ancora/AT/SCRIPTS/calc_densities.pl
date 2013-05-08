#!/usr/bin/perl
#
# calc_densities.pl
# Calculate feature densities on chromosomes. Run script without arguments for usage.
#
# Copyright Boris Lenhard research group, BCCS, UNIFOB AS, University of Bergen, Norway

use warnings;
use strict;
use AT::DB::GenomeAssemblyTwoBit;
use AT::DB::GenomeAssemblyNibs;
use AT::Tools::RangeHandler;
use Getopt::Std;

my $BLK_IDX_MASK = 0x1ff;  # = 511  note: must exceed $WIN_SIZE

# parse args
my %args;
getopts('ac:d:em:n:o:p:s', \%args);
my ($ASM_PATH, $FT_FN, $BLK_SIZE, $WIN_SIZE) = @ARGV;

unless ($WIN_SIZE) {
print STDERR <<EndOfUsage;

Usage: perl $0 [options] <assembly> <feature-file> <step-size> <win-size>

Mandatory arguments:
 <assembly>      .2bit file or directory containing .nib files
 <feature-file>  Features to calculate density of; bed format
 <step-size>     Step size in basepairs
 <win-size>      Window size specified as number of steps; must be <$BLK_IDX_MASK

Options:
 -a              Output absolute score (number of bases covered by features)
                 instead of percentage
 -c <percentage> Set wig track height to clip the N% highest values 
                 E.g. "1" to clip 1% highest scores
 -d <desc>       Description of track (for wig/bed track)
 -e              Do not include empty windows in output (scores of 0)
 -m <mask-file>  Mask-file in bed format
 -n <name>       Name of track (for wig/bed track) or method (for gff)
 -o <format>     Output format (wig, bed or gff); default is wig
 -p <string>     Parameter string for bed/wig track line  
                 E.g. "visibility=full color=255,0,0"
		 Note that name and description should not be part of this
                 string as they are set by other other options.
		 viewLimits can be set in this string to override the values
                 automatically computed by the script.
 -s              Stack overlapping features. Default is to collapse them.

Density values are calculated as the percentage of bases within a window that is
covered by features in the feature file. In stacking mode (-s option), bases that
are covered by several features will be counted several time. E.g. if the entire
window is covered by exactly two features, the density will be 100 in default
mode and 200 in stacking mode. If a mask file is given, density values will be
the percentage of unmasked bases within a window that is covered by features in
the feature-file.

Features in the feature-file and mask-file must not overlap. If they do, the script
will report this and terminate.

The output is given on stdout.

Note that to use assemblies in .2bit format, the executables twoBitToFa and
twoBitInfo must be in your \$PATH.

EndOfUsage
    exit;
}

my $SKIP_EMPTY = $args{'e'};
my $ABS_SCORE = $args{'a'};
my $MASK_FN = $args{'m'};
my $OUTPUT_FORMAT = $args{'o'} || 'wig';
my $GENERIC_NAME = "dens_".($BLK_SIZE/1000)."k_".$WIN_SIZE; 
my $WIG_NAME = $args{'n'} || $GENERIC_NAME;
my $WIG_DESC = $args{'d'} || $GENERIC_NAME;
my $GFF_SOURCE = "calc_densities";
my $GFF_METHOD = $args{'n'} || $GENERIC_NAME;
my $CLIP_PERCENTAGE = $args{'c'} || 0;
my $TRACK_PARAM_STR = $args{'p'} || '';
my $STACK_MODE = $args{'s'};

# Check arguments
if($WIN_SIZE < 1 or $WIN_SIZE >= $BLK_IDX_MASK) {
    die "Invalid win-size\n";
}
if($OUTPUT_FORMAT ne 'wig' and $OUTPUT_FORMAT ne 'gff' and $OUTPUT_FORMAT ne 'bed') {
   die "Invalid output format $OUTPUT_FORMAT\n";
}


# Create a genome assembly object
my $gendb;
if(-d $ASM_PATH) {
    $gendb = AT::DB::GenomeAssemblyNibs->new(dir => $ASM_PATH)
        or die "Could not find .nib files at $ASM_PATH";
}
else {
    $gendb = AT::DB::GenomeAssemblyTwoBit->new(file => $ASM_PATH)
	or die "Could not open .2bit file $ASM_PATH";
}

# Read features to calculate density for
print STDERR "Reading features...\n";
my $fts_by_chr = read_bed($FT_FN);

# Read mask
my $mask_by_chr;
if($MASK_FN) {
    print STDERR "Reading mask...\n";
    $mask_by_chr = read_bed($MASK_FN);
}

# Print track header if bed output format
print "track name=\"$WIG_NAME\" description=\"$WIG_DESC\" $TRACK_PARAM_STR\n" if($OUTPUT_FORMAT eq 'bed');

# Process each chromosome
my %win_scores_by_chr;
foreach my $chr (sort keys %$fts_by_chr) {

    # Get chromosome size
    my $chr_size = $gendb->get_chr_size($chr);
    if($chr_size < $WIN_SIZE*$BLK_SIZE) {
	warn "Skipping $chr ($chr_size bp) because it is shorter than the window size.\n";
	next;
    }
    print STDERR "Processing $chr...\n";

    # Get sorted features (sort by start pos)
    my $fts = [ sort {$a->[0] <=> $b->[0]} @{$fts_by_chr->{$chr} || []} ];
    my $mask = [ sort {$a->[0] <=> $b->[0]} @{$mask_by_chr->{$chr} || []} ];

    # Calculate collapsed features and mask
    my $collapsed_fts = collapse_intervals($fts);
    my $collapsed_mask = collapse_intervals($mask);
    unless($STACK_MODE) {
	$fts = $collapsed_fts;
	$mask = $collapsed_mask;
    }

    # Check that features and mask do not overlap
    my $shared_regions = AT::Tools::RangeHandler->compute_intersection($collapsed_fts, $collapsed_mask);
    if (@$shared_regions) {
	my $int = $shared_regions->[0];
	die "Features overlap with mask at", $int->[0],"-",$int->[1]," (",$int->[1]-$int->[0]+1," bp)";
    }
 
    # Calculate the window scores
    # In case wig format is requested, we save the scores in memory and print them when done
    # This is necessary to set the scale in the wig header
    my $print_and_discard = $OUTPUT_FORMAT eq 'wig' ? 0 : 1;
    my $win_scores = calc_window_scores($chr, '.', 1, $chr_size, $fts, $mask,
					$print_and_discard);
    $win_scores_by_chr{$chr} = $win_scores;
}

if($OUTPUT_FORMAT eq 'wig') {
    print STDERR "Printing result...\n";
    print_wig_track(\%win_scores_by_chr, $WIG_NAME, $WIG_DESC, $TRACK_PARAM_STR);
}

print STDERR "Done!\n";


sub calc_window_scores
{
    my ($chr, $strand, $start, $end, $ivs1, $ivs2, $print) = @_;
    my $offset = int((($WIN_SIZE-1)*$BLK_SIZE)/2+0.5);
    my $max_score = $BLK_SIZE * $WIN_SIZE;

    my @all_win_data;
    my $size = $end - $start + 1;
    my $nr_blocks = int($size / $BLK_SIZE) + ($size % $BLK_SIZE ? 1 : 0);
    my ($prev_win_score1, $prev_win_score2) = (0,0);
    my @blk_scores1 = (0) x ($BLK_IDX_MASK+1);
    my @blk_scores2 = (0) x ($BLK_IDX_MASK+1);
    my ($blk_start, $blk_end) = ($start, $start + $BLK_SIZE-1);
    my ($j,$k) = (0,0);   # $j,$k = interval indices
    for my $i (1..$nr_blocks) {  # $i = block index
	my $blk_score1 = calc_block_score($blk_start, $blk_end, \$j, $ivs1);
	my $blk_score2 = calc_block_score($blk_start, $blk_end, \$k, $ivs2);
	my $block = $i & $BLK_IDX_MASK;
	my $dead_block = ($i - $WIN_SIZE) & $BLK_IDX_MASK;
	my $win_score1 = $prev_win_score1 +$blk_score1 -$blk_scores1[$dead_block];
	my $win_score2 = $prev_win_score2 +$blk_score2 -$blk_scores2[$dead_block];
	if((!$SKIP_EMPTY or $win_score1) and $blk_start >= $start + $offset) {   
	    my $relative_score = 100 * $win_score1/($max_score-$win_score2);

	    # score transformation: to improve signal/noise and improve resolution, the density scores are log-transformed
	    # as follows: S' = 1-(log(2-S)/log(2))   #means 1-"log2" of (2-S)

	    # log-transform score (pushed down low values, improves resolution for values close to 1)
	    # my $transformed_score = 1-(log(2-$relative_score)/log(2));

	    my @win_data = ($chr, $strand, $blk_start-$offset, $blk_end-$offset,
			    $ABS_SCORE ? $win_score1 : $relative_score);
	    if($print) {
		print_win_score(@win_data);
	    }
	    else {
		push @all_win_data, \@win_data;
	    }
	}
	$blk_scores1[$block] = $blk_score1;
	$blk_scores2[$block] = $blk_score2;
	$prev_win_score1 = $win_score1;
	$prev_win_score2 = $win_score2;
	$blk_start += $BLK_SIZE;
	$blk_end += $BLK_SIZE;
    }

    return \@all_win_data;
}


sub calc_block_score
{
    my ($blk_start, $blk_end, $j_ref, $ivs) = @_;
    my $blk_score = 0;
    my $j = $$j_ref;
    my $next_start;
    for(;$j < @$ivs; $j++)
    {
	my ($iv_start, $iv_end) = @{$ivs->[$j]};
	last if($iv_start > $blk_end);
	if($iv_end >= $blk_start) {
	    my $ol_start = ($iv_start > $blk_start) ? $iv_start : $blk_start;
	    my $ol_end = ($iv_end < $blk_end) ? $iv_end : $blk_end;
	    $blk_score += $ol_end - $ol_start + 1;
	    $next_start = $j if(!defined($next_start) and $iv_end > $blk_end);
	}
    }
    $$j_ref = defined($next_start) ? $next_start : $j;
    return $blk_score;
}


sub collapse_intervals
{
    my ($in) = @_;
    # NOTE: input regions must be sorted by start pos

    return [] unless($in and @$in);

    my @out;
    my ($start, $end) = @{$in->[0]};
    for my $i (1..@$in-1) {
	my ($my_start, $my_end) = @{$in->[$i]};
	if($my_start > $end) {
	    push @out, [$start,$end];
	    ($start, $end) = ($my_start, $my_end); 
	}
	else {
	    $end = $my_end if ($end < $my_end);
	}
    }
    push @out, [$start,$end];
    return \@out;
}


sub print_win_score
{
    if($OUTPUT_FORMAT eq 'bed') {
	print_win_score_as_bed(@_);
    }
    elsif($OUTPUT_FORMAT eq 'gff') {
	print_win_score_as_gff(@_);
    }
    else {
	die "unknown output format $OUTPUT_FORMAT";
    }
}


sub print_win_score_as_bed
{
    my ($chr, $strand, $start, $end, $score) = @_;
    $score = sprintf("%.2f",$score) unless($ABS_SCORE);
    #return unless($score);
    print join("\t",
	       $chr,
	       $start-1,
	       $end,
	       '.',
	       $score),"\n";
}


sub print_win_score_as_gff
{
    my ($chr, $strand, $start, $end, $score) = @_;
    $score = sprintf("%.2f",$score) unless($ABS_SCORE);
    #return unless($score);
    print join("\t",
        ($chr,
         $GFF_SOURCE,
         $GFF_METHOD,
         $start,
         $end,
         $score,
         $strand,
         ".",
         "$GFF_METHOD $chr:$GFF_SOURCE")),"\n";
}


sub print_wig_track
{
    my ($windows_by_chr, $name, $desc, $param_str) = @_;

    my @params = ('type=wiggle_0', "name=\"$name\"", "description=\"$desc\"", $param_str);

    unless($param_str =~ /autoScale=/) {
	push @params, 'autoScale=off';
    }

    unless($param_str =~ /viewLimits=/) {
	# Find max window score to show
	my $max_score_shown = 0;
	if($CLIP_PERCENTAGE) {
	    # Clipping mode: clip top $CLIP_PERCENTAGE % of scores.
	    # ALSO should this be per chromosome?
	    my $clip_fraction = 1 - ($CLIP_PERCENTAGE/100);
	    my @all_scores;
	    foreach my $window_list (values %$windows_by_chr) {
		foreach my $win (@$window_list) {
		    my $score = $win->[4];
		    push @all_scores, $score if($score);
		    # Ignore scores of 0 to avoid a situation where all or most positive scores are clipped
		}
	    }
	    @all_scores = sort @all_scores;
	    my $max_score_shown_index = int($clip_fraction * scalar(@all_scores));
	    $max_score_shown = $all_scores[$max_score_shown_index];
	}
	else {
	    # Clipping not activated: just find the highest score
	    foreach my $window_list (values %$windows_by_chr) {
		foreach my $win (@$window_list) {
		    my $score = $win->[4];
		    $max_score_shown = $score if($score > $max_score_shown);
		}
	    }
	}
	# Round max score up to nearest multiple of .1
	$max_score_shown = int(10*$max_score_shown+1)/10;
	push @params, "viewLimits=0:$max_score_shown";
    }

    print join(' ', 'track', @params), "\n";

    if($SKIP_EMPTY) {
	while(my ($chr, $window_list) = each %$windows_by_chr) {
	    print "variableStep chrom=$chr span=$BLK_SIZE\n";
	    foreach my $win (@$window_list) {
		my (undef, undef, $start, undef, $score) = @$win;
		$score = sprintf("%.2f",$score) unless($ABS_SCORE);
		print $start, "\t", $score, "\n";
	    }
	}
    }
    else {
	while(my ($chr, $window_list) = each %$windows_by_chr) {
	    my $pos = $window_list->[0][2];
	    print "fixedStep chrom=$chr start=$pos step=$BLK_SIZE span=$BLK_SIZE\n";
	    foreach my $win (@$window_list) {
		my (undef, undef, $start, undef, $score) = @$win;
		unless($pos eq $start) { die "window step error: chr=$chr pos=$pos start=$start" } 
		$score = sprintf("%.2f",$score) unless($ABS_SCORE);
		print $score, "\n";
		$pos += $BLK_SIZE;
	    }
	}
    }
}


sub read_bed
{
    my ($file) = @_;
    my %loc_by_chr;
    open IN, $file or die "could not open $file";
    while(my $line = <IN>) {
	chomp $line;
	next if(!$line or $line =~/^browser/ or $line =~ /^track/);
	my ($chr, $start, $end) = split /\t/, $line;
	push @{$loc_by_chr{$chr}}, [$start+1,$end];
    }
    close IN;
    return \%loc_by_chr;
}


sub _debug_print_intervals
{
    my ($ivs) = @_;
    foreach my $iv (@$ivs) {
	print STDERR join(' ',@$iv),"\n";
    }
}
