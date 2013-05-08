#!/usr/bin/perl -w
use strict;

use Getopt::Std;
use Class::Struct;
use AT::GFX::Locus;
use AT::GFX::AlignedLoci;
use AT::GFX::LocusConfig;


my $dbuser = "engstrom";
my $dbpass = "cog117";

my %args;
getopts('c:f:s:w:', \%args);
my $MIN_EXON_SIZE = 0;
my $OUTPUT_FORMAT = $args{'f'} || 'png';
my $MAX_NR_SEC_PANELS = $args{'s'};
my $WIDTH = $args{'w'};
my $CONN_HEIGHT = $args{'c'};
my $OUT_FN = pop @ARGV;
my ($CONFIG_FN, @LOC) = @ARGV;
unless(@LOC) {
    print STDERR <<EndOfUsage;

The draw_loci.pl script creates images of genomic features at single loci or multiple aligned loci.

Usage:
perl $0 [options] config-file assembly1:region1 [assembly2[:region2] [assembly3:region3]] out-file

Regions should be specified as chr:start-end.

Options:
 -f png|svg       Output format (default png).
                  Note: svg output is not implemented for aligned loci
 -w <width>       Width of (largest) reference region in pixels
 -c <height>      Height of connector polygons in pixels
 -s <nr>          Max number of secondary panels to show

EndOfUsage
exit;
}

# Check output format
my $OUTPUT_SVG;
if($OUTPUT_FORMAT eq 'svg') {
    $OUTPUT_SVG = 1;
}
elsif($OUTPUT_FORMAT ne 'png') {
    die "Invalid output format $OUTPUT_FORMAT. Valid formats: svg, png.\n";
}

# Read config file
my $CONFIG = AT::GFX::LocusConfig->new(file => $CONFIG_FN);

# Parse loci arguments
my (@assemblies, @loci);
my $i = 0;
foreach my $loc_str (@LOC) {
    my ($asm_name, $chr, $pos) = split ':', $loc_str;
    my $asm = $CONFIG->assemblies->{$asm_name};
    unless($asm) {
	die "Assembly $asm_name given on command line is not defined in $CONFIG_FN\n";
    }
    $assemblies[$i] = $asm;
    if($chr and $pos) {
	my ($start, $end) = split '-', $pos;
	$start =~ s/,//g;
	$end =~ s/,//g;
	push @loci, [$i, $chr, $start, $end];
    }
    $i++;
}
die "No loci specified.\n" unless(@loci);

# connect to databases and set some parameters
my @init_params = map { $_->get_drawing_parameters } @assemblies;

# prepare to draw image: connect to db and construct drawing object
my ($drawer, @draw_params);
if(@assemblies == 1) {
    $drawer = AT::GFX::Locus->new(%{$init_params[0]},
				  width => $WIDTH,
				  min_exon_size => $MIN_EXON_SIZE,
				  svg => $OUTPUT_SVG);
    @draw_params = @{$loci[0]}[1..3];
}
elsif(@assemblies == 2 or @assemblies == 3) {
    $drawer = AT::GFX::AlignedLoci->new(assembly_param => \@init_params,
					width => $WIDTH,
					connector_height => $CONN_HEIGHT,
					min_exon_size => $MIN_EXON_SIZE,
					max_nr_secondary_panels => $MAX_NR_SEC_PANELS,
					svg => $OUTPUT_SVG);
    push @draw_params, \@loci;
}
else {
    die "Please specify 1 - 3 assemblies on the commandline.\n";
}


if($drawer->draw_locus(@draw_params)) {
    open OUT, ">$OUT_FN";
    my $gd = $drawer->gd_image;
    print OUT $OUTPUT_SVG ? $gd->svg : $gd->png;
    close OUT;
}
else {
    print STDERR $drawer->errstr, "\n";
}

