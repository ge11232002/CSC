#!/usr/bin/perl
#
# blat_filter.pl
#
# blat HCNEs against the respective genomes
# remove features according to given maximal count and similarity
# thresholds.
#
# Copyright Boris Lenhard research group, University of Bergen

use warnings;
use strict;
use FindBin '$Bin';
use File::Temp qw/ tempfile tempdir /;
use File::Basename;
use Getopt::Long;
use AT::BlatpslxIO;

my $TWOBIT_DIR = '/export/data/goldenpath';
my $DEF_TEMP_DIR = 'blat_filter.tmp';
my $DEF_CUT_IDENTITY = 90;
my $DEF_COUNT_FN = "$Bin/blat_hit_cutoffs.txt";
my $DEF_BLAT_OPT_WSLO =  '-tileSize=9 -minScore=24 -repMatch=16384';
my $DEF_BLAT_OPT_WSMID = '-tileSize=10 -minScore=28 -repMatch=4096';
my $DEF_BLAT_OPT_WSHI =  '-tileSize=11 -minScore=30 -repMatch=1024';

# parse args
my ($TWOBIT_ASM1, $TWOBIT_ASM2, $PSL_ASM1, $PSL_ASM2, $CUT_COUNT_ASM1, $CUT_COUNT_ASM2);
my $TEMP_DIR = $DEF_TEMP_DIR;
my $CUT_IDENTITY = $DEF_CUT_IDENTITY;
my $COUNT_FN = $DEF_COUNT_FN;
my $BLAT_OPTIONS = 'default';
GetOptions('asm1=s' => \$TWOBIT_ASM1,
	   'asm2=s' => \$TWOBIT_ASM2,
	   'psl1=s' => \$PSL_ASM1,
	   'psl2=s' => \$PSL_ASM2,
	   'tmp=s' => \$TEMP_DIR,
	   'id=f' => \$CUT_IDENTITY,
	   'count1=i' => \$CUT_COUNT_ASM1,
	   'count2=i' => \$CUT_COUNT_ASM2,
	   'countfile=s' => \$COUNT_FN,
           'blatopts=s' => \$BLAT_OPTIONS);

my @CNE_FILES = @ARGV;

usage() unless (@CNE_FILES);

main();
exit 0;


sub usage
{
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;

$cmd

Usage: perl $cmd [options] <CNE files>

Options:
 --asm1 <2bit-file1>  .2bit file for assembly 1 (first species in CNE file)
 --asm2 <2bit-file2>  .2bit file for assembly 2 (second species in CNE file)
 --psl1 <psl-file1>    Existing psl file from Blat run against assembly 1
 --psl2 <psl-file2>    Existing psl file from Blat run against assembly 2
 --tmp <dir>           Directory for temporary files (default $DEF_TEMP_DIR)
 --id <percentage>     %-identity cutoff above which hits are counted
                       (default $DEF_CUT_IDENTITY)
 --count1 <count>      Max. nr of hits above cutoff allowed in asm1
 --count2 <count>      Max. nr of hits above cutoff allowed in asm2
 --countfile <file>    File mapping assembly ids to max. nr of hits
                       (default $DEF_COUNT_FN)
 --blatopts <string>   Options for blat. By default set according to window size
                       of CNE set:
                            win <= 35: '$DEF_BLAT_OPT_WSLO'
                       35 < win <= 45: '$DEF_BLAT_OPT_WSMID'
		       45 < win:       '$DEF_BLAT_OPT_WSHI'

EndOfUsage
    exit;
}


sub main
{

    # Create temp directory unless it exists
    unless(-e $TEMP_DIR) {
	mkdir($TEMP_DIR) or die "Could not create directory $TEMP_DIR";
    }

    # Read count cutoff file unless count cutoff specified
    my $count_cutoff_index;
    unless($CUT_COUNT_ASM1 and $CUT_COUNT_ASM2) {
	$count_cutoff_index = read_count_cutoffs($COUNT_FN);
    }

    # Process each file
    foreach my $cne_fn (@CNE_FILES) {

	# Try to get some info from filename
	my ($cne_basename) = fileparse($cne_fn);
	my ($asm1, $asm2, $cutoff, $win_size) = $cne_basename =~ /cne2w_(\w+)_(\w+)_(\d+)_(\d+)$/;
	my $parsable_filename = ($asm1 and $asm2 and $cutoff and $win_size) ? 1 : 0;
	unless($parsable_filename or
	       (($TWOBIT_ASM1 or $PSL_ASM1 ) and ($TWOBIT_ASM2 or $PSL_ASM2) and $CUT_COUNT_ASM1 and $CUT_COUNT_ASM1)) {
	    die "Could not parse species and cutoff from filename $cne_basename. Need parsable filename or corresponding information as options.\n";
	}

	# Get count and identity cutoffs
	my $cut_count_asm1 = $CUT_COUNT_ASM1 || $count_cutoff_index->{$asm1} or die "No count cutoff defined for $asm1\n";
	my $cut_count_asm2 = $CUT_COUNT_ASM2 || $count_cutoff_index->{$asm2} or die "No count cutoff defined for $asm2\n";

	# Create outfile filename
	my $outfile_fn = $parsable_filename ? "cne2wBf_${asm1}_${asm2}_${cutoff}_${win_size}" : "$cne_fn.Bf";

	# Read CNEs
	my $cnes = read_cnes($cne_fn);

	# Get blat options
	my $blat_opt;
	if($BLAT_OPTIONS ne 'default') {
	    $blat_opt = $BLAT_OPTIONS;
	}
	elsif(!defined($win_size) or $win_size > 45) {
	    $blat_opt = $DEF_BLAT_OPT_WSHI;
	}
	elsif($win_size > 35) {
	    $blat_opt = $DEF_BLAT_OPT_WSMID;
	}
	else {
	    $blat_opt = $DEF_BLAT_OPT_WSLO;
	}

	# Run blat unless we have PSL files
	my $psl_asm1 = $PSL_ASM1 || run_blat($cnes, 1, $TWOBIT_ASM1 || "$TWOBIT_DIR/$asm1/assembly.2bit", $CUT_IDENTITY, $blat_opt);
	my $psl_asm2 = $PSL_ASM2 || run_blat($cnes, 2, $TWOBIT_ASM2 || "$TWOBIT_DIR/$asm2/assembly.2bit", $CUT_IDENTITY, $blat_opt);

	# Parse blat results
	process_psl($psl_asm1, $psl_asm2, $cnes, $CUT_IDENTITY, $cut_count_asm1, $cut_count_asm2, $outfile_fn);
	
	# Remove psl files
	unlink($psl_asm1);
	unlink($psl_asm2);
    }

    # Remove temporary directory if empty
    my @tmp_files = glob("$TEMP_DIR/*");
    if(@tmp_files == 0) {
	rmdir($TEMP_DIR) or warn "WARNING: could not remove temporary directory $TEMP_DIR: $!\n";
    }
}


sub read_cnes
{
    my $cne_fn = shift;
    my @cnes;
    open(CNE, "$cne_fn") or die "Could not open $cne_fn\n";
    while(<CNE>){
	my @cne = split;
	push @cnes, \@cne;
    }
    close CNE;
    return \@cnes;
}


sub read_count_cutoffs
{
    my $fn = shift;
    open IN, $fn or die "Could not open $fn\n";
    my %index;
    while(my $line = <IN>) {
	chomp $line;
	next if($line =~ /^\s*$/ or $line =~ /^\s*#/);	
	my ($id, $cutoff) = split /\s+/, $line;
	die "Error parsing file $fn, line: $line\n" unless($id and $cutoff);
	$index{$id} = $cutoff;
    }
    close IN;
    return \%index;
}


sub run_blat
{
    my ($cnes, $asm_nr, $twobit_fn, $cut_identity, $blat_opt) = @_;
    $cut_identity = int($cut_identity);

    die "invalid asm_nr $asm_nr" unless($asm_nr == 1 or $asm_nr == 2);

    # Create file with CNE coordinates for blat
    my $temp_cne = File::Temp::tempnam($TEMP_DIR, "asm$asm_nr");
    open(OUT, ">$temp_cne") or die;
    my %seen;
    foreach my $cne (@$cnes) {
	my $id;
	if($asm_nr == 1) {
    # The blat takes the start as 0-based, end as 1-based. Here the old script is wrong!
    #$id = $cne->[0].":".($cne->[1]+1)."-".$cne->[2];
    $id = $cne->[0].":".$cne->[1]."-".$cne->[2];
	}
	else {
    #$id = $cne->[3].":".($cne->[4]+1)."-".$cne->[5];
    $id = $cne->[3].":".$cne->[4]."-".$cne->[5];
	}
	next if($seen{$id});
	$seen{$id} = 1;
	print OUT "$twobit_fn:$id\n";
    }
    close OUT;

    # Run blat
    my $temp_blat_out = File::Temp::tempnam($TEMP_DIR, "blat$asm_nr");
    my $blat_cmd =  "blat $blat_opt -minIdentity=$cut_identity $twobit_fn $temp_cne $temp_blat_out\n";
    print STDERR "$blat_cmd\n";
    system($blat_cmd) == 0 or die "Command '$blat_cmd' failed";

    # Remove CNE file 
    unlink($temp_cne);

    # Return name of blat output file
    return $temp_blat_out;
}


sub process_psl
{
    my ($psl_asm1, $psl_asm2, $cnes, $cut_identity, $cut_count_asm1, $cut_count_asm2, $out_fn) = @_;
    $cut_identity = $cut_identity/100;

    # Open psl files
    open (PSL_ASM1, $psl_asm1) or die "Could not open $psl_asm1\n";
    open (PSL_ASM2, $psl_asm2) or die "Could not open $psl_asm2\n";
    my $map_asm1 = AT::BlatpslxIO->new(-fh => \*PSL_ASM1);
    my $map_asm2 = AT::BlatpslxIO->new(-fh => \*PSL_ASM2);

    my (%count_tot_asm1,
	%count_above_id_asm1,
	%count_tot_asm2,
	%count_above_id_asm2); 

    # Read mappings on asm1
    print STDERR "Reading mappings on assembly1..\n";
    while( my $m = $map_asm1->next_mapping() ){
	$count_tot_asm1{$m->qName}++;
	if($m->matches/$m->qSize >= $cut_identity){
	    $count_above_id_asm1{$m->qName}++;
	}
    }

    # Read mappings on asm2
    print STDERR "Reading mappings on assembly2..\n";
    while(my $m = $map_asm2->next_mapping() ){
	$count_tot_asm2{$m->qName}++;
	if($m->matches/$m->qSize >= $cut_identity){
	    $count_above_id_asm2{$m->qName}++;
	}
    }

    # Filter CNEs
    open OUT, ">$out_fn"; # output file
    foreach my $cne (@$cnes) {
	my ($chr1,$start1,$end1,$chr2,$start2,$end2,$score) = @$cne;
	
  my $asm1_loc_id = $chr1.":".($start1+1)."-".$end1;
	my $asm2_loc_id = $chr2.":".($start2+1)."-".$end2;

	my $nr_hits_above_id_asm1 = defined($count_above_id_asm1{$asm1_loc_id}) ? $count_above_id_asm1{$asm1_loc_id} : 0;
	my $nr_hits_above_id_asm2 = defined($count_above_id_asm2{$asm2_loc_id}) ? $count_above_id_asm2{$asm2_loc_id} : 0;
	my $nr_hits_tot_asm1 = defined($count_tot_asm1{$asm1_loc_id}) ? $count_tot_asm1{$asm1_loc_id} : 0;
	my $nr_hits_tot_asm2 = defined($count_tot_asm2{$asm2_loc_id}) ? $count_tot_asm2{$asm2_loc_id} : 0;
	
	# Reject HCNEs that have more than a given nr of hits above a given threshold in either assembly
	if($nr_hits_above_id_asm1 <= $cut_count_asm1 and $nr_hits_above_id_asm2 <= $cut_count_asm2) {
	    print OUT join "\t", (@$cne,
				  $nr_hits_above_id_asm1, $nr_hits_above_id_asm2, $count_tot_asm1{$asm1_loc_id}||0,$count_tot_asm2{$asm2_loc_id}||0), "\n";
	}
    }

    # Close filehandles
    close PSL_ASM1;
    close PSL_ASM2;
    close OUT;
}
