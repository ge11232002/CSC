#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin '$Bin';
use Getopt::Std;
use Class::Struct;
use AT::DB::GenomeAssemblyTwoBit;

struct comparison => [ asm1 => '$', asm2 => '$', thresholds => '@' ];

my $ASM_DIR = "/export/data/goldenpath";
my $CNE_DIR = "/export/data/CNEs";
my $AXT_DIR = "/export/downloads/ucsc/axtNet";
my $TMP_DIR = "detect_cnes.tmp";
my $SCAN_BINARY = "ceScan";
my $MERGE_SCRIPT = "$Bin/merge_twoway_cne_sets.pl";
my $TWOBIT_INFO_BINARY = "twoBitInfo";

my %args;
getopts('bc:f',\%args);
my $MAKE_BED = $args{'b'};
my $KEEP_TEMP_FILES = $args{'f'};
my $COMPARISONS_FILE = $args{'c'} || "$Bin/comparisons.txt";

main();
exit;

sub asm_fn { my $asm = shift; return "$ASM_DIR/$asm/assembly.2bit"; }
sub sizes_fn { my $asm = shift; return "$ASM_DIR/$asm/assembly.sizes"; }
sub filter_fn { my $asm = shift; return "$CNE_DIR/$asm/filters/filter_regions.$asm.bed"; }
sub out_prefix { my ($asm1, $asm2) = @_; return "$TMP_DIR/cne_${asm1}_${asm2}"; }

sub usage
{
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd 

Main script for detecting CNEs.

Usage:
perl $cmd [options] all
perl $cmd [options] <asm1> <asm2> [thresholds]

Options:
-c <file>   Set comparisons file. Default: $COMPARISONS_FILE.
-t          Do not remove temporary files.

EndOfUsage
exit;
}


sub main
{
    # Some think it's silly to have a main function in Perl scripts
    # but it is useful for avoiding too many global variables.

    # Get list of comparisons generate CNEs for
    my $comparisons = get_comparisons(@ARGV);

    # create temporary directory
    unless(-e $TMP_DIR) {
	mkdir($TMP_DIR) or die "ERROR: failed to create directory $TMP_DIR: $!\n";
    }

    # Process each comparison
    foreach my $c (@$comparisons) {
	print STDERR 'Processing '.$c->asm1.' vs '.$c->asm2."...\n";
	# Make sure we have all files
	input_files_exist($c->asm1, $c->asm2);
	input_files_exist($c->asm2, $c->asm1);
	# Run ceScan to find CNEs
	run_scan($c->asm1, $c->asm2, $c->thresholds);
	run_scan($c->asm2, $c->asm1, $c->thresholds);
	# Run merge script to merge CNEs
	run_merge($c->asm1, $c->asm2, $c->thresholds);
	# Remove temporary files 
	unless($KEEP_TEMP_FILES) {
	    remove_tmp_files($c->asm1, $c->asm2, $c->thresholds);
	    remove_tmp_files($c->asm2, $c->asm1, $c->thresholds);
	}
    }

    # Remove temporary directory if empty
    my @tmp_files = glob("$TMP_DIR/*");
    if(@tmp_files == 0) {
	rmdir($TMP_DIR) or warn "WARNING: could not remove temporary directory $TMP_DIR: $!\n";
    }
}


sub file_must_exist
{
    my $fn = shift;
    die "ERROR: file $fn not found" unless(-e $fn);
}

sub input_files_exist
{
    my ($asm1, $asm2) = @_;

    my $asm_fn = asm_fn($asm1);
    my $sizes_fn = sizes_fn($asm1);
    my $filter_fn = filter_fn($asm1);

    # Check that 2bit file exists
    file_must_exist($asm_fn);

    # Check that size file exists, otherwise create it
    unless(-e $sizes_fn) {
	system($TWOBIT_INFO_BINARY, $asm_fn, $sizes_fn) == 0
	    or die "ERROR: failed to create file $sizes_fn\n";
	warn "Created file $sizes_fn\n";
    }

    # Check that filter file exists
    file_must_exist($filter_fn);

    # Check that axt files exist 
    if(-e "$AXT_DIR/$asm1/$asm1.$asm2.net.axt" or -e "$AXT_DIR/$asm1/$asm1.$asm2.net.axt.gz") {
	# single file - ok
    }
    else {
	# multiple files - check that we have one for each chrom that is larger that 10 Mb (an arbitrary threshold)
	my $asm = AT::DB::GenomeAssemblyTwoBit->new(file => $asm_fn);
	my $sizes = $asm->get_all_chr_sizes();
	my $nr_files_found = 0;
	while(my ($chr, $size) = each %$sizes) {    
	    my $axt_fn = "$AXT_DIR/$asm1/$chr.$asm1.$asm2.net.axt";
	    if(-e $axt_fn or -e "$axt_fn.gz") {
		$nr_files_found++;
	    }
	    elsif($size > 1e7) {
		warn "WARNING: alignment file $axt_fn.gz not found\n"
	    }
	}
	die "ERROR: no axtNet files from $asm1 to $asm2\n" unless $nr_files_found;
    }
}


sub run_scan
{
    my ($asm1, $asm2, $thresholds) = @_;

    # Get arguments for ceScan
    my $sizes_fn = sizes_fn($asm2);
    my $filter_fn1 = filter_fn($asm1);
    my $filter_fn2 = filter_fn($asm2);
    my $out_prefix = out_prefix($asm1, $asm2);
    my $threshold_str = join(' ', map { "$_->[0],$_->[1]" } @$thresholds);

    # Make axt file string
    my $axt_files;
    $axt_files = "$AXT_DIR/$asm1/$asm1.$asm2.net.axt";
    unless(-e $axt_files) {
	$axt_files = "$AXT_DIR/$asm1/$asm1.$asm2.net.axt.gz";
	unless(-e $axt_files) {
	    $axt_files = join(' ' , glob("$AXT_DIR/$asm1/*.$asm1.$asm2.net.axt $AXT_DIR/$asm1/*.$asm1.$asm2.net.axt.gz"));
	}
    }

    # Assemble commandline
    my $cmd = "$SCAN_BINARY -tFilter=$filter_fn1 -qFilter=$filter_fn2 -outPrefix=$out_prefix $sizes_fn $threshold_str $axt_files";

    # Scan or die
    system($cmd) == 0 or die "ERROR: failed to execute command: $cmd\n";
}


sub remove_tmp_files
{
    my ($asm1, $asm2, $thresholds) = @_;
    my $out_prefix = out_prefix($asm1, $asm2);
    foreach my $t (@$thresholds) {
	my $fn = "${out_prefix}_$t->[0]_$t->[1]";
	unlink($fn) or warn "WARNING: failed to remove temporary file $fn\n";    
    }
}


sub run_merge
{
    my ($asm1, $asm2, $thresholds) = @_;
    my $out_prefix1 = out_prefix($asm1, $asm2);
    my $out_prefix2 = out_prefix($asm2, $asm1);
    my $bed_flag = $MAKE_BED ? '-b' : '';
    foreach my $t (@$thresholds) {
	my $fn1 = "${out_prefix1}_$t->[0]_$t->[1]";
	my $fn2 = "${out_prefix2}_$t->[0]_$t->[1]";
	my $cmd = "perl $MERGE_SCRIPT $bed_flag -s $fn1 $fn2";
	system($cmd) == 0 or die "ERROR: failed to execute command: $cmd\n";
    }
}

sub get_comparisons
{
    my (@args) = @_;

    my $comparisons;

    if(@args == 1 and $ARGV[0] eq 'all') {
	$comparisons = read_comparisons($COMPARISONS_FILE);
    }
    elsif(@args == 2) {
	my ($asm1, $asm2) = @ARGV;
	($asm1, $asm2) = ($asm2, $asm1) if($asm1 gt $asm2);
	my $all_comparisons = read_comparisons($COMPARISONS_FILE);
	foreach my $cmp (@$all_comparisons) {
	    if($cmp->asm1 eq $asm1 and $cmp->asm2 eq $asm2) {
		$comparisons = [$cmp];
		last;
	    }
	}
	die "Comparison between $asm1 and $asm2 not found in $COMPARISONS_FILE\n" unless($comparisons);
    }
    elsif(@args >= 3) {
	my $asm1 = shift @args;
	my $asm2 = shift @args;
	($asm1, $asm2) = ($asm2, $asm1) if($asm1 gt $asm2);
	my $thresholds = parse_threshold_strings(@args);
	my $cmp = comparison->new(asm1 => $asm1,
				  asm2 => $asm2,
				  thresholds => $thresholds);
	$comparisons = [$cmp];
    }
    else {
	usage();
    }

    return $comparisons;
}


sub read_comparisons
{
    my $fn = shift;
    open IN, $fn or die "Could not open $fn\n";
    my $ln = 0;
    my @comparisons;
    my %asm_pair_index;
    while(my $line = <IN>) {
	$ln++;
	chomp $line;
	$line =~ s/#.*//;
	next if($line =~ /^\s*$/);
	my ($asm1, $asm2, @threshold_strings) = split /\s+/, $line;
	die "Too few columns on line $ln, file $fn\n" unless(@threshold_strings);
	($asm1, $asm2) = ($asm2, $asm1) if($asm1 gt $asm2);
	die "ERROR: duplicate comparison between $asm1 and $asm2 in $fn.\n" if($asm_pair_index{$asm1}{$asm2});
	$asm_pair_index{$asm1}{$asm2} = 1;
	my $thresholds = parse_threshold_strings(\@threshold_strings);
	push @comparisons, comparison->new(asm1 => $asm1,
					   asm2 => $asm2,
					   thresholds => $thresholds);
    }
    return \@comparisons;
}


sub parse_threshold_strings
{
    my $strings = shift;
    my @thresholds;
    foreach my $s (@$strings) {
	my ($min_score, $win_size) = $s =~ /^(\d+),(\d+)$/;
	die "ERROR: illegal threshold string $s\n"
	    unless($min_score and $win_size and $min_score > 0 and $win_size > 1 and $min_score <= $win_size);
	push @thresholds, [$min_score, $win_size];
    }    
    return \@thresholds;
}
