#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Std;
use MyPerlVars;
use CNE::DB;
use AT::GFX::SeqColorIndex;

my $DEF_WIN_SIZE = 300;
my $STEP_SIZE = 1000;
my $DEF_CLIP_PC = 1;
my $COLOR_INDEX = AT::GFX::SeqColorIndex->new();
my $ASM_DIR = '/export/data/goldenpath';
my $DENSITY_SCRIPT = 'perl /opt/www/AT/SCRIPTS/calc_densities.pl';

my %args;
getopts('c:d:h:w:', \%args);
my $DB_HOST = $args{'h'} || 'localhost';
my $DB_NAME = $args{'d'} || 'cne';
my $WIN_SIZE = $args{'w'} || $DEF_WIN_SIZE;
my $CLIP_PC = defined($args{'c'}) ? $args{'c'} : $DEF_CLIP_PC;

main();
exit;

sub usage
{
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

This script creates files with CNE locations (bed format) and CNE densities
(wig format). 

Usage: 
    
$cmd [options] <assembly ids>

Download files will be created for the listed assembly ids. For example:
- Creating download files for assemblies hg18, mm9 and danRer5 
  $cmd hg18 mm9 danRer5
- Creating Dowload files for all assemblies in the database (note: may take long!)
  $cmd all

The CNE data are retrieved from a MySQL database. Creating files for an assembly
generates one bed and one wig file for each CNE set available for that assembly
in the database. There is currently no way of telling the script to only create
files for a particular pairwise comparison.

The files are placed in the subdirectories created under the current working
directory. The directory structure created is identical to that in the Ancora
download section, so that the directories created by this script simply can be
moved to the download section. 

Options:
-d <database name>     Optional. Default: cne.
-h <database host>     Optional. Default: localhost.
-w <window size in kb> For densities. Default: $DEF_WIN_SIZE kb.
-c <clip percentage>   For densities. Default: $DEF_CLIP_PC%.

EndOfUsage

    exit;
}


sub main
{
    my @asm_list = @ARGV;

    usage() unless(@asm_list);

    # Connect to db
    my $db = CNE::DB->connect(dbhost => $DB_HOST,
			      dbname => $DB_NAME,
			      dbuser => $MyPerlVars::sqlUser,
			      dbpass => $MyPerlVars::sqlPass)
	or die "could not connect to db $DB_NAME @ $DB_HOST";
    
    my $tables = $db->get_all_cne_table_names();
    
    # Find out which reference assemblies we should create files for
    my $req_asm_index;
    unless($asm_list[0] eq 'all') {
	$req_asm_index = { map { $_ => 1 } @asm_list };
    }
    
    # Create the files
    foreach my $t (@$tables) {
	make_files_for_table($db, $t, $req_asm_index);
    }
}   

   
sub make_files_for_table
{
    my ($db, $table, $req_asm_index) = @_;

    # get table info
    my $table_info = $db->get_cne_table_info($table);
    my $asm1 = $table_info->{assembly1};
    my $asm2 = $table_info->{assembly2};
    my $id = $table_info->{min_identity};
    my $len = $table_info->{min_length};

    return if($req_asm_index and not ($req_asm_index->{$asm1} or $req_asm_index->{$asm2}));

    # get assembly info
    my $asm1_info = $db->get_assembly_info($asm1);
    my $asm2_info = $db->get_assembly_info($asm2);
    unless($asm1_info) {
	warn "No information about assembly $asm1 in database; skipping table $table";
	return;
    }
    unless($asm2_info) {
	warn "No information about assembly $asm2 in database; skipping table $table";
	return;
    }

    # get all cnes
    my $cnes = $db->get_all_cnes_in_table(table_name => $table);

    # write bed & wig files
    if(!$req_asm_index or $req_asm_index->{$asm1}) {
	make_bed_and_wig_files($cnes, $asm1_info, $asm2_info, $id, $len, 0)
    }
    if(!$req_asm_index or $req_asm_index->{$asm2}) {
	make_bed_and_wig_files($cnes, $asm2_info, $asm1_info, $id, $len, 1)
    }
}


sub make_bed_and_wig_files
{
    my ($cnes, $asm1_info, $asm2_info, $id, $len, $swap) = @_;  
    my $bed_fn = write_bed_file($cnes, $asm1_info, $asm2_info, $id, $len, $swap);
    my $wig_fn = write_wig_file($asm1_info, $asm2_info, $id, $len, $bed_fn);
    print STDERR "Created files $bed_fn and $wig_fn.\n";
    system('gzip', '-f', $bed_fn) == 0 or die "gzip failed";
    system('gzip', '-f', $wig_fn) == 0 or die "gzip failed";
}

sub write_bed_file
{
    my ($cnes, $asm1_info, $asm2_info, $id, $len, $swap) = @_;  
  
    my $asm1 = $asm1_info->{assembly_id};
    my $asm2 = $asm2_info->{assembly_id};

    # Make track name & desc strings
    my $org2 = ucfirst make_organism_string($asm2_info);
    my $track_name = "$org2 HCNEs ".($id*100)."% / $len";
    my $track_desc = "$org2 HCNEs (".($id*100)."% identity over $len columns)";

    # Get file name and open file
    my $file_id = make_hcne_id_string($asm1, $asm2, $id, $len); 
    my $dir = get_dir_name($asm1_info, $asm2_info);
    my $file_name = "$dir/$file_id.bed";
    open OUT, ">$file_name" or die "could not open $file_name for output";

    # Print bed file header
    print OUT "track name=\"$track_name\" description=\"$track_desc\" itemRgb=On\n";

    # Print data
    if($swap) {
	foreach my $cne (@$cnes) {
	    print OUT join("\t",
			   $cne->chr2, $cne->start2-1, $cne->end2,
			   $cne->chr1.':'.$cne->start1.'-'.$cne->end1,
			   0, '.', $cne->start2-1, $cne->end2,
			   join(',',$COLOR_INDEX->get_seq_color($cne->chr1))), "\n";
	}
    }
    else {
	foreach my $cne (@$cnes) {
	    print OUT join("\t",
			   $cne->chr1, $cne->start1-1, $cne->end1,
			   $cne->chr2.':'.$cne->start2.'-'.$cne->end2,
			   0, '.', $cne->start1-1, $cne->end1,
			   join(',',$COLOR_INDEX->get_seq_color($cne->chr2))), "\n";
	}
    }

    # Cloe file
    close OUT;

    return $file_name;
}


sub write_wig_file
{
    my ($asm1_info, $asm2_info, $id, $len, $bed_fn) = @_;  
  
    my $asm1 = $asm1_info->{assembly_id};
    my $asm2 = $asm2_info->{assembly_id};

    # Make track name & desc strings
    my $org2 = ucfirst make_organism_string($asm2_info);
    my $track_name = "$org2 HCNE density ".($id*100)."% / $len";
    my $track_desc = "$org2 HCNE density (".($id*100)."% identity over $len columns; $WIN_SIZE kb sliding window)";

    # Get file name
    my $file_id = make_density_id_string($asm1, $asm2, $id, $len); 
    my $dir = get_dir_name($asm1_info, $asm2_info);
    my $wig_fn = "$dir/$file_id.wig";

    # Get 2bit file name
    my $twoBit_fn = "$ASM_DIR/$asm1/assembly.2bit";
    die "2bit file $twoBit_fn does not exist" unless(-e $twoBit_fn);

    # Get parameters for density calc
    die "Window size ${WIN_SIZE}000 bp is not evenly divisible by step size $STEP_SIZE bp" if(($WIN_SIZE*1000) % $STEP_SIZE);
    my $win_nr_steps = $WIN_SIZE * 1000 / $STEP_SIZE;

    # Run density computation script
    system("$DENSITY_SCRIPT -n '$track_name' -c $CLIP_PC -d '$track_desc' ".
	   "-p 'visibility=full maxHeightPixels=64:64:11' ".
	   "$twoBit_fn $bed_fn $STEP_SIZE $win_nr_steps >$wig_fn");

    return $wig_fn;
}


sub get_dir_name
{
    my($asm1_info, $asm2_info) = @_;

    # Get directory names
    my $asm1 = $asm1_info->{assembly_id};
    my $org2 = make_organism_string($asm2_info);
    my $dir1 = $asm1;
    my $dir2 = "$asm1/vs_$org2";
    $dir2 =~ s/ /_/g;

    # Create directories if they don't exist
    foreach my $dir ($dir1, $dir2) {
	unless(-e $dir) {
	    mkdir($dir) or die "could not create directory $dir";
	}
    }

    return $dir2;
}


sub make_hcne_id_string
{
    my ($asm1, $asm2, $id, $len) = @_;
    $id = $id * 100;
    return "HCNE_${asm1}_${asm2}_${id}pc_${len}col";
}


sub make_density_id_string
{
    my ($asm1, $asm2, $id, $len) = @_;
    $id = $id * 100;
    return "HCNE_density_${asm1}_${asm2}_${id}pc_${len}col";
}


sub make_organism_string
{
    my $asm_info = shift;
    return $asm_info->{organism_common} || $asm_info->{organism_latin};
}
