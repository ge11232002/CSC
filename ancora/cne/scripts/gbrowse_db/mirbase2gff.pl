use warnings;
use strict;
use Class::Struct;
use AT::DB::GenomeAssemblyTwoBit;
use AT::Tools::RangeHandler;
use Getopt::Std;
use Data::Dumper;

my ($TWOBIT_FN, $MIRBASE_FN) = @ARGV;

unless($MIRBASE_FN){

    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Usage: 
    
$cmd <2bit file> <miRBase gff file>

<2bit file>         Name of .2bit file for assembly
<miRBase gff file>  File to convert

convert the gff3 file to the gff2 file which is used by ancora

Output: ancora gbrowse gff file (on standard output)

EndOfUsage

    exit;
}

# Get list of valid chromsomes from 2bit file
my $asm = AT::DB::GenomeAssemblyTwoBit->new(file => $TWOBIT_FN)
  or die "could not open assembly file $TWOBIT_FN";
my %valid_chr = map { $_ => 1 } $asm->get_chr_names;
my %skipped_chr;

open IN, $MIRBASE_FN or die "could not open $MIRBASE_FN";

my %mirnas_by_id;
while(my $line = <IN>) {
    next if($line =~ /^#/);
    chomp $line;
    my @fields = split /\t/, $line;
    
    # Change chr name to UCSC-style and check that it is valid
    my $chr = $fields[0];
    #$chr = 'M' if($chr eq 'MT');
    $chr = 'chrM' if($chr eq 'chrMT');
    unless($valid_chr{$chr}) {
	$skipped_chr{$chr} = 1;
	next;
    }
    #$fields[0] = $chr = "chr$chr";
    $fields[0] = $chr;
    # Parse group field
    my $group = $fields[8];
    #my ($id) = $group =~ /; ID="([^"]+)"/;
    my ($id) = $group =~ /Name=(.+)/;
    $id =~ s/;.+//;
    die "could not parse group field: $group" unless($id);

    push @{$mirnas_by_id{$id}}, \@fields;
}

while (my ($id, $fields_list) = each %mirnas_by_id) {
    my $n = @$fields_list;
    my $i = 1;
    foreach my $fields (@$fields_list) {
	    my $group = $n == 1 ? "Gene \"$id\"" : "Gene \"$id mapping $i\"";
	    if($id =~ /^...-/) {
	        $group .= '; Alias "'.substr($id,4).'"';
	    }
	    else { warn "Could not parse short name from $id\n"} 
	    $fields->[8] = $group;
	    # output gff
	    print join("\t", @$fields), "\n";
	    $i++;
    }
}
