use warnings;
use strict;
use Class::Struct;
use AT::DB::GenomeAssemblyTwoBit;
use AT::Tools::RangeHandler;
use Getopt::Std;


my ($TWOBIT_FN, $REDFLY_FN, $FT_TYPE) = @ARGV;

unless($FT_TYPE){

    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Usage: 
    
$cmd <2bit file> <REDfly gff file> <type>

<2bit file>         Name of .2bit file for assembly
<REDfly gff file>   File to convert, in GFF v.3 format
<type>              Feature type, usually CRM or TFBS

Output: ancora gbrowse gff file (on standard output)

EndOfUsage

    exit;
}

# Get list of valid chromsomes from 2bit file
my $asm = AT::DB::GenomeAssemblyTwoBit->new(file => $TWOBIT_FN)
  or die "could not open assembly file $TWOBIT_FN";
my %valid_chr = map { $_ => 1 } $asm->get_chr_names;
my %skipped_chr;

open IN, $REDFLY_FN or die "could not open $REDFLY_FN";

my %features_by_id;

while(my $line = <IN>) {
    next if($line =~ /^#/);
    chomp $line;
    my @fields = split /\t/, $line;
    
    # Change chr name to UCSC-style and check that it is valid
    my $chr = $fields[0];
    $chr =~ s/^\s+//; 
    $chr =~ s/\s+$//; 
    $chr = 'M' if($chr eq 'MT');
    unless($valid_chr{"chr$chr"}) {
	$skipped_chr{$chr} = 1;
	next;
    }
    $fields[0] = $chr = "chr$chr";
    $fields[1] = "REDfly";
    $fields[2] = $FT_TYPE;

    # Parse group field
    my $group = $fields[8];
    my ($id) = $group =~ /ID="([^"]+)"/;
    print STDERR "The id is ", $id, "\n";
    #my ($redfly_internal_id) = $group =~ /"REDfly=([\d+]+)[",]/;
    my ($redfly_internal_id) = $group =~ /REDfly="([^"]+)"/;
    print STDERR "The internal id is ", $redfly_internal_id, "\n";
    die "could not parse group field: $group" unless($id and $redfly_internal_id);

    push @{$features_by_id{$id}}, [$redfly_internal_id, @fields];
}

print STDERR "Skipped genes on sequences: ".join(", ", keys %skipped_chr)." beacuse they are not in the 2bit file.\n"
  if(scalar keys %skipped_chr);

while (my ($id, $fields_list) = each %features_by_id) {
    my $n = @$fields_list;
    my $i = 1;
    foreach my $fields (@$fields_list) {
	#my $group = $n == 1 ? "Gene \"$id\"" : "Gene \"$id mapping $i\"";
	my $redfly_internal_id = shift @$fields;
	my $group = "REDfly \"$id\" ; REDflyID $redfly_internal_id";
	$fields->[8] = $group;
	# output gff
	print join("\t", @$fields), "\n";
	$i++;
    }
}

# Make group as REFfly symbol; ID id
# source=REDfly; method = TFBS or CRM
