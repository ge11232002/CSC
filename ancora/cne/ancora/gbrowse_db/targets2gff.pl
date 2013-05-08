use warnings;
use strict;

my ($TARGET_FN) = @ARGV;

unless($TARGET_FN){

    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Usage: 
    
perl ucsc2gff.pl -d <db name> refSeq | $cmd <target bed file>

<db name>          Name of UCSC database   
<target bed file>  4-column bed file listing target genes

Output: ancora gbrowse gff file (on standard output)

EndOfUsage

    exit;
}

my %targets;
my $i = 1;
open TARGETS, $TARGET_FN or die "Could not open $TARGET_FN"; 
while(my $line = <TARGETS>) {
    chomp $line;
    my ($chr, $start, $end, $symbol) = split '\s+', $line;
    next unless($chr);
    $start++;
    $targets{$symbol} = [$chr, $start, $end, $i];
    $i++;
}

my %targets_found;
while(my $line = <STDIN>) {
    chomp $line;
    my ($chr, $source, $type, $start, $end, undef, $strand, undef, $group) = split /\t/, $line;
    next unless($type eq 'gene');
    my ($symbol) = $group =~ /Alias "([^"]+)"/;
    next unless($symbol);
    my $target = $targets{$symbol};
    next unless($target);
    my ($t_chr, $t_start, $t_end, $i) = @$target;
    if($chr ne $t_chr) {
	#warn "Gene $symbol is on $chr in database and on $t_chr in target file";
	next;
    }
    my $my_group = "Gene target$i ; Alias \"$symbol\" ; Note \"Putative HCNE target gene\"";
    print join("\t", $chr, 'Target', 'gene', $start, $end, '.', $strand, '.', $my_group), "\n";
    $targets_found{$symbol}++;
}

foreach my $symbol (keys %targets) {
    unless($targets_found{$symbol}) {
	warn "Target gene $symbol not found.\n";
    }
}
