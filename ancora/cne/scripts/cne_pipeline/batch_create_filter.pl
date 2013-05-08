use strict;

my @species = qw(
		 canFam2 
		 mm9
		 galGal3
		 xenTro2
		 danRer5
		 tetNig1
		 fr2
		 gasAcu1
		 oryLat1
		 ci2
		 );

system("perl ~/projects/cne/scripts/create_filter_set.pl -b HS_MAR06 -r -a /export/data/goldenpath/hg18/assembly.2bit > /export/data/CNEs/hg18/filters/repeat_regions.hg18.bed");
foreach my $s (@species){
    print "Processing $s..\n";
#    unless (-e  "${s}/filters/${s}_exon_regions.bed"){
#	warn("${s}/filters/${s}_exon_regions.bed does not exist..");
#	next;
#    }

#    system("perl ~/projects/cne/scripts/create_filter_set.pl -b ${s} -f /export/data/CNEs/${s}/filters/${s}_exon_regions.bed -r -a /export/data/goldenpath/${s}/assembly.2bit > /export/data/CNEs/${s}/filters/filter_regions.${s}.bed");
    system("perl ~/projects/cne/scripts/create_filter_set.pl -b ${s} -r -a /export/data/goldenpath/${s}/assembly.2bit > /export/data/CNEs/${s}/filters/repeat_regions.${s}.bed");
}
