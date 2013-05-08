#!/usr/bin/perl -w

use Bio::SeqIO;

my $seqin = Bio::SeqIO->new(-fh => \*STDIN, -format => 'fasta');

while (my $seq = $seqin->next_seq) {
    my $first10 = $seq->subseq(1,10);
    my $last10 = $seq->subseq($seq->length-9, $seq->length);
    my @ts = $first10 =~ /t/gi;
    my @as = $last10 =~ /a/gi;
    print $seq->id, "\t";
    if(@ts >= 8 and @as < 8) {
	print 'T';
    }
    elsif(@as >= 8 and @ts < 8) {
	print 'A';
    }
    else {
	print '-';
    }
    print "\n";
}
