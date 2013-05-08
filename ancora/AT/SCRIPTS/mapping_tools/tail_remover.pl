#!/usr/bin/perl -w
use strict;

# This script removes leading poly-T stretches and  trailing poly-A stretches
# from sequences. Sequences (fasta-format) are read from STDIN and written to
# STDOUT. Currently, the script takes no options.

use AT::Tools::PolyATailFinder;
use Bio::SeqIO;
use Bio::Seq;

my $seq_in = Bio::SeqIO->new(-fh => \*STDIN, -format => 'fasta');
my $seq_out = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'fasta');

while (my $seq = $seq_in->next_seq) {
    my ($tlen, $alen) = AT::Tools::PolyATailFinder->find_tail_in_seq($seq);
    my $trunc_seq = Bio::Seq->new
		(-id => $seq->id,
		 -description => $seq->description,
		 -seq => $seq->subseq($tlen+1, $seq->length-$alen));
    $seq_out->write_seq($trunc_seq);
}
