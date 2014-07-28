# Whole genome pairwise alignment with last
This documentation describes the pipeline of running genome pairwise alignment
with the aligner called [last](http://last.cbrc.jp/).
This pipeline is not the default one used by UCSC. 
However, the advantage of using thiis `last` aligner is that
it can deal with genome which is not assembled quite well (in tens of thousands of contigs/scaffolds).

## last configuration and usage
This program finds local alignments between query sequences, and
reference sequences that have been prepared using `lastdb`.

```sh
lastdb humanDb humanChromosome*.fasta
lastal humanDb dna*.fasta > myalns.maf
# or pipe
zcat seqs.fasta.gz | lastal humanDb > myalns.maf
```

## run pipeline
follow the pipeline `lastPipeline.R`.



