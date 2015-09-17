# This script adds the chr length lines to gff3 for gbrowse
# In the gff3 file for gbrowse, you need a line for definition of chromosom name with Name="chr*" 
library(rtracklayer)
library(Biostrings)
library(GenomicRanges)

## TODO:
gffFn <- "refGene.danRer10.gff3"
fastaFn <- "danRer10.fa"

## fasta and gff
gff <- import.gff3(gffFn)

fastaGRanges <- GRanges(seqnames=names(fasta.seqlengths(fastaFn)),
                        ranges=IRanges(start=1L,
                                       end=fasta.seqlengths(fastaFn)),
                        strand="*",
                        Name=names(fasta.seqlengths(fastaFn))
                        )
stopifnot(all(seqnames(gff) %in% seqnames(fastaGRanges)))
export.gff3(fastaGRanges, con="chroms.gff3")




