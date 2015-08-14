# This script adds the chr length lines to gff3 for gbrowse
# In the gff3 file for gbrowse, you need a line for definition of chromosom name with Name="chr*" 
library(rtracklayer)
library(Biostrings)
library(GenomicRanges)

gffFn <- "/var/lib/gbrowse/databases/danRer10/refGene.danRer10.gff"
fastaFn <- "/var/lib/gbrowse/databases/danRer10_fasta/danRer10.fa"

### fasta
fasta <- readDNAStringSet(fastaFn)
gff <- import.gff(gffFn)

fastaGRanges <- GRanges(seqnames=names(fasta),
                        ranges=IRanges(start=1L,
                                       end=width(fasta)),
                        strand="*",
                        Name=names(fasta)
                        )
stopifnot(all(seqnames(gff) %in% seqnames(fastaGRanges)))
export.gff3(fastaGRanges, con="chroms.gff")




