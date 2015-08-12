## This script extracts the repeats coordinates from the ENSEMBL softmasked sequences.
library(Biostrings)
library(GenomicRanges)
library(rtracklayer)

foo = readBStringSet("/export/data/goldenpath/droMoj28/droMoj28.fa")
names(foo) = sapply(strsplit(names(foo), " "), "[", 1)
foo3 = lapply(lapply(strsplit(as.character(foo),"") ,"%in%", 
                     c("a","c","g","t")), Rle)
foo4 = lapply(foo3, as, "IRanges")
foo5 = GRanges(seqnames=Rle(names(foo4), sapply(foo4, length)),
    ranges=IRanges(start=unlist(sapply(foo4, start)),
             end=unlist(sapply(foo4, end))),
    strand="+")

export.bed(foo5, "repeats.bed")


