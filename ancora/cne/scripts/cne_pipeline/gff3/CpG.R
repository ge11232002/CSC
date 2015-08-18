
library(rtracklayer)
library(Biostrings)

# fasta
fastaFn <- "danRer10.fa"
seqlengths <- fasta.seqlengths(fastaFn)
## We only consider the assembled normal chroms
seqlengths <- seqlengths[grep("_", names(seqlengths), invert=TRUE)]

# CpG
mySession <- browserSession("UCSC")
genome(mySession) <- "danRer10"
cpg.grange <- GRanges(seqnames=names(seqlengths),
                      ranges=IRanges(start=1L,
                                     end=seqlengths))
tbl.cpg <- getTable(ucscTableQuery(mySession, track="cpgIslandExt",
                                   range=cpg.grange, table="cpgIslandExt"))
cpgGRangs <- GRanges(seqnames=tbl.cpg$chrom,
                     ranges=IRanges(start=tbl.cpg$chromStart+1L,
                                    end=tbl.cpg$chromEnd),
                     strand="+",
                     source=Rle("UCSC", nrow(tbl.cpg)),
                     type=Rle("CpG", nrow(tbl.cpg))
                     )
export.gff3(cpgGRangs, con="CpG.gff3")


