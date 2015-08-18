
library(rtracklayer)
library(Biostrings)

# fasta
## TODO: 
fastaFn <- "danRer10.fa"
seqlengths <- fasta.seqlengths(fastaFn)
## We only consider the assembled normal chroms
seqlengths <- seqlengths[grep("_", names(seqlengths), invert=TRUE)]

# CpG
mySession <- browserSession("UCSC")
## TODO:
genome(mySession) <- "danRer10"
cpg.grange <- GRanges(seqnames=names(seqlengths),
                      ranges=IRanges(start=1L,
                                     end=seqlengths))
tbl.cpg <- getTable(ucscTableQuery(mySession, track="rmsk",
                                   range=cpg.grange, table="rmsk"))
cpgGRangs <- GRanges(seqnames=tbl.cpg$genoName,
                     ranges=IRanges(start=tbl.cpg$genoStart+1L,
                                    end=tbl.cpg$genoEnd),
                     strand="+",
                     source=Rle("UCSC", nrow(tbl.cpg)),
                     type=Rle("rmsk", nrow(tbl.cpg)),
                     Name=tbl.cpg$repName
                     )
export.gff3(cpgGRangs, con="rmsk.gff3")


