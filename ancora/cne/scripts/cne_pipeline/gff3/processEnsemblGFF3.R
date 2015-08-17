
library(Biostrings)
library(rtracklayer)
scripts <- list.files(path="/home/gtan/Repos/Bitbucket/project_specific/gtan/GRBEdges/scripts", pattern=".*\\.R", ignore.case=TRUE, full.names=TRUE)
for(script in scripts){message(script);source(script);}

# Set
gffFn <- "Danio_rerio.GRCz10.81.gff3"
fastaFn <- "danRer10.fa"

# fasta and gff
gff <- import.gff3(gffFn)
gff <- txdbENS2UCSC(gff)

gff <- gff[seqnames(gff) %in% names(fasta.seqlengths(fastaFn))]
gff <- gff[which(gff$type != "chromosome")]

export.gff3(gff, con=sub("\\gff3$", "gbrowse.gff3", gffFn))

