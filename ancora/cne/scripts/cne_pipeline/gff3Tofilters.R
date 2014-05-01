
library(rtracklayer)

## set
gffFn <- "../annotation/Glossina-morsitans-Yale_BASEFEATURES_GmorY1.3.gff3_0"
filtersFn <- "filter_regions.GmorY1.bed"
chroms <- seqnames(seqinfo(TwoBitFile("/export/data/goldenpath/GmorY1/assembly.2bit")))

## compute
gff <- read.table(gffFn, header=FALSE, sep="\t", stringsAsFactors=FALSE, quote="")
filterFeatures = c("exon", "repeat_region", "long_terminal_repeat")

gff = gff[gff$V3 %in% filterFeatures, c("V1", "V4", "V5")]

gff = transform(gff, V4=V4-1)
dim(gff)
gff = gff[gff$V1 %in% chroms, ]
dim(gff)

write.table(gff, file=filtersFn, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

