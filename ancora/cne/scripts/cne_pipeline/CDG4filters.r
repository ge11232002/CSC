
library(rtracklayer)
CDGGff = list.files("../annotation/", "*.gff", full.names=TRUE)
stopifnot(length(CDGGff) == 1)

filtersFn = "filter_regions.CBS51388_real.bed"

chroms = seqnames(seqinfo(TwoBitFile("/export/data/goldenpath/CBS51388/assembly.2bit")))

gff = read.table(CDGGff, header=FALSE, sep="\t", stringsAsFactors=FALSE, quote="")

filterFeatures = c("exon", "repeat_region", "long_terminal_repeat")

gff = gff[gff$V3 %in% filterFeatures, c("V1", "V4", "V5")]

gff = transform(gff, V4=V4-1)
dim(gff)
gff = gff[gff$V1 %in% chroms, ]
dim(gff)

write.table(gff, file=filtersFn, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

