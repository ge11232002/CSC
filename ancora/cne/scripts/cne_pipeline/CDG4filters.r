

CDGGff = "../annotation/C_parapsilosis_CDC317_version_s01-m03-r07_features.gff"
filtersFn = "filter_regions.C_parapsilosis_CDC317.bed"

gff = read.table(CDGGff, header=FALSE, sep="\t", stringsAsFactors=FALSE)

filterFeatures = c("exon", "repeat_region", "long_terminal_repeat")

gff = gff[gff$V3 %in% filterFeatures, c("V1", "V4", "V5")]

write.table(gff, file=filtersFn, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

