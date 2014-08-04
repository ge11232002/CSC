## This script is used to mask the reference genome sequence with the repeats coordinates from gff file downloaded from http://www.candidagenome.org/download/gff/C_albicans_SC5314/

gffFn = "C_albicans_SC5314_version_A22-s02-m01-r05_features.gff"

library(rtracklayer)
library(GenomicRanges)
gff = read.table(gffFn, sep="\t", stringsAsFactors=FALSE)
repeats = subset(gff, V3 %in% c("repeat_region", "long_terminal_repeat"), 
                 select=c("V1", "V4", "V5", "V7"))
colnames(repeats) = c("chrom", "start", "end", "strand")
repeats <- GRanges(seqnames=repeats$chrom,
                   ranges=IRanges(start=repeats$start,
                                  end=repeats$end),
                   strand=repeats$strand)
export.bed(repeats, con="repeats.bed")

## cmd
### twoBitMask -type=.bed SC5314A22.unmasked.2bit /export/data/CNEs/SC5314A22/annotation/repeats.bed SC5314A22.masked.2bit

