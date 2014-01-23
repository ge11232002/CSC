## This script is used to mask the reference genome sequence with the repeats coordinates from gff file downloaded from http://www.candidagenome.org/download/gff/C_albicans_SC5314/

gffFn = "/export/data/CNEs/CBS51388/annotation/A_niger_CBS_513_88_version_s01-m06-r09_features.gff"
twoBitFn = "/export/data/goldenpath/CBS51388/assembly.2bit"

library(rtracklayer)
library(GenomicRanges)
reference = import.2bit(twoBitFn)
gff = read.table(gffFn, sep="\t", stringsAsFactors=FALSE)
chromSizes = subset(gff, V3 %in% c("contig", "chromosome") , select=c("V1", "V4", "V5"))
colnames(chromSizes) = c("chrom", "start", "end")
rownames(chromSizes) = chromSizes$chrom

repeats = subset(gff, V3 %in% c("repeat_region", "long_terminal_repeat"), 
                 select=c("V1", "V4", "V5", "V7"))
colnames(repeats) = c("chrom", "start", "end", "strand")
if(length(repeats) == 0L)
  stop("no need to proceed")

## process the coordinates of repeats on the negative strand
indexNegative = repeats$strand == "-"
newStarts = chromSizes[repeats$chrom[indexNegative], "end"] - repeats$end[indexNegative] + 1
newEnds = chromSizes[repeats$chrom[indexNegative], "end"] - repeats$start[indexNegative] + 1
repeats$start[indexNegative] = newStarts
repeats$end[indexNegative] = newEnds
repeats$strand = "+"
repeats = GRanges(seqnames=repeats$chrom, 
                  ranges=IRanges(start=repeats$start, end=repeats$end),
                  strand=repeats$strand)
repeats = reduce(repeats)

## mask the reference seq
referenceNames = names(reference)
reference = as.character(reference) 
oneRepeat = repeats[1]
for(i in 1:length(repeats)){
  oneRepeat = repeats[i]
  subseq(reference[as.character(seqnames(oneRepeat))], 
         start=start(oneRepeat), end=end(oneRepeat)) = 
  tolower(subseq(reference[as.character(seqnames(oneRepeat))],
                 start=start(oneRepeat), end=end(oneRepeat)))
}
library(seqinr)
reference = strsplit(reference, "")
write.fasta(reference, names=names(reference), file.out=sub("2bit$", "fa", twoBitFn))

## make it back to 2bit
unlink(twoBitFn)
cmd = paste("faToTwoBit", sub("2bit$", "fa", twoBitFn), twoBitFn)
system(cmd)
#unlink(sub("2bit$", "fa", twoBitFn))

