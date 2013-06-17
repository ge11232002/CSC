
library(rtracklayer)

build = "ce10"
ensemblGeneFile = "../annotation/ensembl_genes.txt"
assembly = file.path("/export/data/goldenpath", build, "assembly.2bit")

ensemblGene = read.table(ensemblGeneFile, sep="\t", header=TRUE, quote="", colClass="character")

bed = ensemblGene[, c("Chromosome.Name", "Exon.Chr.Start..bp.", "Exon.Chr.End..bp.")]

## for hg19, mm10, mm9, tetNig2, canFam3, galGal3, anoCar2, equCab2, ornAna1, ce10
bed = transform(bed, Chromosome.Name=paste("chr", Chromosome.Name, sep=""))
bed = transform(bed, Chromosome.Name=sub("chrMT", "chrM", Chromosome.Name))

## for petMar2
# do nothing

## for xenTro3
#bed = transform(bed, Chromosome.Name=sub("\\.1", "", Chromosome.Name))
#indexTransform = grep("^GL", bed[, "Chromosome.Name"], invert=TRUE)
#bed[indexTransform, "Chromosome.Name"] = paste("chr", bed[indexTransform, "Chromosome.Name"], sep="")
#bed = transform(bed, Chromosome.Name=sub("^chrMT$", "chrM", Chromosome.Name))

## for danRer7
#indexTransform = grep("^Zv", bed[, "Chromosome.Name"], invert=TRUE)
#bed[indexTransform, "Chromosome.Name"] = paste("chr", bed[indexTransform, "Chromosome.Name"], sep="")
#bed = transform(bed, Chromosome.Name=sub("^chrMT$", "chrM", Chromosome.Name))

## for fr2
#indexTransform = grep("^scaffold", bed[, "Chromosome.Name"], invert=TRUE)
#bed[indexTransform, "Chromosome.Name"] = paste("chr", bed[indexTransform, "Chromosome.Name"], sep="")
#bed = transform(bed, Chromosome.Name=sub("^chrMT$", "chrM", Chromosome.Name))

validChromsIndex = bed[, "Chromosome.Name"] %in% seqnames(seqinfo(TwoBitFile(assembly)))
message("before:", nrow(bed), "   after:", sum(validChromsIndex), " deleted:", nrow(bed)-sum(validChromsIndex))

bed = bed[validChromsIndex, ]
write.table(bed, file="ensembleExons.txt", quote=FALSE, sep="\t", row.names=FALSE,
            col.names=FALSE)


