
library(Biostrings)
library(rtracklayer)
scripts <- list.files(path="/home/gtan/Repos/Bitbucket/project_specific/gtan/GRBEdges/scripts", pattern=".*\\.R", ignore.case=TRUE, full.names=TRUE)
for(script in scripts){message(script);source(script);}

# TODO:
gffFn <- "Danio_rerio.GRCz10.81.gff3"
fastaFn <- "danRer10.fa"

# fasta and gff
gff <- import.gff3(gffFn)
gff <- txdbENS2UCSC(gff)

gff <- gff[seqnames(gff) %in% names(fasta.seqlengths(fastaFn))]
gff <- gff[which(gff$type != "chromosome")]
gff$biotype <- NULL
gff$description <- NULL
gff$havana_gene <- NULL
gff$havana_version <- NULL
gff$logic_name <- NULL
gff$version <- NULL
gff$havana_transcript <- NULL
gff$gene_id <- NULL
gff$transcript_id <- NULL
gff$constitutive <- NULL
gff$ensembl_end_phase <- NULL
gff$ensembl_phase <- NULL
gff$exon_id <- NULL
gff$rank <- NULL
gff$protein_id <- NULL


export.gff3(gff, con=sub("\\gff3$", "gbrowse.gff3", gffFn))

