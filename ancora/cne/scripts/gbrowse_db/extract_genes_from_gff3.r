selfDir = "~/Repositories/genomics"
selfScripts = list.files(path=selfDir, pattern='.*\\.r', full.names=TRUE, recursive=TRUE)
for(rs in selfScripts){message(rs);source(rs)}

################start of Read the gff3 from flybase and make the gff for gbrowse
get_parents = function(extras, filed=NULL){
  parentAttr = my.getGffAttributeField(extras, "Parent")
  parentsArray = strsplit(parentAttr, split=",")[[1]]
  return(parentsArray)
}
read_flybase_gff = function(gff_fn){
  gff = my.gffRead(gff_fn)
  gff = gff[gff[ ,"source"]=="FlyBase", ]
  indexType = gff[, "type"] %in% c("gene", "CDS", "exon", "pseudogene")
  indexType = indexType | grepl("RNA$", gff[, "type"])
  gff = gff[indexType, ]
  gc()
  gff$ID = my.getGffAttributeField(gff$attributes, "ID")
  gff$Name = my.getGffAttributeField(gff$attributes, "Name")
  gff$Parent = my.getGffAttributeField(gff$attributes, "Parent")
  Parent2 = my.getGffAttributeField(gff$attributes, "Derives_from")
  gff$Parent[is.na(gff$Parent)] = Parent2[is.na(gff$Parent)]
  rm(Parent2)
  gc()
## Handle gene feature
  indexGene = gff[ ,"type"] == "gene"
  genes = gff[indexGene, c("seqid", "start", "end", "strand", "Name")]
  colnames(genes) = c("chr", "start", "end", "strand", "gene_symbol")
  rownames(genes) = gff$ID[indexGene]
## Handle CDS, exon, transcript / pseudogene feature (coding or noncoding)
## The regex is supposed to capture: mRNA, tRNA, ncRNA, rRNA, miRNA, snoRNA, scRNA, snoRNA, snRNA
  indexTranscripts = grepl("RNA$", gff[, "type"]) | gff[, "type"] == "pseudogene"
  transcripts = data.frame(row.names=gff$ID[indexTranscripts])
  colnames(transcripts) = c("chr", "type", "start", "end", "strand", "transcript_symbol", "gene_id")
  transcripts[ , c("chr", "type", "start", "end", "strand", "transcript_symbol", "gene_id")] = gff[indexTranscripts, c("seqid", "type", "start", "end", "strand", "Name", "Parent")]
  # stop if any transcript without gene
  stopifnot(!any(is.na(transcripts[ ,"gene_id"])))
  # stop if any transcript has more than one gene parent.
  stopifnot(!any(sapply(strsplit(transcripts[ ,"gene_id"], ","), length) > 1))
# CDS
  indexCDS = gff[ ,"type"] == "CDS"
  parentsCDS = strsplit(gff[indexCDS, "Parent"], ",")
  startCDS = sapply(split(rep(gff[indexCDS, "start"], sapply(parentsCDS, length)), unlist(parentsCDS)), min)
  endCDS = sapply(split(rep(gff[indexCDS, "end"], sapply(parentsCDS, length)), unlist(parentsCDS)), max)
  transcripts[ , c("cds_start", "cds_end")] = cbind(startCDS[rownames(transcripts)], endCDS[rownames(transcripts)])
  rm(startCDS, endCDS, parentsCDS)
# exon
  indexExon = gff[, "type"] == "exon"
  parentsExon = strsplit(gff[indexExon, "Parent"], ",")
  startExon = split(rep(gff[indexExon, "start"], sapply(parentsExon, length)), unlist(parentsExon))
  #startExon = sapply(startExon, paste, collapse=",")
  endExon = split(rep(gff[indexExon, "end"], sapply(parentsExon, length)), unlist(parentsExon))
  #endExon = sapply(endExon, paste, collapse=",")
  transcripts[ , c("exons_start", "exons_end")] = cbind(startExon[rownames(transcripts)], endExon[rownames(transcripts)])
  transcripts$exons_start = startExon[rownames(transcripts)]
  transcripts$exons_end   = endExon[rownames(transcripts)]
  return(list("genes"=genes, "transcripts"=transcripts))
}

validate_transcripts = function(transcripts){
  require(GenomicRanges)
  message("Validating transcripts...")
  indexWarning = !is.na(transcripts[ , "cds_start"]) & transcripts[ , "type"] != "mRNA"
  if(any(indexWarning)){
    warning("The following transcripts have CDS coordinates but are not classified as mRNA", rownames(transcripts)[indexWarning])
  }
  indexWarning = is.na(transcripts[ , "cds_start"]) & transcripts[ , "type"] == "mRNA"
  if(any(indexWarning)){
    warning("The following transcripts do not have CDS coordinates but are classified as mRNA", rownames(transcripts)[indexWarning])
  }
## collapse the overlapping exons
  expandNumber = sapply(transcripts[ , "exons_start"], length)
  # the GRanges use the strand "+" "-" "*"
  collapsed = reduce(GRanges(seqnames=rep(rownames(transcripts), expandNumber),
          ranges=IRanges(start=unlist(transcripts[ , "exons_start"]),
                         end=unlist(transcripts[ ,"exons_end"])),
          strand=rep(sub("\\.", "*", transcripts[ , "strand"]), expandNumber)))
  transcripts$exons_start = split(start(collapsed), as.vector(seqnames(collapsed)))[rownames(transcripts)]
  transcripts$exons_end   = split(end(collapsed), as.vector(seqnames(collapsed)))[rownames(transcripts)]
## check the exons and transcripts boundaries
  indexExon = !sapply(transcripts$exons_start, is.null)
  indexWarning = transcripts$start[indexExon] != sapply(transcripts$exons_start[indexExon], min) | transcripts$end[indexExon] != sapply(transcripts$exons_end[indexExon], max)
  if(any(indexWarning)){
    warning("The following transcripts' boundaries do not match exons ", rownames(transcripts)[indexExon][indexWarning])
  }
  return(transcripts)
}

write_gbrowse_gff = function(genes, transcripts, source="FlyBase", output="output.gff"){
  message("Start writing the gbrowse gff genes")
  genes$group = paste("Gene ", rownames(genes), sep="")
  indexGeneSymbol = !is.na(genes$gene_symbol)
  genes$group[indexGeneSymbol] = paste(genes$group[indexGeneSymbol],
                                       " ; Alias \"", genes$gene_symbol[indexGeneSymbol], "\"", sep="")
  genes$group = paste(genes$group, " ; Note \"", source, " Gene\"", sep="")
  genes$chr = paste("chr", genes$chr, sep="")
  outputGff = cbind(genes$chr, source, "gene", genes$start, genes$end, ".", genes$strand, ".", genes$group)
  write.table(outputGff, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)

  message("Start writing the gbrowse gff transcripts")
  transcripts$chr = paste("chr", transcripts$chr, sep="")
  transcripts$group = paste("Transcript ", rownames(transcripts), sep="")
  indexTranscriptsSymbol = !is.na(transcripts$transcript_symbol)
  transcripts$group[indexTranscriptsSymbol] = paste(transcripts$group[indexTranscriptsSymbol],
                                    " ; Symbol \"", transcripts$transcript_symbol[indexTranscriptsSymbol], "\"", sep="")
## write the mRNA
  indexmRNA = transcripts$type == "mRNA"
  outputGff = cbind(transcripts[indexmRNA, "chr"], source, "mRNA",
                    transcripts[indexmRNA, "start"], transcripts[indexmRNA, "end"],
                    ".", transcripts[indexmRNA, "strand"], ".",
                    transcripts[indexmRNA, "group"])
  write.table(outputGff, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
## write CDS and UTR
  indexCDS = !is.na(transcripts$cds_start)
  expandNumber = sapply(transcripts[indexCDS, "exons_start"], length)
  exonsGRanges = GRanges(seqnames=rep(rownames(transcripts)[indexCDS], expandNumber),
          ranges=IRanges(start=unlist(transcripts[indexCDS, "exons_start"]),
                         end=unlist(transcripts[indexCDS,"exons_end"])),
          strand=rep(sub("\\.", "*", transcripts[indexCDS, "strand"]), expandNumber))
  cdsGRanges = GRanges(seqnames=rownames(transcripts)[indexCDS],
                       ranges=IRanges(start=transcripts[indexCDS, "cds_start"],
                                      end=transcripts[indexCDS, "cds_end"]),
                       strand=sub("\\.", "*", transcripts[indexCDS, "strand"]))
  fputrGRanges = GRanges(seqnames=rownames(transcripts)[indexCDS],
                         ranges=IRanges(start=transcripts[indexCDS, "start"],
                                        end=transcripts[indexCDS, "cds_start"]-1),
                         strand=sub("\\.", "*", transcripts[indexCDS, "strand"]))
  tputrGranges = GRanges(seqnames=rownames(transcripts)[indexCDS],
                    #some cds_end is larger than end. Is it possible?
                         ranges=IRanges(start=pmin(transcripts[indexCDS, "cds_end"]+1, transcripts[indexCDS, "end"]),
                                        end=transcripts[indexCDS, "end"]),
                         strand=sub("\\.", "*", transcripts[indexCDS, "strand"]))
  cdsGRangesIntersect = intersect(exonsGRanges, cdsGRanges)
  fputrGRangesIntersect = intersect(exonsGRanges, fputrGRanges)
  tputrGrangesIntersect = intersect(exonsGRanges, tputrGranges)
  fputrGRangesAdjusted = c(fputrGRangesIntersect[strand(fputrGRangesIntersect) != "-"],
                           tputrGrangesIntersect[strand(tputrGrangesIntersect) == "-"])
  tputrGrangesAdjusted =c(tputrGrangesIntersect[strand(tputrGrangesIntersect) != "-"],
                          fputrGRangesIntersect[strand(fputrGRangesIntersect) == "-"])
  rm(fputrGRangesIntersect, tputrGrangesIntersect)
# write CDS
  outputGff = cbind(transcripts[as.vector(seqnames(cdsGRangesIntersect)), "chr"], source,
                    "CDS", start(cdsGRangesIntersect), end(cdsGRangesIntersect),
                    ".", as.vector(strand(cdsGRangesIntersect)), ".",
                    transcripts[as.vector(seqnames(cdsGRangesIntersect)), "group"]
                    )
  write.table(outputGff, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
# write fputr
  outputGff = cbind(transcripts[as.vector(seqnames(fputrGRangesAdjusted)), "chr"], source,
                    "5'-UTR", start(fputrGRangesAdjusted), end(fputrGRangesAdjusted),
                    ".", as.vector(strand(fputrGRangesAdjusted)), ".",
                    transcripts[as.vector(seqnames(fputrGRangesAdjusted)), "group"]
                    )
  write.table(outputGff, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
# write tputr
  outputGff = cbind(transcripts[as.vector(seqnames(tputrGrangesAdjusted)), "chr"], source,
                    "3'-UTR", start(tputrGrangesAdjusted), end(tputrGrangesAdjusted),
                    ".", as.vector(strand(tputrGrangesAdjusted)), ".",
                    transcripts[as.vector(seqnames(tputrGrangesAdjusted)), "group"]
                    )
  write.table(outputGff, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
## write the other transcripts and exons
  indexNOCDS = is.na(transcripts$cds_start)
# write transcipt
  outputGff = cbind(transcripts[indexNOCDS, "chr"], source, "transcript",
                    transcripts[indexNOCDS, "start"], transcripts[indexNOCDS, "end"],
                    ".", transcripts[indexNOCDS, "strand"], ".",
                    transcripts[indexNOCDS, "group"])
  write.table(outputGff, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
# write the exons
  transciptsSubset = transcripts[indexNOCDS, ]
  transciptsSubset = transciptsSubset[!sapply(transciptsSubset$exons_start, is.null), ]
  expandNumber = sapply(transciptsSubset[ , "exons_start"], length)
  outputGff = cbind(rep(transciptsSubset[, "chr"], expandNumber), source, "exon",
                    unlist(transciptsSubset[ , "exons_start"]),
                    unlist(transciptsSubset[ , "exons_end"]),
                    ".", rep(transciptsSubset[ , "strand"], expandNumber),
                    ".", rep(transciptsSubset[ , "group"], expandNumber)
                    )
  write.table(outputGff, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
}

################## end of flybase gff


################## start of wormbase gff
read_wormbase_gff = function(gff_fn, genesID_fn){
  gff = my.gffRead(gff_fn)
  genesID = read.table(genesID_fn, header=FALSE, sep=",", quote="", as.is=TRUE)
  colnames(genesID) = c("gene_id", "symbol", "alias")
  gff = gff[gff$source %in% c("Coding_transcript", "Non_coding_transcript", "Pseudogene") | grepl("RNA$", gff$source), ]
  gc()
  gff$ID = my.getGffAttributeField(gff$attributes, "ID")
  gff$Parent = my.getGffAttributeField(gff$attributes, "Parent")
## Handle gene / pseudogene feature
  indexGene = gff$type %in% c("gene", "Pseudogene")
  genes = gff[indexGene, c("seqid", "start", "end", "strand", "ID")]
  genes$gene_symbol = sub("^Gene:", "", genes$ID)
  symbolCorrection = grep("^WBGene", genes$gene_symbol)
  genes$gene_symbol[symbolCorrection] = genesID[match(genes$gene_symbol[symbolCorrection], genesID$gene_id), "alias"]
  genes$gene_id = genesID[match(genes$gene_symbol, genesID$alias), "gene_id"]
  genes$ID = NULL
  colnames(genes) = c("chr", "start", "end", "strand", "gene_symbol", "gene_id")
  rownames(genes) = genes$gene_id
  genes$gene_id = NULL
  genes$chr = sub("^CHROMOSOME_", "", genes$chr)
## Handle transcript feature (coding or noncoding)
  indexTranscripts = gff$type == "mRNA" | gff$type == "ncRNA"
  transcripts = data.frame(row.names=sub("^Transcript:", "", gff$ID[indexTranscripts]))
  transcripts[ , c("chr", "type", "start", "end", "strand", "gene_id")] = gff[indexTranscripts, c("seqid", "type", "start", "end", "strand", "Parent")]
  transcripts$gene_id = sub("^Gene:", "", transcripts$gene_id)
  idCorrection = grep("^WBGene", transcripts$gene_id, invert=TRUE)
  transcripts$gene_id[idCorrection] = genesID[match(transcripts$gene_id[idCorrection][idCorrection], genesID$alias), "gene_id"]
  transcripts$transcript_symbol = genesID[match(transcripts$gene_id, genesID$gene_id), "symbol"]
  # warning if any transcript without gene
  warning(any(is.na(transcripts[ ,"gene_id"])))
  transcripts = transcripts[!is.na(transcripts$gene_id), ]
  # stop if any transcript has more than one gene parent.
  stopifnot(!any(sapply(strsplit(transcripts[ ,"gene_id"], ","), length) > 1))
## CDS
  indexCDS = gff[ ,"type"] == "CDS"
  CDS = gff[indexCDS, ]
  CDS$Parent = gsub("Transcript:", "", CDS$Parent)
  parentsCDS = strsplit(CDS$Parent, ",")
  startCDS = sapply(split(rep(CDS$start, sapply(parentsCDS, length)), unlist(parentsCDS)), min)
  endCDS = sapply(split(rep(CDS$end, sapply(parentsCDS, length)), unlist(parentsCDS)), max)
  transcripts[ , c("cds_start", "cds_end")] = cbind(startCDS[rownames(transcripts)], endCDS[rownames(transcripts)])
  rm(startCDS, endCDS, parentsCDS, CDS)
## exon
  indexExon = gff[, "type"] == "exon" | gff[, "type"] == "five_prime_UTR" | gff[, "type"] == "three_prime_UTR"
  exons = gff[indexExon, ]
  exons$Parent = gsub("Transcript:", "", exons$Parent)
  parentsExon = strsplit(exons$Parent, ",")
  startExon = split(rep(exons$start, sapply(parentsExon, length)), unlist(parentsExon))
  endExon = split(rep(exons$end, sapply(parentsExon, length)), unlist(parentsExon))
  transcripts[ , c("exons_start", "exons_end")] = cbind(startExon[rownames(transcripts)], endExon[rownames(transcripts)])
  transcripts$exons_start = startExon[rownames(transcripts)]
  transcripts$exons_end   = endExon[rownames(transcripts)]
  transcripts$chr = sub("^CHROMOSOME_", "", transcripts$chr)
  return(list("genes"=genes, "transcripts"=transcripts))
}




## Implementation

# FlyBase
gff_fn = "/export/data/CNEs/dm3/annotation/dmel-all-r5.51.gff"
gff = read_flybase_gff(gff_fn)
gff$transcripts = validate_transcripts(gff$transcripts)
write_gbrowse_gff(gff$genes, gff$transcripts, source="FlyBase", output="flybase.gff")

# WormBase
gff_fn = "/export/data/CNEs/ce10/annotation/c_elegans.WS220.annotations.gff3"
genesID_fn = "/export/data/CNEs/ce10/annotation/geneIDs.WS220"

gff = read_wormbase_gff(gff_fn, genesID_fn)
gff$transcripts = validate_transcripts(gff$transcripts)

write_gbrowse_gff(gff$genes, gff$transcripts, source="WormBase", output="wormbase.gff")

