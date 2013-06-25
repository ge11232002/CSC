

readAxt = function(axtFile){
  axtContent = readLines(axtFile)
  axtContent = axtContent[grep("^#", axtContent, invert=TRUE)]
  axtContent = axtContent[axtContent != ""]
  axtContent = split(axtContent, c("header", "target", "query"))
  header = strsplit(axtContent$header, split=" ", fixed=TRUE)
  myAxt = axt(targetNames=Rle(sapply(header, "[",2)),
              targetRanges=IRanges(start=as.integer(sapply(header, "[",3)), end=as.integer(sapply(header, "[",4)), names=sapply(header, "[",1)),
              targetSeqs=DNAStringSet(axtContent$target),
              queryNames=Rle(sapply(header, "[",5)),
              queryRanges=IRanges(start=as.integer(sapply(header, "[",6)), end=as.integer(sapply(header, "[",7)), names=sapply(header, "[",1)),
              querySeqs=DNAStringSet(axtContent$query),
              strand=Rle(sapply(header, "[",8)),
              score=Rle(as.integer(sapply(header, "[",9)))
              )
  myAxt2 = axt2(targetRanges=GRanges(seqnames=Rle(sapply(header, "[",2)),
                                     ranges=IRanges(start=as.integer(sapply(header, "[",3)), end=as.integer(sapply(header, "[",4)), names=sapply(header, "[",1)),
                                     strand=Rle("+")),
                targetSeqs=DNAStringSet(axtContent$target),
                queryRanges=GRanges(seqnames=Rle(sapply(header, "[",5)),
                                    ranges=IRanges(start=as.integer(sapply(header, "[",6)), end=as.integer(sapply(header, "[",7)), names=sapply(header, "[",1)),
                                    strand=Rle(sapply(header, "[",8))),
                querySeqs=DNAStringSet(axtContent$query),
                score=as.integer(sapply(header, "[",9))
                )
  rm(axtContent)
  rm(header)
  gc()
  return(myAxt)
}

# tFilterFile = "/export/data/CNEs/danRer7/filters/filter_regions.danRer7.bed"
# qFilterFile = "/export/data/CNEs/hg19/filters/filter_regions.hg19.bed"
readBed = function(bedFile){
## This GRanges is in 1-based.
  require(rtracklayer)
  bed = import(bedFile, asRangedData = FALSE)
  strand(bed) = "+"
  bed = reduce(bed)
}
# system.time(foo<-readBed(bedFile))
#    user  system elapsed
#    51.890   0.757  52.659

readBedToGRanges = function(bedFile){
  require(GenomicRanges)
  dyn.load("~/Repos/CSC/CNEr/src/io.so")
  bed = .Call("myReadBed", bedFile)
  bed = GRanges(seqnames=Rle(bed[[1]]),
                ranges=IRanges(start=bed[[2]], end=bed[[3]]),
                strand=factor("+"))
  return(bed)
}
# system.time(foo<-readBedToGRanges(bedFile))
#    user  system elapsed
#    2.272   0.133   2.414



