

readAxt_obsolete = function(axtFile){
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

#############################C version######################################
readBedToGRanges = function(bedFile=NULL){
  require(GenomicRanges)
  if(is.null(bedFile)){
    return(NULL)
  }
  if(!file.exists(bedFile)){
    stop("No such file ", bedFile) 
  }
  dyn.load("~/Repos/CSC/CNEr/src/io.so")
  bed = .Call("myReadBed", bedFile)
  bed = GRanges(seqnames=Rle(bed[[1]]),
                ranges=IRanges(start=bed[[2]], end=bed[[3]]),
                strand=factor("+"))
  return(bed)
}
# bedHuman = readBedToGRanges("/export/data/CNEs/hg19/filters/filter_regions.hg19.bed")
# bedZebrafish = readBedToGRanges("/export/data/CNEs/danRer7/filters/filter_regions.danRer7.bed")
# system.time(foo<-readBedToGRanges(bedFile))
#    user  system elapsed
#    2.272   0.133   2.414

# axtFiles = c("/export/downloads/ucsc/axtNet/hg19/chr2.hg19.danRer7.net.axt.gz", "/export/downloads/ucsc/axtNet/hg19/chr3.hg19.danRer7.net.axt.gz", "/export/downloads/ucsc/axtNet/hg19/chr4.hg19.danRer7.net.axt.gz")
readAxt = function(axtFiles){
  # Read axt files into R axt object
  index_noexists = !file.exists(axtFiles)
  if(any(index_noexists)){
    stop("No such file ", paste(axtFiles[index_noexists], sep=" "))
  }
  require(GenomicRanges)
  require(Biostrings)
  dyn.load("~/Repos/CSC/CNEr/src/io.so")
  myAxt = .Call("myReadAxt", axtFiles)
  axt = axt(targetRanges=GRanges(seqnames=Rle(myAxt[[1]]),
                                 ranges=IRanges(start=myAxt[[2]],
                                                end=myAxt[[3]]),
                                 strand=Rle(myAxt[[4]])),
            targetSeqs=DNAStringSet(myAxt[[5]]),
            queryRanges=GRanges(seqnames=Rle(myAxt[[6]]),
                                ranges=IRanges(start=myAxt[[7]],
                                               end=myAxt[[8]]),
                                strand=Rle(myAxt[[9]])),
            querySeqs=DNAStringSet(myAxt[[10]]),
            score=myAxt[[11]],
            symCount=myAxt[[12]]
            )
  return(axt)
}
