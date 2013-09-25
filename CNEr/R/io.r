# tFilterFile = "/export/data/CNEs/danRer7/filters/filter_regions.danRer7.bed"
# qFilterFile = "/export/data/CNEs/hg19/filters/filter_regions.hg19.bed"
readBed = function(bedFile){
## This GRanges is in 1-based.
  require(rtracklayer)
  bed = import(bedFile, asRangedData = FALSE)
  strand(bed) = "+"
  bed = reduce(bed)
}

#############################C version######################################
readBedToGRanges = function(bedFile=NULL){
## This GRanges have the different coordinates system with the original bed file. i.e. with 1-based start end coordinates.
  if(is.null(bedFile)){
    return(NULL)
  }
  if(!file.exists(bedFile)){
    stop("No such file ", bedFile) 
  }
  #dyn.load("~/Repos/CSC/CNEr/src/CNEr.so")
  bed = .Call2("myReadBed", bedFile, PACKAGE="CNEr")
  bed = GRanges(seqnames=Rle(bed[[1]]),
                ranges=IRanges(start=bed[[2]], end=bed[[3]]),
                strand=factor("+"))
  return(bed)
}

readAxt = function(axtFiles){
  # Read axt files into R axt object.
  # The coordinates are 1-based for start and end.
  index_noexists = !file.exists(axtFiles)
  if(any(index_noexists)){
    stop("No such file ", paste(axtFiles[index_noexists], sep=" "))
  }
  require(GenomicRanges)
  require(Biostrings)
  #dyn.load("~/Repos/CSC/CNEr/src/CNEr.so")
  myAxt = .Call2("readAxt", axtFiles, PACKAGE="CNEr")
  #myAxtInfo = .Call("axt_info", axtFiles)
  axts = axt(targetRanges=GRanges(seqnames=Rle(myAxt[[1]]),
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
  return(axts)
}

axtInfo = function(axtFiles){
  #dyn.load("~/Repos/CSC/CNEr/src/CNEr.so")
  ans = .Call2("axt_info", axtFiles, PACKAGE="CNEr")
  return(ans)
}

