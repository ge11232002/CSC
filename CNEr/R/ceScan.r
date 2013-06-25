

# tFilter = bedHuman
# qFilter = bedZebrafish
ceScan = function(axt, tFilter=NULL, qFilter=NULL, qSizes=NULL, thresholds=c("30,40", "40,50")){
  dyn.load("~/Repos/CSC/CNEr/src/ceScanBranch1.so")
  winSize = as.integer(sapply(strsplit(thresholds, ","), "[", 2))
  minScore = as.integer(sapply(strsplit(thresholds, ","), "[", 1))
  CNE = .Call("myCeScan", as.vector(seqnames(tFilter)), start(tFilter), end(tFilter),
              as.vector(seqnames(qFilter)), start(qFilter), end(qFilter),
              as.vector(qSizes[[1]]), as.vector(qSizes[[2]]), 
              as.vector(seqnames(targetRanges(axt))), start(targetRanges(axt)), end(targetRanges(axt)), as.vector(strand(targetRanges(axt))), as.vector(targetSeqs(axt)),
              as.vector(seqnames(queryRanges(axt))), start(queryRanges(axt)), end(queryRanges(axt)), as.vector(strand(queryRanges(axt))), as.vector(querySeqs(axt)),
              score(axt), symCount(axt), winSize, minScore
              )
  return(CNE)
}

ceMerge = function(cne1, cne2){

}

# axtFiles = list.files(path="/export/downloads/ucsc/axtNet/hg19", pattern=".*hg19\\.mm10*", full.names=TRUE)
# axt1 = readAxt(axtFiles)
# axtFiles = list.files(path="/export/downloads/ucsc/axtNet/mm10", pattern=".*mm10\\.hg19.*", full.names=TRUE)
# axt2 = readAxt(axtFiles)
# thresholds=c("30,40", "40,50")
# sizes1 = data.frame(seqnames=seqnames(seqinfo(TwoBitFile("/export/data/goldenpath/hg19/assembly.2bit"))), seqlengths=seqlengths(seqinfo(TwoBitFile("/export/data/goldenpath/hg19/assembly.2bit"))), stringsAsFactors=FALSE)
# sizes2 = data.frame(seqnames=seqnames(seqinfo(TwoBitFile("/export/data/goldenpath/mm10/assembly.2bit"))), seqlengths=seqlengths(seqinfo(TwoBitFile("/export/data/goldenpath/mm10/assembly.2bit"))), stringsAsFactors=FALSE)
# filter1 = readBedToGRanges("/export/data/CNEs/hg19/filters/filter_regions.hg19.bed")
# filter2 = readBedToGRanges("/export/data/CNEs/mm10/filters/filter_regions.mm10.bed")
detectCNEs = function(axt1, filter1=NULL, sizes1, axt2, filter2=NULL, sizes2, thresholds=NULL){
  CNE1 = ceScan(axt1, filter1, filter2, sizes2, thresholds)
  CNE2 = ceScan(axt2, filter2, filter1, sizes1, thresholds)
  ceMerge(CNE1, CNE2)

}
