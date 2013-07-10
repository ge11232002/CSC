

# tFilter = bedHuman
# qFilter = bedZebrafish
ceScan = function(axts, tFilter=NULL, qFilter=NULL, qSizes=NULL, thresholds=c("30,40", "40,50")){
  dyn.load("~/Repos/CSC/CNEr/src/CNEr.so")
  if(!is.null(qFilter))
    if(is.null(qSizes) || !is(qSizes, "Seqinfo"))
      stop("qSizes must exist and be a Seqinfo object when qFilter exists")
  
  winSize = as.integer(sapply(strsplit(thresholds, ","), "[", 2))
  minScore = as.integer(sapply(strsplit(thresholds, ","), "[", 1))
  resFiles = tempfile(pattern = paste(minScore, winSize, "ceScan", sep="-"), tmpdir = tempdir(), fileext = "")
  .Call("myCeScan", as.vector(seqnames(tFilter)), start(tFilter), end(tFilter),
              as.vector(seqnames(qFilter)), start(qFilter), end(qFilter),
              as.vector(seqnames(qSizes)), as.vector(seqlengths(qSizes)), 
              as.vector(seqnames(targetRanges(axts))), start(targetRanges(axts)), end(targetRanges(axts)), as.vector(strand(targetRanges(axts))), as.vector(targetSeqs(axts)),
              as.vector(seqnames(queryRanges(axts))), start(queryRanges(axts)), end(queryRanges(axts)), as.vector(strand(queryRanges(axts))), as.vector(querySeqs(axts)),
              score(axts), symCount(axts), winSize, minScore, as.vector(resFiles)
              )
  CNE = lapply(resFiles, 
               function(x){res=read.table(x, header=FALSE, sep="\t")
               colnames(res)=c("tName", "tStart", "tEnd", "qName", "qStart", "qEnd", "strand", "score", "cigar")
               return(res)})
  names(CNE) =  paste(minScore, winSize, sep="_")
  unlink(resFiles)
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
