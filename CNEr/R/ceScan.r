

# tFilter = bedHuman
# qFilter = bedZebrafish
ceScan = function(axts, tFilter=NULL, qFilter=NULL, qSizes=NULL, thresholds=c("30,40", "40,50")){
  ## Here the returned tStart and qStart are in 1-based coordinates.
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

# cne1 = read.table("/mnt/biggley/home/gtan/debug/ceScan_C_Filter/11-07-2013/cne_hg19_danRer7_30_40", header=FALSE, sep="\t")
# cne2 = read.table("/mnt/biggley/home/gtan/debug/ceScan_C_Filter/11-07-2013/cne_danRer7_hg19_30_40", header=FALSE, sep="\t")
# colnames(cne1) = c("tName", "tStart", "tEnd", "qName", "qStart", "qEnd", "strand", "score", "cigar")
# colnames(cne2) = c("tName", "tStart", "tEnd", "qName", "qStart", "qEnd", "strand", "score", "cigar")
ceMerge = function(cne1, cne2){
  require(GenomicRanges)
  ## clean the +1 leater as in R, we already use the 1-based start coordinates.
  cne1T = GRanges(seqnames=cne1$tName, ranges=IRanges(start=cne1$tStart+1, end=cne1$tEnd), strand="+")
  cne1Q = GRanges(seqnames=cne1$qName, ranges=IRanges(start=cne1$qStart+1, end=cne1$qEnd), strand=cne1$strand)
  cne2T = GRanges(seqnames=cne2$qName, ranges=IRanges(start=cne2$qStart+1, end=cne2$qEnd), strand="+")
  cne2Q = GRanges(seqnames=cne2$tName, ranges=IRanges(start=cne2$tStart+1, end=cne2$tEnd), strand=cne2$strand)
  cneT = c(cne1T, cne2T)
  cneQ = c(cne1Q, cne2Q)
  cneT_overlap = findOverlaps(cneT, type="within", ignoreSelf=TRUE, ignoreRedundant=TRUE)
  #cneT_overlap1 = findOverlaps(cneT, type="equal", ignoreSelf=TRUE, ignoreRedundant=TRUE)
  #cneT_overlap2 = findOverlaps(cneT, type="any", ignoreSelf=TRUE, ignoreRedundant=TRUE)
  cneQ_overlap = findOverlaps(cneQ, type="within", ignoreSelf=TRUE, ignoreRedundant=TRUE)
  #cneQ_overlap1 = findOverlaps(cneQ, type="equal", ignoreSelf=TRUE, ignoreRedundant=TRUE)
  #cneQ_overlap2 = findOverlaps(cneQ, type="any", ignoreSelf=TRUE, ignoreRedundant=TRUE)
  redundance = intersect(cneT_overlap, cneQ_overlap)
  res = rbind(cne1, cne2)[-queryHits(redundance), ] 
  return(res)
  #stopifnot(length(setdiff (cneT[queryHits(redundance)],  cneT[subjectHits(redundance)])) != 0)
  #stopifnot(length(setdiff (cneQ[queryHits(redundance)],  cneQ[subjectHits(redundance)])) != 0)
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
