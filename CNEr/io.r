

readAxt = function(axtFile){
  axtContent = readLines(axtFile)
  axtContent = axtContent[grep("^#", axtContent, invert=TRUE, fixed=TRUE)]
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
  rm(axtContent)
  rm(header)
  gc()
  return(myAxt)
}



