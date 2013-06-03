

readAxt = function(axtFile){
  axt = readLines(axtFile)
  axt = axt[grep("^#", axt, invert=TRUE, fixed=TRUE)]
  axt = axt[axt != ""]
  axt = split(axt, c("header", "target", "query"))
  header = strsplit(axt$header, split=" ", fixed=TRUE)
  myAxt = new("axt", targetNames=Rle(sapply(header, "[",2)),
              targetRanges=IRanges(start=as.integer(sapply(header, "[",3)), end=as.integer(sapply(header, "[",4))),
              targetSeqs=DNAStringSet(axt$target),
              queryNames=Rle(sapply(header, "[",5)),
              queryRanges=IRanges(start=as.integer(sapply(header, "[",6)), end=as.integer(sapply(header, "[",7))),
              querySeqs=DNAStringSet(axt$query),
              strand=Rle(sapply(header, "[",8)),
              score=rle(as.integer(sapply(header, "[",9)))
              )
  rm(c(axt, header))
  gc()
}
