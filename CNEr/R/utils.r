
compDNAStringSet = function(DNAStringSet1, DNAStringSet2){
  tmp = cbind(strsplit(as.character(DNAStringSet1), ""), 
              strsplit(as.character(DNAStringSet2), ""))
  apply(tmp, 1, function(x){x[[1]]==x[[2]]})
}

#system.time(foo<-compDNAStringSet(targetSeqs(myAxt), querySeqs(myAxt)))
#system.time(foo1<-RleList(foo))

makeReversedFilter = function(qFilter, chromSizes){
  revFilterBed = GRanges(seqnames=seqnames(qFilter),
                         ranges=IRanges(start=chromSizes[as.vector(seqnames(qFilter))] - end(qFilter),
                                        end=chromSizes[as.vector(seqnames(qFilter))] - start(qFilter)
                                        ),
                         strand=Rle("-"))
  return(revFilterBed)
}

seqToAlignment = function(DNAStringSet){
  foo = strsplit(as.character(DNAStringSet), "")
  foo = lapply(foo, function(x){grep("-", x, invert=TRUE)})
  return(foo)
}

reverseCigar = function(cigar){
  require(GenomicRanges)
  cigar = sapply(splitCigar(cigar), function(x){
                 paste0(rev(x[[2]]), rev(rawToChar(x[[1]], multiple=TRUE)), 
                        collapse="")
                         }
  )
  return(cigar)
}

my.system = function(cmd, echo=TRUE, intern=FALSE, ...){
  if (echo){
    message(cmd)
  }
  res = system(cmd, intern=intern, ...)
  if (!intern){
    stopifnot(res == 0)
  }
  return(res)
}
