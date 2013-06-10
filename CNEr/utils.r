
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

