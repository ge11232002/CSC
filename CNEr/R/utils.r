
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

### rever the cigar string. i.e. 20M15I10D will be reversed to 10D15I20M.
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

binFromCoordRange = function(starts, ends){
  bins = .Call("bin_from_coord_range", as.integer(starts), as.integer(ends))
  return(bins)
}


get_cne_ranges_in_region = function(CNE, whichAssembly=c(1,2), chr, CNEstart, CNEend, min_length){
 ## This CNE data.frame does not have the bin column yet. I am not sure whether it is necessary to add this column in R since it's quiet fast to select the cnes which meet the criteria (~0.005 second).
  if(whichAssembly == 1)
    res = subset(CNE, chr1==chr & start1>=CNEstart & end1<=CNEend & end1-start1+1>=min_length & end2-start2+1>=min_length, select=c("start1", "end1"))
  else if(whichAssembly == 2)
    res = subset(CNE, chr2==chr & start2>=CNEstart & end2<=CNEend & end1-start1+1>=min_length & end2-start2+1>=min_length, select=c("start2", "end2"))
  else
    stop("whichAssembly should be 1 or 2")
  # Here we return a IRanges object to store the start and end
  res = IRanges(start=res[ ,1], end=res[ ,2])
  return(res)
}


get_cne_ranges_in_region_partitioned_by_other_chr = function(CNE, whichAssembly=c(1,2), chr, CNEstart, CNEend, min_length){
  if(whichAssembly == 1)
    res = subset(CNE, chr1==chr & start1>=CNEstart & end1<=CNEend & end1-start1+1>=min_length & end2-start2+1>=min_length, select=c("chr2", "start1", "end1"))
  else if(whichAssembly == 1)
    res = subset(CNE, chr2==chr & start2>=CNEstart & end2<=CNEend & end1-start1+1>=min_length & end2-start2+1>=min_length, select=c("chr1", "start2", "end2"))
  else
    stop("whichAssembly should be 1 or 2")
  # Here we return a GRanges object.
  res = GRanges(seqnames=res[ ,1], ranges=IRanges(start=res[ ,2], end=res[ ,3]))
  return(res)
}

saveCNEToSQLite = function(CNE, dbName, tableName, overwrite=FALSE){
  require(RSQLite)
  CNE$bin1 = binFromCoordRange(CNE$start1, CNE$end1)
  CNE$bin2 = binFromCoordRange(CNE$start2, CNE$end2)
  # reorder it
  CNE = CNE[ ,c("bin1", "chr1", "start1", "end1", "bin2", "chr2", "start2", "end2", "strand", "similarity", "cigar")]
  drv = dbDriver("SQLite")
  dbName = dbName
  con = dbConnect(drv, dbname=dbName)
  dbWriteTable(con, tableName, CNE, row.names=FALSE, overwrite=overwrite)
  dbDisconnect(con)
}

readCNEFromSQLite = function(dbName, tableName, whichAssembly=c(1,2))

