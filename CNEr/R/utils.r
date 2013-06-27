
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


###
### We register the old-style (a.k.a. S3) class below as a formally defined
### class (a.k.a. S4) because we are using it in some method signatures.
### Note that dispatch still works without this registration but causes
### 'R CMD INSTALL' to (gently) complain.
###

setOldClass("probetable")


###
### Some low-level (not exported) helper functions.
###

isNumericOrNAs <- function(x)
{
    is.numeric(x) || (is.atomic(x) && is.vector(x) && all(is.na(x)))
}

normargUseNames <- function(use.names)
{
    if (is.null(use.names))
        return(TRUE)
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    use.names
}

### Returns an integer vector.
pow.int <- function(x, y)
{
    if (!is.numeric(x))
        stop("'x' must be a numeric vector")
    if (!is.integer(x))
        x <- as.integer(x)
    ans <- rep.int(1L, length(x))
    for (i in seq_len(y))
        ans <- ans * x
    ans
}








