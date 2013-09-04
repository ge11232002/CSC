


calculate_conservation = function(aln1, aln2, windowSize, which="1"){
  ## This function is used to calculate the conservation profiles for a pairwise alignment.
  # x: a DNAStringSet with length 2 holds the alignment
  # windowSize: the smooth window size
  # which: which seq in the alignment is computed.
  if(windowSize %% 2 == 0){
    warning("windows size is not even, turned into odd by -1")
    windowSize = windowSize - 1L
  }
  which = match.arg(which, c("1", "2"))
  if(nchar(aln1) != nchar(aln2))
    stop("'aln1' and 'aln2' must have the same number of characters")

  if(which == "1"){
    alignedSeq1 = aln1
    alignedSeq2 = aln2
  }else{
    alignedSeq1 = aln2
    alignedSeq2 = aln1
  }
  alignedSeq1 = strsplit(as.character(alignedSeq1), "")[[1]]
  alignedSeq2 = strsplit(as.character(alignedSeq2), "")[[1]]
  indexGap = alignedSeq1 == "-" | alignedSeq1 == "." | alignedSeq1 == "_"
  alignedSeq1 = alignedSeq1[!indexGap]
  alignedSeq2 = alignedSeq2[!indexGap]
  matches = alignedSeq1 == alignedSeq2
  conservations = runmean(matches, k=windowSize, alg="C", endrule="mean")
  return(conservations)
}

setMethod("calConservation", signature(aln1="DNAString", aln2="DNAString"),
          function(aln1, aln2, windowSize=51L, which="1"){
            calculate_conservation(as.character(aln1), as.character(aln2),
                                   windowSize=windowSize, which=which)
          }
          )
setMethod("calConservation", signature(aln1="DNAStringSet", aln2="missing"),
          function(aln1, aln2, windowSize=51L, which="1"){
            if(length(aln1) != 2)
              stop("'aln1' must be of length 2 when 'aln2' is missing")
            calculate_conservation(as.character(aln1[1]), 
                                   as.character(aln1[2]), 
                                   windowSize=windowSize, which=which)
          }
          )
setMethod("calConservation", signature(aln1="character", aln2="missing"),
          function(aln1, aln2, windowSize=51L, which="1"){
            if(length(aln1) != 2)
              stop("'aln1' must be of length 2 when 'aln2' is missing")
            calculate_conservation(aln1[1], aln1[2], windowSize=windowSize, which=which)
          }
          )
setMethod("calConservation", signature(aln1="character", aln2="character"),
          function(aln1, aln2, windowSize=51L, which="1"){
            calculate_conservation(aln1, aln2, windowSize=windowSize, which=which)
          }
          )



