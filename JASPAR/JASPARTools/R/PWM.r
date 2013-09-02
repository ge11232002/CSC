
### ------------------------------------------------------------------------
### The "PWM" generic and methods. This is a bit different from the implementation of Biostrings.


setGeneric("toPWM", signature="x",
           function(x, pseudocounts=NULL,
                    bg_probabilities=c(A=0.25, C=0.25, G=0.25, T=0.25))
             standardGeneric("toPWM")
           )

setMethod("toPWM", "character",
          function(x, pseudocounts=NULL, 
                   bg_probabilities=c(A=0.25, C=0.25, G=0.25, T=0.25)){
            if(is.null(pseudocounts))
              pseudocounts = 0.8
            dnaset = DNAStringSet(x)
            toPWM(dnaset, pseudocounts=pseudocounts,
                  bg_probabilities=bg_probabilities)
          }
          )
setMethod("toPWM", "DNAStringSet",
          function(x, pseudocounts=NULL,
                   bg_probabilities=c(A=0.25, C=0.25, G=0.25, T=0.25)){
            if(!isConstant(width(x)))
              stop("'x' must be rectangular (i.e. have a constant width)")
            if(is.null(pseudocounts))
              pseudocounts = 0.8
            pfm = consensusMatrix(x)
            toPWM(pfm, pseudocounts=pseudocounts,
                  bg_probabilities=bg_probabilities)
          }
          )
setMethod("toPWM", "PFMatrix",
          function(x, pseudocounts=NULL){
            if(is.null(pseudocounts))
              pseudocounts = 0.8
            pwmMatrix = toPWM(Matrix(x), pseudocounts=pseudocounts,
                              bg_probabilities=bg(x))
            pwm = PWMatrix(ID=ID(x), name=name(x), matrixClass=matrixClass(x),
                           strand=strand(x), bg=bg(x), matrix=pwmMatrix,
                           pseudocounts=pseudocounts)
          }
          )

### Assumes 'x' is a Position *Frequency* Matrix (PFM) and computes the
### corresponding Position *Weight* Matrix (PWM).
setMethod("toPWM", "matrix",
    ## This is validated by the TFBS perl module version.
          function(x, pseudocounts=NULL,
                   bg_probabilities=c(A=0.25, C=0.25, G=0.25, T=0.25)){
            x = Biostrings:::.normargPfm(x)
            bg_probabilities = Biostrings:::.normargPriorParams(bg_probabilities)
            nseq = colSums(x)
            priorN = sum(bg_probabilities)
            if(is.null(pseudocounts))
              pseudocounts = 0.8
            if(length(pseudocounts) == 1)
              p = sweep(x + bg_probabilities*pseudocounts, MARGIN=2, nseq + pseudocounts, "/")
              #p = (x + bg_probabilities*pseudocounts) / (nseq + pseudocounts)
            else
              #p = (x + bg_probabilities %*% t(pseudocounts)) / (nseq + pseudocounts)
              p = sweep(x + bg_probabilities %*% t(pseudocounts), MARGIN=2, nseq + pseudocounts, "/")
            prior.probs = bg_probabilities / priorN
            #ans = log2(p / prior.probs)
            #Here ans's colSums is 1s. Need to be adapted for seq logo maybe later.
            ans = log2(sweep(p, MARGIN=1, prior.probs, "/"))
            return(ans)
          }
          )
### ---------------------------------------------------------------------
### searchSeq: scans a nucleotide sequence with the pattern represented by the PWM
### Currently we make it as a normal function. Is it necessary to make it a setMethod?
searchSeq = function(pwm, subject, min.score="80%", ...){
  pwmMatrix = unitScale(Matrix(pwm))
  matchPWM(pwm, subject, min.score=min.score)
}

### ----------------------------------------------------------------------
### searchAln: Scans a pairwise alignment of nucleotide sequences with the pattern represented by the PWM: it reports only those hits that are present in equivalent positions of both sequences and exceed a specified threshold score in both, AND are found in regions of the alignment above the specified
### Let's make it a setMethod function for taking different subject (alignment).
setGeneric("searchAln", signature="x",
           function(pwm, x, min.score="80%", windowSize=50L, cutoff="70%",
                    conservation=NULL, ...)
             standardGeneric("searchAln")
           )

calculate_conservation = function(x, windowSize, which=c("1", "2")){
  ## This function is used to calculate the conservation profiles for a pairwise alignment.
  # x: a DNAStringSet with length 2 holds the alignment
  # windowSize: the smooth window size
  # which: which seq in the alignment is computed.
  if(windowSize %% 2 == 0){
    warning("windows size is not even, turned into odd by -1")
    windowSize = windowSize - 1L
  }
  which = match.arg(which)
  stopifnot(length(x) == 2)
  if(which == "1"){
    alignedSeq1 = x[1]
    alignedSeq2 = x[2]
  }else{
    alignedSeq1 = x[2]
    alignedSeq2 = x[1]
  }
  stopifnot(nchar(alignedSeq1) == nchar(alignedSeq2))
  alignedSeq1 = strsplit(as.character(alignedSeq1), "")[[1]]
  alignedSeq2 = strsplit(as.character(alignedSeq2), "")[[1]]
  indexGap = alignedSeq1 == "-" | alignedSeq1 == "."
  alignedSeq1 = alignedSeq1[!indexGap]
  alignedSeq2 = alignedSeq2[!indexGap]
  matches = alignedSeq1 == alignedSeq2
  conservations = runmean(matches, k=windowSize, alg="C", endrule="mean")
}

do_sitesearch = function(pwm, x, min.score, windowSize, cutoff, conservation){
}

setMethod("searchAln", "DNAStringSet",
          function(pwm, x, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            windowSize = as.integer(windowSize)
            if(cutoff > 1 || cutoff < 0)
              stop("cutoff must be from 0 to 1.")
            do_sitesearch(pwm, x, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
          }
          )

}
