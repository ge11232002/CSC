

### ------------------------------------------------------------------------
### The "ICM" generic and methods
setGeneric("ICM", signature="x",
           function(x, pseudocounts=NULL,
                    bg_probabilities=c(A=0.25, C=0.25, G=0.25, T=0.25))
             standardGeneric("ICM")
           )

setMethod("ICM", "character",
          function(x, pseudocounts=NULL,
                   bg_probabilities=c(A=0.25, C=0.25, G=0.25, T=0.25)){
            dnaset = DNAStringSet(x)
            ICM(dnaset, bg_probabilities=bg_probabilities)
          }
          )
setMethod("ICM", "DNAStringSet",
          function(x, pseudocounts=NULL,
                   bg_probabilities=c(A=0.25, C=0.25, G=0.25, T=0.25)){
            if(!isConstant(width(x)))
              stop("'x' must be rectangular (i.e. have a constant width)")
            pfm = consensusMatrix(x)
            ICM(pfm, bg_probabilities=bg_probabilities)
          }
          )

### Assumes 'x' is a Position *Frequency* Matrix (PFM) and computes the
### corresponding Information *Content* Matrix (ICM).
setMethod("ICM", "matrix",
          function(x, pseudocounts=0.8, ## This is the recommended value from http://nar.oxfordjournals.org/content/37/3/939.long.
                   bg_probabilities=c(A=0.25, C=0.25, G=0.25, T=0.25)){
            x = Biostrings:::.normargPfm(x)
            ## From here 'x' is guaranteed to have at least 1 column and to have
            ## all its columns sum to the same value.
            bg_probabilities = Biostrings:::.normargPriorParams(bg_probabilities)
            nseq = sum(x[ ,1L])
            p = (x + bg_probabilities*pseudocounts) / (nseq + pseudocounts)
            D = log2(nrow(x)) + colSums(p * log2(p))
            ICM = t(t(p) * D)
            return(ICM)
          }
          )


