
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
            dnaset = DNAStringSet(x)
            toPWM(dnaset, bg_probabilities=bg_probabilities)
          }
          )
setMethod("toPWM", "DNAStringSet",
          function(x, pseudocounts=NULL,
                   bg_probabilities=c(A=0.25, C=0.25, G=0.25, T=0.25)){
            if(!isConstant(width(x)))
              stop("'x' must be rectangular (i.e. have a constant width)")
            pfm = consensusMatrix(x)
            toPWM(pfm, bg_probabilities=bg_probabilities)
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
              p = (x + bg_probabilities*pseudocounts) / (nseq + pseudocounts)
            else
              p = (x + bg_probabilities %*% t(pseudocounts)) / (nseq + pseudocounts)
            prior.probs = bg_probabilities / priorN
            #ans = log2(p / prior.probs)
            #Here ans's colSums is 1s. Need to be adapted for seq logo maybe later.
            ans = log2(sweep(p, MARGIN=1, prior.probs, "/"))
            return(ans)
          }
          )



