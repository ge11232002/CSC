
### ------------------------------------------------------------------------
### The "PWM" generic and methods. This is a bit different from the implementation of Biostrings.
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
### Currently we make it as a normal function. Is it necessary to make it a setMethod? Yes. It's necessary to make it a setMethod.
setMethod("searchSeq", "PWMatrix",
# scans a nucleotide sequence with the pattern represented by the PWM.
          function(x, subject, seqname="Unknown", strand="*", min.score="80%"){
            ans_views = matchPWM(unitScale(Matrix(x)), subject, min.score=min.score)
            #score = rep(0, length(ans_views)) # fix the score issue....
            score = PWMscoreStartingAt(unitScale(Matrix(x)), subject(ans_views),
                                       start(ans_views))
            # The score here from PWMscoreStartingAt is the unitscaled score. Let's make it into original one, synced with TFBS module. This is validated!
            score = score * (maxScore(Matrix(x)) - minScore(Matrix(x))) + minScore(Matrix(x))
            stopifnot(strand %in% c("+", "-", "*")) # need to ask Boris strand.
            if(length(strand) == 1)
              strand = rep(strand, length(ans_views))
            stopifnot(length(strand) == length(ans_views))
            ans_site = newSite(views=ans_views, seqname=seqname,
                               score=score, strand=strand, 
                               sitesource="TFBS", primary="TF binding site",
                               pattern=x
                               )
          }
          )

setMethod("searchSeq", "PWMatrixList",
# scans a nucleotide sequence with all patterns represented stored in $matrixset;
          function(x, subject, seqname="Unknown", strand="*", min.score="80%"){
            #pwms = lapply(Matrix(x), unitScale)
            #ans = lapply(pwms, matchPWM, subject, min.score)
            ans_list = lapply(x, searchSeq, subject=subject, seqname=seqname, 
                              strand=strand, min.score=min.score)
            ans = SiteList(ans_list)
            return(ans)
          }
          )

### ----------------------------------------------------------------------
### searchAln: Scans a pairwise alignment of nucleotide sequences with the pattern represented by the PWM: it reports only those hits that are present in equivalent positions of both sequences and exceed a specified threshold score in both, AND are found in regions of the alignment above the specified
### Let's make it a setMethod function for taking different subject (alignment).

setMethod("searchAln", "PWMatrix",
          function(x, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            doSiteSearch(x, aln1, aln2, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
          }
          )

setMethod("searchAln", "PWMatrixList",
          function(x, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            #ans = lapply(x, doSiteSearch, subject, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
            ans_list = lapply(x, searchAln, aln1, aln2, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
            ans = SitePairList(ans_list)
            return(ans)
          }
          )


