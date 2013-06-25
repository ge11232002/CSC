library(Biostrings)
library(GenomicRanges)
library(rtracklayer)

setClass(Class="axt",
         representation(targetRanges="GRanges",
                        targetSeqs="DNAStringSet",
                        queryRanges="GRanges",
                        querySeqs="DNAStringSet",
                        score="integer",
                        symCount="integer"
                        )
         )



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Slot getters and setters.
###
setGeneric("targetRanges", function(x) standardGeneric("targetRanges"))
setMethod("targetRanges", "axt", function(x) x@targetRanges)
setGeneric("targetSeqs", function(x) standardGeneric("targetSeqs"))
setMethod("targetSeqs", "axt", function(x) x@targetSeqs)
setGeneric("queryRanges", function(x) standardGeneric("queryRanges"))
setMethod("queryRanges", "axt", function(x) x@queryRanges)
setGeneric("querySeqs", function(x) standardGeneric("querySeqs"))
setMethod("querySeqs", "axt", function(x) x@querySeqs)
setMethod("score", "axt", function(x) x@score)
setGeneric("symCount", function(x) standardGeneric("symCount"))
setMethod("symCount", "axt", function(x) x@symCount)
setMethod("length", "axt", function(x) length(targetRanges(x)))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###
axt = function(targetRanges=GRanges(), targetSeqs=DNAStringSet(),
               queryRanges=GRanges(), querySeqs=DNAStringSet(),
               score=integer(), symCount=integer()){
  new("axt", targetRanges=targetRanges, targetSeqs=targetSeqs,
      queryRanges=queryRanges, querySeqs=querySeqs,
      score=score, symCount=symCount)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Updating and cloning.
###
### An object is either 'update'd in place (usually with a replacement
### method) or 'clone'd (copied), with specified slots/fields overridden.

### For an object with a pure S4 slot representation, these both map to
### initialize. Reference classes will want to override 'update'. Other
### external representations need further customization.

setMethod("update", "axt",
          function(object, ..., check=TRUE){
            initialize(object, ...)
          }
)

setGeneric("clone", function(x, ...) standardGeneric("clone"))  # not exported
setMethod("clone", "ANY",  # not exported
    function(x, ...)
    {
        if (nargs() > 1L)
            initialize(x, ...)
        else
            x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting and combining.
###

setMethod("[", "axt",
          function(x, i, ..., drop){
            if(length(list(...)) > 0L)
              stop("invalid subsetting")
            if(missing(i))
              return(x)
            i = IRanges:::normalizeSingleBracketSubscript(i, x)
            ans_targetRanges = targetRanges(x)[i]
            ans_targetSeqs = targetSeqs(x)[i]
            ans_queryRanges = queryRanges(x)[i]
            ans_querySeqs = querySeqs(x)[i]
            ans_score = score(x)[i]
            ans_symCount = symCount(x)[i]
            clone(x, targetRanges=ans_targetRanges, targetSeqs=ans_targetSeqs,
                  queryRanges=ans_queryRanges, querySeqs=ans_querySeqs,
                  score=ans_score, symCount=ans_symCount)
          }
          )

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method.
###

#showAxt = function(x, margin=""){
#  lx = length(x)
#  cat(class(x), " with ",
#      lx, " ", ifelse(lx == 1L, "alignment pair", "alignment pairs"),
#                  ":\n", sep="")
#    out = makePrettyMatrixForCompactPrintingAxt(x, .makeNakedMatFromAxt)
#    if(nrow(out) != 0L)
#          rownames(out) = paste0(margin, rownames(out))
#      print(out, quote=FALSE, right=TRUE)
#}
#
#setMethod("show", "axt",
#          function(object)
#            showAxt(object, margin="  ")
#)
