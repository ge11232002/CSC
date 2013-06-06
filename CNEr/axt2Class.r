

setClass(Class="axt2",
         representation(targetRanges="GRanges",
                        targetSeqs="DNAStringSet",
                        queryRanges="GRanges",
                        querySeqs="DNAStringSet",
                        score="integer"
                        )
         )



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Slot getters and setters.
###
setMethod("targetRanges", "axt2", function(x) x@targetRanges)
setMethod("targetSeqs", "axt2", function(x) x@targetSeqs)
setMethod("queryRanges", "axt2", function(x) x@queryRanges)
setMethod("querySeqs", "axt2", function(x) x@querySeqs)
setMethod("target", "axt2", function(x) x@target)
setMethod("score", "axt2", function(x) x@score)
setMethod("length", "axt2", function(x) length(targetRanges(x)))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###
axt2 = function(targetRanges=GRanges(), targetSeqs=DNAStringSet(),
                queryRanges=GRanges(), querySeqs=DNAStringSet(),
                score=integer()){
  new("axt2", targetRanges=targetRanges, targetSeqs=targetSeqs,
      queryRanges=queryRanges, querySeqs=querySeqs,
      score=score)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Updating and cloning.
###
### An object is either 'update'd in place (usually with a replacement
### method) or 'clone'd (copied), with specified slots/fields overridden.

### For an object with a pure S4 slot representation, these both map to
### initialize. Reference classes will want to override 'update'. Other
### external representations need further customization.

setMethod("update", "axt2",
          function(object, ..., check=TRUE){
            initialize(object, ...)
          }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting and combining.
###

setMethod("[", "axt2",
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
            clone(x, targetRanges=ans_targetRanges, targetSeqs=ans_targetSeqs,
                  queryRanges=ans_queryRanges, querySeqs=ans_querySeqs,
                  score=ans_score)
          }
          )

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method.
###

showAxt = function(x, margin=""){
  lx = length(x)
  cat(class(x), " with ",
      lx, " ", ifelse(lx == 1L, "alignment pair", "alignment pairs"),
                  ":\n", sep="")
    out = makePrettyMatrixForCompactPrintingAxt(x, .makeNakedMatFromAxt)
    if(nrow(out) != 0L)
          rownames(out) = paste0(margin, rownames(out))
      print(out, quote=FALSE, right=TRUE)
}

setMethod("show", "axt2",
          function(object)
            showAxt2(object, margin="  ")
)
