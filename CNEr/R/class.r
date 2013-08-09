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
            #i = IRanges:::normalizeSingleBracketSubscript(i, x)
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
setGeneric("subAxt", function(x, searchGRanges, select=c("target", "query"), type=c("any", "within")) standardGeneric("subAxt"))
setMethod("subAxt", "axt",
## This is to fetch the axts within the specific chrs, starts, ends based on target sequences.
          function(x, searchGRanges, select=c("target", "query"),
                   type=c("any", "widthin")){
            type = match.arg(type)
            select = match.arg(select)
            if(length(searchGRanges) == 0)
              return(x)
            if(select == "target"){
              index = which(!is.na(findOverlaps(targetRanges(x), searchGRanges, type=type, select="first")))
            }else{
              index = which(!is.na(findOverlaps(queryRanges(x), searchGRanges, type=type, select="first")))
            }
            return(x[index])
          }
          )

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method.
###
### 'x' must be an XString or MaskedXString object.
toSeqSnippet <- function(x, width)
{
    if (width < 7L)
        width <- 7L
    seqlen <- length(x)
    if (seqlen <= width) {
        as.character(x)
    } else {
        w1 <- (width - 2) %/% 2
        w2 <- (width - 3) %/% 2
        paste(as.character(subseq(x, start=1, width=w1)),
              "...",
              as.character(subseq(x, end=seqlen, width=w2)),
              sep="")
    }
}


.axt.show_frame_line = function(x, i, iW, tNameW, tStartW, tEndW, qNameW, qStartW, qEndW, scoreW){
  cat(format(i, width=iW, justify="right"), " ",
      format(as.vector(seqnames(targetRanges(x)[i])), width=tNameW, justify="right"), " ",
      format(start(targetRanges(x)[i]), width=tStartW, justify="right"), " ",
      format(end(targetRanges(x)[i]), width=tEndW, justify="right"), " ",
      format(as.vector(seqnames(queryRanges(x)[i])), width=qNameW, justify="right"), " ",
      format(start(queryRanges(x)[i]), width=qStartW, justify="right"), " ",
      format(end(queryRanges(x)[i]), width=qEndW, justify="right"), " ",
      format(as.vector(strand(queryRanges(x))[i]), width=1, justify="right"), " ",
      format(score(x)[i], width=scoreW, justify="right"), " ",
      sep=""
      )
  cat("\n")
  snippetWidth = getOption("width")
  seq_snippet = toSeqSnippet(targetSeqs(x)[[i]], snippetWidth)
  cat(seq_snippet)
  cat("\n")
  seq_snippet = toSeqSnippet(querySeqs(x)[[i]], snippetWidth)
  cat(seq_snippet)
  cat("\n")
}

showAxt = function(x, margin="", half_nrow=5L){
  lx = length(x)
  if(is.null((head_nrow = getOption("showHeadLines"))))
    head_nrow = half_nrow
  if(is.null((tail_nrow = getOption("showTailLines"))))
    tail_nrow = half_nrow
  iW = nchar(as.character(lx))
  if(lx < (2*half_nrow+1L) | (lx < (head_nrow+tail_nrow+1L))) {
    tNameW = max(nchar(as.vector(seqnames(targetRanges(x)))))
    tStartW = max(nchar(as.character(start(targetRanges(x)))))
    tEndW = max(nchar(as.character(end(targetRanges(x)))))
    qNameW = max(nchar(as.vector(seqnames(queryRanges(x)))))
    qStartW = max(nchar(as.character(start(queryRanges(x)))))
    qEndW = max(nchar(as.character(end(queryRanges(x)))))
    scoreW = max(nchar(as.character(score(x))))
    for(i in seq_len(lx))
      .axt.show_frame_line(x, i, iW, tNameW, tStartW, tEndW, qNameW, qStartW, qEndW, scoreW)
  }else{
    tNameW = max(nchar(as.vector(seqnames(targetRanges(x)[c(1:head_nrow, (lx-tail_nrow+1L):lx)]))))
    tStartW = max(nchar(as.character(start(targetRanges(x)[c(1:head_nrow, (lx-tail_nrow+1L):lx)]))))
    tEndW = max(nchar(as.character(end(targetRanges(x)[c(1:head_nrow, (lx-tail_nrow+1L):lx)]))))
    qNameW = max(nchar(as.vector(seqnames(queryRanges(x)[c(1:head_nrow, (lx-tail_nrow+1L):lx)]))))
    qStartW = max(nchar(as.character(start(queryRanges(x)[c(1:head_nrow, (lx-tail_nrow+1L):lx)]))))
    qEndW = max(nchar(as.character(end(queryRanges(x)[c(1:head_nrow, (lx-tail_nrow+1L):lx)]))))
    scoreW = max(nchar(as.character(score(x)[c(1:head_nrow, (lx-tail_nrow+1L):lx)])))
    if(head_nrow > 0)
      for(i in 1:head_nrow)
        .axt.show_frame_line(x, i, iW, tNameW, tStartW, tEndW, qNameW, qStartW, qEndW, scoreW)
    cat(format("...", width=iW, justify="right"),
        format("...", width=tNameW, justify="right"),
        format("...", width=tStartW, justify="right"),
        format("...", width=tEndW, justify="right"),
        format("...", width=qNameW, justify="right"),
        format("...", width=qStartW, justify="right"),
        format("...", width=qEndW, justify="right"),
        format("...", width=scoreW, justify="right")
        )
    cat("\n")
    if(tail_nrow > 0)
      for(i in (lx-tail_nrow+1L):lx)
        .axt.show_frame_line(x, i, iW, tNameW, tStartW, tEndW, qNameW, qStartW, qEndW, scoreW)
  }
}
    #out = makePrettyMatrixForCompactPrintingAxt(x, .makeNakedMatFromAxt)
    #if(nrow(out) != 0L)
    #      rownames(out) = paste0(margin, rownames(out))
    #  print(out, quote=FALSE, right=TRUE)

setMethod("show", "axt",
          function(object){
            lx = length(object)
            cat(" A ", class(object), " with ", length(object), " ", 
                ifelse(lx == 1L, "alignment pair", "alignment pairs"), ":\n", sep="")
            if(length(object) != 0){
              showAxt(object, margin="  ")
            }
          }
)
