
### ---------------------------------------------------------------------
### Site object: a nucleotide sequence feature object representing (possibly putative) transcription factor binding site. Different from TFBS perl module, here one Site object contains multiple sites.
###

setClass("Site",
         slots=c(views="XStringViews",
                 score="numeric",    # vector 
                 strand="character",  ## make it Rle() later.
                 seqname="character", # length 1
                 sitesource="character", # length 1
                 primary="character",  # length 1
                 pattern="PWMatrix"   # length 1
                 )
         )

### -------------------------------------------------------------------
### The accessor-like method
###
setGeneric("views", signature="x", function(x) standardGeneric("views"))
setMethod("views", "Site", function(x) x@views)

setMethod("score", "Site", function(x) x@score)
setMethod("strand", "Site", function(x) x@strand)

setGeneric("seqname", signature="x", function(x) standardGeneric("seqname"))
setMethod("seqname", "Site", function(x) x@seqname)

setGeneric("sitesource", signature="x", function(x) standardGeneric("sitesource"))
setMethod("sitesource", "Site", function(x) x@sitesource)

setGeneric("primary", signature="x", function(x) standardGeneric("primary"))
setMethod("primary", "Site", function(x) x@primary)

setGeneric("pattern", signature="x", function(x) standardGeneric("pattern"))
setMethod("pattern", "Site", function(x) x@pattern)

setMethod("length", "Site", function(x) length(views(x)))
### -------------------------------------------------------------------
### The constructor
###
newSite = function(views, score, strand="*",
                seqname="Unknown",
                sitesource="TFBS", primary="TF binding site",
                pattern){
  new("Site", views=views, seqname=seqname, score=score, strand=strand, 
      sitesource=sitesource, primary=primary, pattern=pattern)
}

### -----------------------------------------------------------------
### Methods
###
setGeneric("writeGFF3", signature="x", function(x) standardGeneric("writeGFF3"))
setMethod("writeGFF3", "Site",
          function(x){
            gff = list(seqname=seqname(x),
                       source=sitesource(x),
                       feature=sitesource(x),
                       start=start(views(x)),
                       end=end(views(x)),
                       score=score(x),
                       strand=strand(x),
                       frame=".",
                       attributes=paste(
                                        paste("TF", name(pattern(x)), sep="="),
                                        paste("class", matrixClass(pattern(x)), sep="="),
                                        paste("sequence", as.character(views(x)), sep="="),
                                        sep=";")
                       )
            gff = as.data.frame(gff)
            return(gff)
          }
          )
setGeneric("writeGFF2", signature="x", function(x) standardGeneric("writeGFF2"))
setMethod("writeGFF2", "Site",
          function(x){
            gff = list(seqname=seqname(x),
                       source=sitesource(x),
                       feature=sitesource(x),
                       start=start(views(x)),
                       end=end(views(x)),
                       score=score(x),
                       strand=strand(x),
                       frame=".",
                       attributes=paste(
                                        paste("TF", paste0("\"", name(pattern(x)), "\""), sep=" "),
                                        paste("class", paste0("\"", matrixClass(pattern(x)), "\""), sep=" "),
                                        paste("sequence", paste0("\"", as.character(views(x)), "\""), sep=" "),
                                        sep="; ")
                       )
            gff = as.data.frame(gff)
            return(gff)
          }
          )

          
### -----------------------------------------------------------------
### The "show" method
### Perhaps it is not a bad idea to show them in gff format.
setMethod("show", "Site",
          function(object){
            gff = writeGFF3(object)
            cat("An object of class", class(object), "with", 
                length(object), "site", 
                ifelse(length(object)==1, "sequence", "sequences"))
            cat("\n")
            print(gff)
          }
          )


