
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
setMethod("views", "Site", function(x) x@views)

setMethod("score", "Site", function(x) x@score)
setMethod("strand", "Site", function(x) x@strand)

setMethod("seqname", "Site", function(x) x@seqname)

setMethod("sitesource", "Site", function(x) x@sitesource)

setMethod("primary", "Site", function(x) x@primary)

setMethod("pattern", "Site", function(x) x@pattern)

setMethod("length", "Site", function(x) length(views(x)))
### -------------------------------------------------------------------
### The constructor
###
Site = function(views, score, strand="*",
                seqname="Unknown",
                sitesource="TFBS", primary="TF binding site",
                pattern){
  new("Site", views=views, seqname=seqname, score=score, strand=strand, 
      sitesource=sitesource, primary=primary, pattern=pattern)
}

### ------------------------------------------------------------------
### The getters
###
setMethod("[", "Site", 
          function(x, i){
            if(missing(i))
              return(x)
            ans_views = views(x)[i]
            ans_score = score(x)[i]
            ans_strand = strand(x)[i]
            clone(x, views=ans_views, score=ans_score, strand=ans_strand)
          }
          )

### -----------------------------------------------------------------
### Methods
###
setMethod("writeGFF3", "Site",
          function(x){
            if(length(x) == 0)
              return(data.frame())
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
setMethod("writeGFF2", "Site",
          function(x){
            if(length(x) == 0)
              return(data.frame())
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

setMethod("relScore", "Site",
          function(x){
          # Luckliy, the maxScore, minScore implementation is same with TFBS perl module. Validated!
            ans = (score(x) - minScore(Matrix(pattern(x)))) / (maxScore(Matrix(pattern(x))) - minScore(Matrix(pattern(x))))
            return(ans)
          }
          )
          
### -----------------------------------------------------------------
### The "show" method
### Perhaps it is not a bad idea to show them in gff format.
setMethod("show", "Site",
          function(object){
            cat("An object of class", class(object), "with", 
                length(object), "site", 
                ifelse(length(object)==1, "sequence", "sequences"))
            cat("\n")
            if(length(object) > 10000)
              object = object[1:10000]
            gff = writeGFF3(object)
            print(gff)
          }
          )


