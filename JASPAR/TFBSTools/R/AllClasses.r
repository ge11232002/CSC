### -----------------------------------------------------------------
### The position frequency matrix (PWM) class
###
setClass("PFMatrix",
         slots=c(ID="character",
                 name="character",
                 matrixClass="character",
                 strand="character",
                 bg="numeric",
                 tags="list",
                 matrix="matrix")
         )

setClass("PWMatrix", contains="PFMatrix",
         slots=c(pseudocounts="numeric")
         )
setClass("ICMatrix", contains="PWMatrix",
         slots=c(
                 schneider="logical"
                 )
         )
setClassUnion("XMatrix", c("PFMatrix", "ICMatrix", "PWMatrix"))


### ----------------------------------------------------------------------
### The PFMatrix constructor
###
PFMatrix = function(ID="Unknown", name="Unknown", matrixClass="Unknown",
                    strand="*", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), tags=list(), matrix=matrix()){
  new("PFMatrix", ID=ID, name=name, matrixClass=matrixClass, strand=strand, bg=bg,
      tags=tags,
      matrix=matrix)
}

ICMatrix = function(ID="Unknown", name="Unknown", matrixClass="Unknown",
                    strand="*", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), tags=list(), matrix=matrix(),
                    pseudocounts=numeric(), schneider=logical()){
  new("ICMatrix", ID=ID, name=name, matrixClass=matrixClass, strand=strand, bg=bg,
      tags=tags,
      matrix=matrix, pseudocounts=pseudocounts, schneider=schneider)
}
PWMatrix = function(ID="Unknown", name="Unknown", matrixClass="Unknown",
                    strand="*", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), tags=list(), matrix=matrix(),
                    pseudocounts=numeric()){
  new("PWMatrix", ID=ID, name=name, matrixClass=matrixClass, strand=strand, bg=bg,
      tags=tags,
      matrix=matrix, pseudocounts=pseudocounts)
}

### -----------------------------------------------------------------
### The set of PWM class
###
setClass("PFMatrixList",
         contains="SimpleList",
         representation(
                        ),
         prototype(
                   elementType="PFMatrix"
                   )
         )

setClass("PWMatrixList",
         contains="SimpleList",
         representation(
                        ),
         prototype(
                   elementType="PWMatrix"
                   )
         )

setClass("ICMatrixList",
         contains="PWMatrixList",
         representation(
                        ),
         prototype(
                   elementType="ICMatrix"
                   )
         )

setClassUnion("XMatrixList", c("PFMatrixList", "ICMatrixList", "PWMatrixList"))

### -----------------------------------------------------------------------
### XMatrixList() constructor.
###

setMethod("XMatrixList", "list",
          function(x, use.names=TRUE, type, ...){
            ok = sapply(x, is, "XMatrix")
            if(!all(ok))
              stop("XMatrixList() only accepts XMatrix objects!")
            if(!use.names)
              names(x) = NULL
            IRanges:::newList(type, x)
          }
          )

PFMatrixList = function(..., use.names=TRUE){
  listData = list(...)
  XMatrixList(listData, use.names=use.names, type="PFMatrixList")
}

PWMatrixList = function(..., use.names=TRUE){
  listData = list(...)
  XMatrixList(listData, use.names=use.names, type="PWMatrixList")
}

ICMatrixList = function(..., use.names=TRUE){
  listData = list(...)
  XMatrixList(listData, use.names=use.names, type="ICMatrixList")
}

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
### The Site constructor
###
Site = function(views, score, strand="*",
                seqname="Unknown",
                sitesource="TFBS", primary="TF binding site",
                pattern){
  new("Site", views=views, seqname=seqname, score=score, strand=strand,
      sitesource=sitesource, primary=primary, pattern=pattern)
}


### --------------------------------------------------------------
### SiteList objects
###

setClass("SiteList",
         contains="SimpleList",
         representation(
                        ),
         prototype(
                   elementType="Site"
                   )
         )
### -------------------------------------------------------------
### SiteList() constructor
###
SiteList = function(..., use.names=TRUE){
  listData = list(...)
  if(is(listData[[1]], "list"))
    listData = listData[[1]]
  ok = sapply(listData, is, "Site")
  if(!all(ok))
    stop("SiteList() only accepts Site objects!")
  if(!use.names)
    names(listData) = NULL
  IRanges:::newList("SiteList", listData)
}

### ----------------------------------------------------------------
### Methods
###
setMethod("writeGFF3", "SiteList",
          function(x){
            ans = do.call(rbind, lapply(x, writeGFF3))
            return(ans)
          }
          )
setMethod("writeGFF2", "SiteList",
           function(x){
             ans = do.call(rbind, lapply(x, writeGFF2))
             return(ans)
           }
           )


### -------------------------------------------------------------------
### SitePair object: a nucleotide sequence feature object representing (possibly putative) transcription factor binding site from A alignment

setClass("SitePair",
         slots=c(site1="Site",
                 site2="Site"
                 )
         )

### -----------------------------------------------------------------
### The SitePair constructor
###
SitePair = function(site1, site2){
  new("SitePair", site1=site1, site2=site2)
}


### ----------------------------------------------------------------
### SitePairList obejct: holds the list of SitePair.
###

setClass("SitePairList",
         contains="SimpleList",
         representation(
                        ),
         prototype(
                   elementType="SitePair"
                   )
         )

### ------------------------------------------------------------
### SitePairList constructor
###
SitePairList = function(..., use.names=TRUE){
  listData = list(...)
  if(is(listData[[1]], "list")) # This is pretty ugly. better solution?
    listData = listData[[1]]
  ok = sapply(listData, is, "SitePair")
  if(!all(ok))
    stop("SitePairList() only accepts SitePair objects!")
  if(!use.names)
    names(listData) = NULL
  IRanges:::newList("SitePairList", listData)
}

