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

