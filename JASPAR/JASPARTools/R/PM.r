
### ------------------------------------------------------------------------
### The generic position matrix objects.
setClass("XMatrix", contains=c("matrix", "VIRTUAL"),
         slots=c(ID="character",
                 name="character",
                 matrixClass="character",
                 strand="character",
                 bg="numeric"
                 ))

### XMatrix subclasses without additional slots
setClass("PFMatrix", contains="XMatrix")
setClass("ICMatrix", contains="XMatrix")
setClass("PWMatrix", contains="XMatrix")


### ----------------------------------------------------------------------
### The accessor-like method
setGeneric("ID", signature="x", function(x) standardGeneric("ID"))
setMethod("ID", "XMatrix", function(x) x@ID)
setGeneric("name", signature="x", function(x) standardGeneric("name"))
setMethod("name", "XMatrix", function(x) x@name)
setGeneric("matrixClass", signature="x", function(x) standardGeneric("matrixClass"))
setMethod("matrixClass", "XMatrix", function(x) x@matrixClass)
setMethod("strand", "XMatrix", function(x) x@strand)
setGeneric("bg", signature="x", function(x) standardGeneric("bg"))
setMethod("bg", "XMatrix", function(x) x@bg)

setReplaceMethod("ID", "XMatrix", 
                 function(x, value){
                   x@ID = value
                   return(x)
                 }
                 )
setReplaceMethod("name", "XMatrix",
                 function(x, value){
                   x@name = value
                   return(x)
                 }
                 )
setReplaceMethod("matrixClass", "XMatrix",
                 function(x, value){
                   x@matrixClass = value
                   return(x)
                 }
                 )
setReplaceMethod("strand", "XMatrix",
                 function(x, value){
                   x@strand = value
                   return(x)
                 }
                 )
setReplaceMethod("bg", "XMatrix",
                 function(x, value){
                   x@bg = value
                   return(x)
                 }
                 )

### -----------------------------------------------------------------------
### The "show" method.
###

setMethod("show", "XMatrix",
          function(object){
            
          }
          )


