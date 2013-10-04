### -----------------------------------------------------------------
### Getters
###
setMethod("ID", "XMatrix", function(x) x@ID)
setMethod("name", "XMatrix", function(x) x@name)
setMethod("matrixClass", "XMatrix", function(x) x@matrixClass)
setMethod("Matrix", "XMatrix", function(x) x@matrix)
setMethod("strand", "XMatrix", function(x) x@strand)
setMethod("bg", "XMatrix", function(x) x@bg)
setMethod("tags", "XMatrix", function(x) x@tags)
setMethod("matrixType", "PFMatrix", function(x) "PFM")
setMethod("matrixType", "ICMatrix", function(x) "ICM")
setMethod("matrixType", "PWMatrix", function(x) "PWM")
setMethod("pseudocounts", "PWMatrix", function(x) x@pseudocounts)
setMethod("schneider", "ICMatrix", function(x) x@schneider)

### -----------------------------------------------------------------
### Setters
###
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
setReplaceMethod("Matrix", "XMatrix",
                 function(x, value){
                   x@matrix = value
                   return(x)
                 }
                 )
setReplaceMethod("pseudocounts", "PWMatrix",
                 function(x, value){
                   x@pseudocounts = value
                   return(x)
                 }
                 )
setReplaceMethod("schneider", "ICMatrix",
                 function(x, value){
                   x@schneider = value
                   return(x)
                 }
                 )

