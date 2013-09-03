### -------------------------------------------------------------------
### XMatrixList objects
### Not clear why compressedList is not working. check later.

#setClass("PFMatrixList",
#         contains="CompressedList",
#         representation(
#                        unlistData="PFMatrix"
#                        ),
#         prototype(
#                   elementType="PFMatrix"
#                   )
#         )
setClass("PFMatrixList",
         contains="SimpleList",
         representation(
                        ),
         prototype(
                   elementType="PFMatrix"
                   )
         )


#setClass("PWMatrixList",
#         contains="CompressedList",
#         representation(
#                        unlistData="PWMatrix"
#                        ),
#         prototype(
#                   elementType="PWMatrix"
#                   )
#         )

setClass("PWMatrixList",
         contains="SimpleList",
         representation(
                        ),
         prototype(
                   elementType="PWMatrix"
                   )
         )
#setClass("ICMatrixList",
#         contains="CompressedList",
#         representation(
#                        unlistData="ICMatrix"
#                        ),
#         prototype(
#                   elementType="ICMatrix"
#                   )
#         )

setClass("ICMatrixList",
         contains="SimpleList",
         representation(
                        ),
         prototype(
                   elementType="ICMatrix"
                   )
         )

setClassUnion("XMatrixList", c("PFMatrixList", "ICMatrixList", "PWMatrixList"))

### ----------------------------------------------------------------------
### The accessor-like method
###


### -----------------------------------------------------------------------
### XMatrixList() constructor.
###

setGeneric("XMatrixList", signature="x",
           function(x, use.names=TRUE, ...)
             standardGeneric("XMatrixList")
           )
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





