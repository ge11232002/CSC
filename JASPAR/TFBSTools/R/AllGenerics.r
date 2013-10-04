## Accessors
setGeneric("ID", signature="x", function(x) standardGeneric("ID"))
setGeneric("ID<-", signature="x", function(x, value) standardGeneric("ID<-"))
setGeneric("name", signature="x", function(x) standardGeneric("name"))
setGeneric("name<-", signature="x", function(x, value) standardGeneric("name<-"))
setGeneric("matrixClass", signature="x", function(x) standardGeneric("matrixClass"))
setGeneric("matrixClass<-", signature="x", function(x, value) standardGeneric("matrixClass<-"))
setGeneric("Matrix", signature="x", function(x) standardGeneric("Matrix"))
setGeneric("Matrix<-", signature="x", function(x, value) standardGeneric("Matrix<-"))
setGeneric("bg", signature="x", function(x) standardGeneric("bg"))
setGeneric("bg<-", signature="x", function(x, value) standardGeneric("bg<-"))
setGeneric("tags", signature="x", function(x) standardGeneric("tags"))

setGeneric("matrixType", signature="x", function(x) standardGeneric("matrixType"))
setGeneric("pseudocounts", signature="x", function(x) standardGeneric("pseudocounts"))
setGeneric("pseudocounts<-", signature="x", function(x, value) standardGeneric("pseudocounts<-"))
setGeneric("schneider", signature="x", function(x) standardGeneric("schneider"))
setGeneric("schneider<-", signature="x", function(x, value) standardGeneric("schneider<-"))

## Constructors
setGeneric("XMatrixList", signature="x",
           function(x, use.names=TRUE, ...)
             standardGeneric("XMatrixList")
           )


# JASPAR DB
setGeneric("getMatrixByID", signature="x", function(x, ID) standardGeneric("getMatrixByID"))
setGeneric("getMatrixByName", signature="x", function(x, name)  standardGeneric("getMatrixByName"))
setGeneric("getMatrixSet", signature="x", function(x, opts) standardGeneric("getMatrixSet"))
setGeneric("storeMatrix",
           function(x, pfmList) standardGeneric("storeMatrix")
           )
setGeneric("initializeJASPARDB", signature="x",
           function(x) standardGeneric("initializeJASPARDB"))
setGeneric("deleteMatrixHavingID", signature="x",
           function(x, IDs) standardGeneric("deleteMatrixHavingID"))


