
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


