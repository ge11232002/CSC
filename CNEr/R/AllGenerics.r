
### -----------------------------------------------------------------
### Axt object related
###
setGeneric("targetRanges", function(x) standardGeneric("targetRanges"))
setGeneric("targetSeqs", function(x) standardGeneric("targetSeqs"))
setGeneric("queryRanges", function(x) standardGeneric("queryRanges"))
setGeneric("querySeqs", function(x) standardGeneric("querySeqs"))
setGeneric("symCount", function(x) standardGeneric("symCount"))
setGeneric("subAxt", function(x, chr, start, end, #strand=c("+", "-", "*"),
                              select=c("target", "query"),
                              type=c("any", "within"),
                              qSize=NULL) 
                      standardGeneric("subAxt")
                      )
### -----------------------------------------------------------------
### general
### Not Exported!
setGeneric("clone", function(x, ...) standardGeneric("clone"))

### -----------------------------------------------------------------
### ceScan
### Exported!
setGeneric("ceScan", function(axts, tFilter, qFilter, qSizes, thresholds="49,50")
           standardGeneric("ceScan"))


