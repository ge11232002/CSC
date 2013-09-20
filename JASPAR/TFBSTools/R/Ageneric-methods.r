
setGeneric("toICM", signature="x",
           function(x, pseudocounts=NULL, schneider=FALSE,
                    bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
             standardGeneric("toICM")
           )

setGeneric("toPWM", signature="x",
           function(x, type="log2probratio", pseudocounts=NULL,
                    bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
             standardGeneric("toPWM")
           )

setGeneric("plotLogo", signature="x",
           function(x, ic.scale = TRUE, xaxis = TRUE, yaxis = TRUE,
                    xfontsize = 15, yfontsize = 15) standardGeneric("plotLogo")
           )
setGeneric("total_ic", signature="x",
           function(x) standardGeneric("total_ic"))

setGeneric("searchSeq", signature="x",
           function(x, subject, seqname="Unknown", strand="*", min.score="80%")
             standardGeneric("searchSeq"))

setGeneric("searchAln", 
           function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                    conservation=NULL)
             standardGeneric("searchAln")
           )
setGeneric("doSiteSearch",
           function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                    conservation=NULL)
             standardGeneric("doSiteSearch")
           )

setGeneric("calConservation", 
           function(aln1, aln2, windowSize=51L, which="1")
             standardGeneric("calConservation")
           )

setGeneric("clone", signature="x", function(x, ...) standardGeneric("clone"))
setGeneric("revcom", signature="x", function(x) standardGeneric("revcom"))
setGeneric("writeGFF3", signature="x", function(x) standardGeneric("writeGFF3"))
setGeneric("writeGFF2", signature="x", function(x) standardGeneric("writeGFF2"))
setGeneric("relScore", signature="x", function(x) standardGeneric("relScore"))

## DB
setGeneric("get_Matrix_by_ID", signature="x", function(x, ID, type="PFM") standardGeneric("get_Matrix_by_ID"))
setGeneric("get_Matrix_by_name", signature="x", function(x, name, type="PFM") standardGeneric("get_Matrix_by_name"))
setGeneric("get_MatrixSet", signature="x", function(x, opts) standardGeneric("get_MatrixSet"))
setGeneric("store_Matrix",
           function(x, pfmList) standardGeneric("store_Matrix")
           )
setGeneric("initializeJASPARDB", signature="x",
           function(x) standardGeneric("initializeJASPARDB"))
setGeneric("delete_Matrix_having_ID", signature="x",
           function(x, IDs) standardGeneric("delete_Matrix_having_ID"))

## wrappers
setGeneric("runMEME", signature="x", function(x, binary="meme", seqtype="DNA", arguments="", tmpdir=tempdir()) standardGeneric("runMEME"))
setGeneric("sitesSeq", signature="x", function(x, n=10, type="none") standardGeneric("sitesSeq"))


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
setGeneric("views", signature="x", function(x) standardGeneric("views"))
setGeneric("seqname", signature="x", function(x) standardGeneric("seqname"))
setGeneric("sitesource", signature="x", function(x) standardGeneric("sitesource"))
setGeneric("primary", signature="x", function(x) standardGeneric("primary"))
setGeneric("site1", signature="x", function(x) standardGeneric("site1"))
setGeneric("site2", signature="x", function(x) standardGeneric("site2"))
setGeneric("alignments", signature="x", function(x) standardGeneric("alignments"))
setGeneric("conservation1", signature="x", function(x) standardGeneric("conservation1"))
setGeneric("seqlength", signature="x", function(x) standardGeneric("seqlength"))
setGeneric("alnlength", signature="x", function(x) standardGeneric("alnlength"))

setGeneric("XMatrixList", signature="x",
           function(x, use.names=TRUE, ...)
             standardGeneric("XMatrixList")
           )

