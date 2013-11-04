
### -----------------------------------------------------------------
### searchMatrix method. compare two position frequency matrix.
###
compareMatrix = function(pfmSubject, pfmQuery, openPenalty, extPenalty){
  # The true aligning engine. Taking two ordinary matrixs.
  pfmSubject = normargPfm(pfmSubject)
  pfmQuery = normargPfm(pfmQuery)
  ans = .Call("matrixAligner", pfmSubject, pfmQuery, openPenalty, extPenalty)
  return(ans)
}

setMethod("searchMatrix", signature(pfmSubject="matrix", pfmQuery="matrix"),
         function(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01){
           score = compareMatrix(pfmSubject, pfmQuery, openPenalty=openPenalty,
                                 extPenalty=extPenalty)
           relScore = 100 * score / max(ncol(pfmSubject), ncol(pfmQuery)) / 2
           return(c(score=score, relScore=relScore))
         }
         )

setMethod("searchMatrix", signature(pfmSubject="PFMatrix", pfmQuery="PFMatrix"),
          function(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01){
            ans = searchMatrix(Matrix(pfmSubject), Matrix(pfmQuery), 
                               openPenalty=openPenalty, extPenalty=extPenalty)
            return(ans)
          }
          )

setMethod("searchMatrix", signature(pfmSubject="PFMatrix", pfmQuery="matrix"),
          function(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01){
            ans = searchMatrix(Matrix(pfmSubject), pfmQuery,
                               openPenalty=openPenalty, extPenalty=extPenalty)
            return(ans)
          }
          )

setMethod("searchMatrix", signature(pfmSubject="matrix", pfmQuery="PFMatrix"),
          function(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01){
            ans = searchMatrix(pfmSubject, Matrix(pfmQuery),
                               openPenalty=openPenalty, extPenalty=extPenalty)
            return(ans)
          }
          )

setMethod("searchMatrix", signature(pfmSubject="PFMatrixList", pfmQuery="matrix"),
          function(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01){
            ans = lapply(pfmSubject, searchMatrix, pfmQuery,
                            openPenalty=openPenalty, extPenalty=extPenalty)
            return(ans)
          }
          )

setMethod("searchMatrix", signature(pfmSubject="PFMatrixList", pfmQuery="PFMatrix"),
          function(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01){
            ans = lapply(pfmSubject, searchMatrix, pfmQuery,
                         openPenalty=openPenalty, extPenalty=extPenalty)
            return(ans)
          }
          )


