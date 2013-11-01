
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
         function(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01,
                  max.results=10, min.percent_score=NULL, min.score=NULL){
           score = compareMatrix(pfmSubject, pfmQuery, openPenalty=openPenalty,
                                 extPenalty=extPenalty)
           relScore = 100 * score / max(ncol(pfmSubject), ncol(pfmQuery)) / 2
           return(c(score=score, relScore=relScore))
         }
         )

setMethod("searchMatrix", signature(pfmSubject="PFMatrix", pfmQuery="PFMatrix"),
          function(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01,
                   max.results=10, min.percent_score=NULL, min.score=NULL){
            score = compareMatrix(Matrix(pfmSubject), Matrix(pfmQuery), 
                                  openPenalty=openPenalty, extPenalty=extPenalty)
            relScore = 100 * score / max(ncol(Matrix(pfmSubject)), ncol(Matrix(pfmQuery))) / 2
            return(c(score=score, relScore=relScore))
          }
          )

