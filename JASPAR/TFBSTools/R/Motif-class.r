
### ---------------------------------------------------------
### Motifs object
### 
setClass("Motifs",
         slots=c(
                 motifList="GRangesList",
                 motifEvalues="numeric",
                 subjectSeqs="DNAStringSet"
                 )
         )

### ---------------------------------------------------------
### The constructor function
###
Motifs = function(motifList=GRangesList(), motifEvalues=numeric(), subjectSeqs=DNAStringSet()){
  new("Motifs", motifList=motifList, motifEvalues=motifEvalues, 
      subjectSeqs=subjectSeqs)
}



### -------------------------------------------------------
### Methods
###
setMethod("sitesSeq", "Motifs",
          function(x){

          }
          )



