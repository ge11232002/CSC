
### ---------------------------------------------------------
### Motifs object
### 
setClass("Motifs",
         slots=c(
                 motifList="GRangesList",
                 subjectSeqs="DNAStringSet"
                 )
         )

### ---------------------------------------------------------
### The constructor function
###
Motifs = function(motifList=GRangesList(), subjectSeqs=DNAStringSet()){
  new("Motifs", motifList=motifList, subjectSeqs=subjectSeqs)
}

