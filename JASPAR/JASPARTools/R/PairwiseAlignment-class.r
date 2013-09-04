### ----------------------------------------------------------------
### The PairwiseAlignmentTFBS object
### Let's define our PairwiseAlignmentTFBS object, based on the PairwiseAlignments object from Biostrings. For simplicity, make it as a slot.


setClass("PairwiseAlignmentTFBS",
         slots=c(alignments="PairwiseAlignments",
                 seqname1="character",
                 seqname2="character",
                 conservation1="numeric",
                 #conservation2="numeric", # because conservation2 is never used.
                 windowSize="integer",
                 cutoff="numeric",
                 seq1length="integer",
                 seq2length="integer"
                 )
         )

### ---------------------------------------------------------------
### The accessor-like methods
###
setMethod("alignments", "PairwiseAlignmentTFBS",
          function(x) x@alignments)
setMethod("seqname", "PairwiseAlignmentTFBS",
          function(x) c(x@seqname1, x@seqname2))
setMethod("conservation1", "PairwiseAlignmentTFBS",
          function(x) x@conservation1)
setMethod("seqlength", "PairwiseAlignmentTFBS",
          function(x) c(x@seq1length, x@seq2length))
setMethod("alnlength", "PairwiseAlignmentTFBS",
          function(x) nchar(alignments(x)))

### ----------------------------------------------------------------
### The constructor
###
PairwiseAlignmentTFBS = function(pattern, subject, type="global", 
                                 substitutionMatrix=NULL, gapOpening=0,
                                 gapExtension=-1,
                                 seqname1="Unknown", seqname2="Unknown",
                                 windowSize=51L, cutoff=0.7){
  alignments = PairwiseAlignments(pattern, subject, type=type, 
                                  substitutionMatrix=substitutionMatrix,
                                  gapOpening=gapOpening, gapExtension=gapExtension)
  conservation1 = calConservation(as.character(pattern(alignments)), as.character(subject(alignments)), windowSize=windowSize)
  seq1length = nchar(gsub("(-|_|\\.)", "", as.character(pattern(alignments))))
  seq2length = nchar(gsub("(-|_|\\.)", "", as.character(subject(alignments))))
  new("PairwiseAlignmentTFBS", alignments=alignments, seqname1=seqname1,
      seqname2=seqname2, conservation1=conservation1, windowSize=windowSize,
      cutoff=cutoff, seq1length=seq1length, seq2length=seq2length)
}

### ---------------------------------------------------------------
### The "show" method
### Add later... what is the pretty way?



