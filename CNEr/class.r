

setClass(Class="axt",
         representation(targetNames="Rle",
                        targetRanges="IRanges",
                        targetSeqs="DNAStringSet",
                        queryNames="Rle",
                        queryRanges="IRanges",
                        querySeqs="DNAStringSet",
                        strand="Rle",
                        score="rle"
                        ),
         contains="Vector",
         prototype(
                   targetNames=Rle(factor()),
                   queryNames=Rle(factor()),
                   strand=Rle(strand())
                   )
         )

showGenomicRanges = function(x, margin="",
                             with.classinfo=FALSE, print.seqlengths=FALSE){
  lx <- length(x)
  nc <- ncol(mcols(x))
  cat(class(x), " with ", lx, " ", ifelse(lx == 1L, "alignment", "alignments"),
      ":\n", sep="")

}

setMethod("show", "axt",
          function(object)
            showAxt(object, margin="  ",
                              with.classinfo=TRUE, print.seqlengths=TRUE)
)

