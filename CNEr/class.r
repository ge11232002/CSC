

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



