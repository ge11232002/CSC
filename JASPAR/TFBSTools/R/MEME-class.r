
### ----------------------------------------------------------------
### The MEME object which holds the result of a MEME run
###
setClass("MEME", 
         slots=c(
                 version="character",
                 alphabet="character",
                 command="character",
                 motifs="Motifs"
                 )
         )

### -----------------------------------------------------------------
### The constructor
###
MEME = function(version=character(), alphabet=c("A", "C", "G", "T"),
                command=character(), motifs){
  new("MEME", version=version, alphabet=alphabet, command=command, motifs=motifs)
}


