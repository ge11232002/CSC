

### ---------------------------------------------------------------
### The real wrapper function for MEME
###
run_MEME = function(inputFastaFn, binary="meme", arguments=""){
   

}



### -------------------------------------------------------------
### The MEME method
###
setMethod("MEME", "character",
          function(x, binary="meme", arguments="", tmpdir=tempdir()){
            run_MEME(x, binary=binary, arguments=arguments)
          }
          )

setMethod("MEME", "DNAStringSet",
          function(x, binary="meme", arguments="", tmpdir=tempdir()){
            tmpFile = tempfile(pattern="MEME_", tmpdir=tmpdir, fileext = ".fasta")
            writeXStringSet(x, filepath=tmpFile, format="fasta")
            MEME(tmpFile, binary=binary, arguments=arguments, tmpdir=tmpdir)
          }
          )

