

### ---------------------------------------------------------------
### The real wrapper function for MEME
###
.read_version = function(memeOutput, fromLine){
  for(i in fromLine:length(memeOutput)){
    if(grepl("^MEME version", memeOutput[i]))
      break
  }
  memeVersion = strsplit(memeOutput[i], " ")[[1]][3]
  return(list(memeVersion, i+1))
}
.read_command = function(memeOutput, fromLine){
  for(i in fromLine:length(memeOutput)){
    if(grepl("^command: ", memeOutput[i]))
      break
  }
  memeCommand = gsub("^command: ", "", memeOutput[i])
  return(list(memeCommand, i +1))
}

# inputFastaFn = "/Users/gtan/src/meme_4.9.1/tests/crp0.s"

run_MEME = function(inputFastaFn, binary="meme", seqtype="DNA", arguments=""){
  cmd = paste(binary, inputFastaFn, "-text", ifelse(seqtype=="DNA", "-dna", "-protein"), arguments, "2>/dev/null")
  memeOutput = my.system(cmd, intern=TRUE)
  #conMemeOutput = pipe(cmd, open="rt")
  #on.exit(close(conMemeOutput))
  
  # get the version of meme used.
  oneLine = memeOutput[grep("^MEME version", memeOutput)]
  version = strsplit(oneLine, " ")[[1]][3]

  # get the command
  oneLine = memeOutput[grep("^command:", memeOutput)]
  command = gsub("^command: ", "", oneLine)

  
  # get the motifs information
  revcomp = grepl("revcomp", command)
  oneLines = memeOutput[grep("^MOTIF  \\d+", memeOutput)]
  splittedLines = strsplit(oneLines, "[[:blank:]]+")
  motifWidths = as.integer(sapply(splittedLines, "[", 5))
  motifOccurrences = as.integer(sapply(splittedLines, "[", 8))
  motifEvalues = as.numeric(sapply(splittedLines, "[", 14))
  
  # get the motifs names
  indexNames = grep("sorted by position p-value", memeOutput)
  oneLines = memeOutput[indexNames]
  motifNames = lapply(strsplit(oneLines, "[[:blank:]]+"), "[", c(2,3))
  motifNames = sapply(motifNames, paste0, collapse=" ")

  # get the motifs ranges
  motifList = list()
  for(i in seq_len(length(indexNames))){
    oneLines = memeOutput[seq(from=indexNames[i]+4, to=indexNames[i]+4+motifOccurrences[i]-1)]
    splittedLines = strsplit(oneLines, "[[:blank:]]+")
    oneRange = GRanges(seqnames=sapply(splittedLines, "[", 1), 
                       ranges=IRanges(start=as.integer(sapply(splittedLines, "[", ifelse(revcomp, 3, 2))), width=motifWidths[i]),
                       strand=ifelse(revcomp, sapply(splittedLines, "[", 2), "+")
                       )
    motifList = c(motifList, oneRange)
  }
  motifList = GRangesList(motifList)
  ans = list(motifList=motifList, motifEvalues=motifEvalues) 
  return(ans)
}



### -------------------------------------------------------------
### The MEME method
###
setMethod("runMEME", "character",
          function(x, binary="meme", seqtype="DNA", arguments="", tmpdir=tempdir()){
            seqtype = match.arg(seqtype, c("DNA", "AA"))
            ans = run_MEME(x, binary=binary, seqtype=seqtype, arguments=arguments)
            subjectSeqs = switch(seqtype,
                                 "DNA"=readDNAStringSet(filepath=x, format="fasta",),
                                 "AA"=readAAStringSet(filepath=x, format="fasta")
                                 )
            ans = MotifSet(motifList=ans[["motifList"]], motifEvalues=ans[["motifEvalues"]], subjectSeqs=subjectSeqs)
            return(ans)
          }
          )

setMethod("runMEME", "DNAStringSet",
          function(x, binary="meme", seqtype="DNA", arguments="", tmpdir=tempdir()){
            tmpFile = tempfile(pattern="MEME_", tmpdir=tmpdir, fileext = ".fasta")
            writeXStringSet(x, filepath=tmpFile, format="fasta")
            on.exit(unlink(tmpFile))
            ans = run_MEME(tmpFile, binary=binary, seqtype=seqtype, arguments=arguments, tmpdir=tmpdir)
            ans = MotifSet(motifList=ans[["motifList"]], motifEvalues=ans[["motifEvalues"]], subjectSeqs=x)
            return(ans)
          }
          )




