

lastz = function(assemblyTarget, assemblyQuery, chrsTarget=NULL, chrsQuery=NULL, distance="medium", format="lav", echoCommand=FALSE){
  require(rtracklayer)
  require(R.utils)
# This matrix is taken from http://genomewiki.ucsc.edu/index.php/GorGor3_conservation_lastz_parameters. Default HOXD70 is medium. HoxD55 is far. human-chimp.v2 is close.
  lastzMatrix = list(medium = matrix(c(91, -114, -31, -123,
                                      -114, 100, -125, -31,
                                      -31, -125, 100,-114,
                                      -123, -31, -114, 91),
                                      nrow=4, ncol=4,
                                      dimnames=list(c("A", "C", "G", "T"),
                                                    c("A", "C", "G", "T"))
                                  ),
                     far = matrix(c(91, -90, -25, -100,
                                    -90, 100, -100, -25,
                                    -25, -100, 100, -90,
                                    -100, -25, -90, 91),
                                    nrow=4, ncol=4,
                                    dimnames=list(c("A", "C", "G", "T"),
                                                  c("A", "C", "G", "T"))
                                  ),
                     near = matrix(c(90, -330, -236, -356,
                                    -330, 100, -318, -236,
                                    -236, -318, 100, -330,
                                    -356, -236, -330, 90),
                                    nrow=4, ncol=4,
                                    dimnames=list(c("A", "C", "G", "T"),
                                                  c("A", "C", "G", "T"))
                                  )
                     )
  matrixFile = paste(Sys.getpid(), "-lastzMatrix.dat", sep="")
## The options used here is taken from RunLastzChain_sh.txt genomewiki.ucsc.edu. http://genomewiki.ucsc.edu/images/9/93/RunLastzChain_sh.txt.
  # B=0, --strand=plus: Search the forward strand only (the one corresponding to the query specifier). By default, both strands are searched.
  # C=0, --nochain: Skip the chaining stage. By default, the chaining stage is skipped.
  # E=30,150;O=600,400: Set the score penalties for opening and extending a gap.
  # H=0, 2000; --inner=<score>: Perform additional alignment between the gapped alignment blocks, using (presumably) more sensitive alignment parameters. By default this is not performed.
  # K=4500,3000,2200; --hspthresh=<score>:  Set the score threshold for the x-drop extension method; HSPs scoring lower are discarded. By default, use the entropy adjustment.
  # L=3000,6000; --gappedthresh=<score>: Set the threshold for gapped extension; alignments scoring lower than <score> are discarded.  By default gapped extension is performed, and alignment ends are trimmed to the locations giving the maximum score. 
  # M=254,50,50;--masking=<count>:  Dynamically mask the target sequence by excluding any positions that appear in too many alignments from further consideration for seeds. By default, a step of 1 is used, no words are removed from the target seed word position table, dynamic masking is not performed, and no target capsule or segment file is used.
  # T=1,2;--seed=12of19:  Seeds require a 19-bp word with matches in 12 specific positions (1110100110010101111). By default the 12-of-19 seed is used, one transition is allowed (except with quantum DNA), the hits are not filtered, twins are not required, and hash collisions are not recovered. 
  #  Y=15000,9400,3400; --ydrop=<dropoff>: Set the threshold for terminating gapped extension; this restricts the endpoints of each local alignment by limiting the local region around each anchor in which extension is performed.
  lastzOptions = list(near=paste("C=0 E=150 H=0 K=4500 L=3000 M=254 O=600 T=2 Y=15000 Q=", matrixFile, sep=""),
                      medium="C=0 E=30 H=0 K=3000 L=3000 M=50 O=400 T=1 Y=9400",
                      far=paste("C=0 E=30 H=2000 K=2200 L=6000 M=50 O=400 T=2 Y=3400 Q=", matrixFile, sep="")
                      )
  write.table(lastzMatrix[[distance]], file=matrixFile, quote=FALSE,
             sep=" ", row.names=FALSE, col.names=TRUE)
  message("lastz")
  if(is.null(chrsTarget)){
    chrsTarget = seqnames(seqinfo(TwoBitFile(assemblyTarget)))
  }else{
    stopifnot(all(chrsTarget %in% seqnames(seqinfo(TwoBitFile(assemblyTarget)))))
  }
  if(is.null(chrsQuery)){
    chrsQuery = seqnames(seqinfo(TwoBitFile(assemblyQuery)))
  }else{
    stopifnot(all(chrsQuery %in% seqnames(seqinfo(TwoBitFile(assemblyQuery)))))
  }
  chrTarget = chrsTarget[1]
  chrQuery = chrsQuery[1]
  outputToReturn = c()
  for(chrTarget in chrsTarget){
    for(chrQuery in chrsQuery){
      ## Deal the "|" in chr name
      output = paste(gsub("|", "\\|", chrTarget, fixed=TRUE), 
                     ".", sub("\\..*$", "", basename(assemblyTarget)),
                     "-", gsub("|", "\\|", chrQuery, fixed=TRUE),
                     ".", sub("\\..*$", "", basename(assemblyQuery)),
                     ".", format, sep="")
      outputToReturn = c(outputToReturn, output)
      cmd = paste("lastz", " ", assemblyTarget, "/", 
                  gsub("|", "\\|", chrTarget, fixed=TRUE), " ", 
                  assemblyQuery, "/", 
                  gsub("|", "\\|", chrQuery, fixed=TRUE), " ", 
                  lastzOptions[[distance]], 
                  " --format=", format, 
                  " --output=", output,
                  " --markend",
                  sep="")
      if(echoCommand){
        message(cmd)
      }else{
        if(!file.exists(output)){
          my.system(cmd)
        }else{
          warning("The output ", output, " already exists! Skipping..")
        }
      }
      #if(countLines(output) == 15L){
        ## delete the empty files
      #  unlink(output)
      #}
    }
  }
  unlink(matrixFile)
  return(outputToReturn)
}

validateLastz = function(lavs){
  filesNotCompleted = c()
  lav = lavs[1]
  for(lav in lavs){
    lastLine = my.system(paste("tail -n 1", lav), intern=TRUE)
    if(lastLine != "# lastz end-of-file"){
      filesNotCompleted = c(filesNotCompleted, lav)
    }
  }
  return(filesNotCompleted)
}


lavToPsl = function(lavs, psls=sub("\\.lav$", ".psl", lavs, ignore.case=TRUE), 
                    removeLav=TRUE){
  #for(i in 1:length(lavs)){
  #  cmd = paste("lavToPsl", lavs[i], psls[i])
  #  my.system(cmd)
  #}
  tempFile <- tempfile(pattern="lavToPsl", tmpdir=".")
  writeLines(paste("lavToPsl", lavs, psls), con=tempFile)
  my.system(paste("sh", tempFile))
  unlink(tempFile)
  if(removeLav){
    unlink(lavs)
  } 
  return("success")
}

axtChain = function(inputs, assemblyTarget, assemblyQuery, format="axt", 
                    outputs=sub(paste("\\.", format, "$", sep=""), ".chain", inputs, ignore.case=TRUE), distance="near", removePsl=TRUE){
  chainOptions =list(near="-minScore=5000 -linearGap=medium",
                   medium="-minScore=3000 -linearGap=medium",
                   far="-minScore=5000 -linearGap=loose"
                   )
  lastzMatrix = list(medium = matrix(c(91, -114, -31, -123,
                           -114, 100, -125, -31,
                           -31, -125, 100,-114,
                           -123, -31, -114, 91),
                        nrow=4, ncol=4,
                        dimnames=list(c("A", "C", "G", "T"),
                                      c("A", "C", "G", "T"))
                        ),
                     far = matrix(c(91, -90, -25, -100,
                       -90, 100, -100, -25,
                       -25, -100, 100, -90,
                       -100, -25, -90, 91),
                     nrow=4, ncol=4,
                     dimnames=list(c("A", "C", "G", "T"),
                                   c("A", "C", "G", "T"))
                     ),
                     near = matrix(c(90, -330, -236, -356,
                        -330, 100, -318, -236,
                        -236, -318, 100, -330,
                        -356, -236, -330, 90),
                      nrow=4, ncol=4,
                      dimnames=list(c("A", "C", "G", "T"),
                                    c("A", "C", "G", "T"))
                      )
                     )
  matrixFile = paste(Sys.getpid(), "-lastzMatrix.dat", sep="")
  write.table(lastzMatrix[[distance]], file=matrixFile, quote=FALSE, 
              sep=" ", row.names=FALSE, col.names=TRUE)
  #for(i in 1:length(inputs)){
  #  cmd = "axtChain"
  #  if(format == "psl"){
  #    cmd = paste(cmd, "-psl")
  #  }
  #  cmd = paste(cmd, chainOptions[[distance]]) 
  #  cmd = paste(cmd, " -scoreScheme=", matrixFile, sep="")
  #  cmd = paste(cmd, inputs[i], assemblyTarget, assemblyQuery, outputs[i])
  #  my.system(cmd)
  #}
  cmd = paste0("axtChain -psl ", chainOptions[[distance]], " -scoreScheme=",
               matrixFile, " ", inputs, " ", assemblyTarget, 
               " ", assemblyQuery, " ", outputs)
  tempFile <- tempfile(pattern="axtChain", tmpdir=".")
  writeLines(cmd, con=tempFile)
  my.system(paste("sh", tempFile))
  unlink(tempFile)
  unlink(matrixFile)
  if(removePsl){
    unlink(inputs)
  }
  return(outputs)
}

chainMergeSort = function(path="chain", assemblyTarget, assemblyQuery, removeChains=TRUE){
  allChain = paste(sub("\\.2bit$", "", basename(assemblyTarget), ignore.case=TRUE), 
                   ".",
                   sub("\\.2bit$", "", basename(assemblyQuery), ignore.case=TRUE),
                   ".all.chain", sep=""
                   )
  cmd = paste0("find ", path, " -name \"*.chain\" | chainMergeSort -inputList=stdin > ", allChain)
  my.system(cmd)
  if(removeChains){
    unlink(chains)
  }
  return(allChain)
}

genomeSizesFrom2Bit = function(twoBitFile, output){
  require(rtracklayer)
  genome.sizes = cbind(seqnames(seqinfo(TwoBitFile(twoBitFile))),
                       seqlengths(seqinfo(TwoBitFile(twoBitFile))))
  write.table(genome.sizes, file=output, row.names=FALSE, col.names=FALSE, quote=FALSE)
}


chainPreNet = function(allChain, assemblyTarget, assemblyQuery, removeAllChain=FALSE){
  target.sizesFile = paste(Sys.getpid(), "-target.sizes.txt", sep="")
  genomeSizesFrom2Bit(assemblyTarget, target.sizesFile)
  query.sizesFile = paste(Sys.getpid(), "-query.sizes.txt", sep="")
  genomeSizesFrom2Bit(assemblyQuery, query.sizesFile)
  allPreChain = paste(sub("\\.2bit$", "", basename(assemblyTarget), ignore.case=TRUE),
                      ".",
                      sub("\\.2bit$", "", basename(assemblyQuery), ignore.case=TRUE),
                      ".all.pre.chain", sep=""
                      )
  cmd = "chainPreNet"
  cmd = paste(cmd, allChain, target.sizesFile, query.sizesFile, allPreChain)
  my.system(cmd)
  unlink(c(target.sizesFile, query.sizesFile))
  if(removeAllChain){
    unlink(allChain)
  }
  return(allPreChain)
}

chainNetSyntenic = function(allPreChain, assemblyTarget, assemblyQuery){
  cmd = "chainNet"
  target.sizesFile = paste(Sys.getpid(), "-target.sizes.txt", sep="")
  genomeSizesFrom2Bit(assemblyTarget, target.sizesFile)
  query.sizesFile = paste(Sys.getpid(), "-query.sizes.txt", sep="")
  genomeSizesFrom2Bit(assemblyQuery, query.sizesFile)
  target.net = paste(Sys.getpid(), "-target.net", sep="")
  query.net  = paste(Sys.getpid(), "-query.net", sep="")
  cmd = paste(cmd, allPreChain, target.sizesFile, query.sizesFile, target.net, query.net)
  my.system(cmd)
  unlink(c(target.sizesFile, query.sizesFile))

  cmd = "netSyntenic"
  netSyntenicFile = paste(sub("\\.2bit$", "", basename(assemblyTarget), ignore.case=TRUE),
                          ".",                                                      
                          sub("\\.2bit$", "", basename(assemblyQuery), ignore.case=TRUE),            
                          ".noClass.net", sep=""
                          )
  cmd = paste(cmd, target.net, netSyntenicFile)
  my.system(cmd)
  unlink(c(target.net, query.net))
  return(netSyntenicFile)
}

netToAxt = function(in.net, in.chain, assemblyTarget, assemblyQuery, removeFiles=FALSE){
  cmd = "netToAxt"
  axtFile = paste(sub("\\.2bit$", "", basename(assemblyTarget), ignore.case=TRUE),
                  ".",
                  sub("\\.2bit$", "", basename(assemblyQuery), ignore.case=TRUE),
                  ".net.axt", sep=""
                  )
  cmd = paste(cmd, in.net, in.chain, assemblyTarget, assemblyQuery, "stdout |",
              "axtSort", "stdin", axtFile)
  my.system(cmd)
  if(removeFiles){
    unlink(c(in.net, in.chain))
  }
  return(axtFile)
}

axtNetWholePipeline = function(lavs, assemblyTarget, assemblyQuery, distance, removeFiles=TRUE){
  ## lastz might exit without any errors. check the files whether are complete.
  filesNotCompleted = validateLastz(lavs)
  stopifnot(is.null(filesNotCompleted))
### step 2: Chaining
  psls = lavToPsl(lavs, removeLav=removeFiles)
  chains = axtChain(psls, assemblyTarget, assemblyQuery, format="psl", distance=distance, removePsl=removeFiles)
  allChain = chainMergeSort(chains, assemblyTarget, assemblyQuery, removeChains=removeFiles)
# now , this allChain is what we get from UCSC download e.g. hg19.danRer7.all.chain.gz
### step 3: Netting
  allPreChain = chainPreNet(allChain, assemblyTarget, assemblyQuery, removeAllChain=removeFiles)
# combine the chains into nets, add the synteny information to the net:
  netSyntenicFile = chainNetSyntenic(allPreChain, assemblyTarget, assemblyQuery)

#  Add classification information using the database tables:
# actually this step is not necessary in this pipeline according to http://blog.gmane.org/gmane.science.biology.ucscgenome.general/month=20130301. The class information will only be used for Genome Browser.
# Since it needs some specific modification of the table names for certain species, we skip this step now.
# If this step is done, then the generated class.net is the gzipped net file that you see in UCSC Downloads area.

### step 4: axtNet
  axtFile = netToAxt(netSyntenicFile, allPreChain, assemblyTarget, assemblyQuery, removeFiles=removeFiles)
### gzip the axtNet file
  my.system(paste("gzip", axtFile))
  return("Success")
}

monitorBatchJobs = function(reg, sleep=3600){
  require(BatchJobs)
  statusBatchJobs = showStatus(reg)
## we check whetehr the jobs are finished or not every hour.
  while(statusBatchJobs[1,"running"] != 0 || statusBatchJobs[1,"expired"] != 0 || statusBatchJobs[1,"n"] != statusBatchJobs[1,"started"]){
    message(statusBatchJobs[1,"submitted"] - statusBatchJobs[1,"done"], " jobs to do!")
    if(statusBatchJobs[1,"expired"] != 0){
      resubmitIds = findExpired(reg)
      submitJobs(reg, resubmitIds)
    }
    Sys.sleep(sleep)
    statusBatchJobs = showStatus(reg)
    message(statusBatchJobs[1,"submitted"] - statusBatchJobs[1,"done"], " jobs to do!")
  }
  return("Success")
}

