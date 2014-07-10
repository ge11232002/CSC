

### ------------------------------------------------------------------ 
### The RepeatMaske rwrapper function for splitted assembly
###
repeatMasker <- function(inLstFn, species="danio", engine="crossmatch",
                         binary="/usr/local/bin/RepeatMasker/RepeatMasker"){
  engine <- match.arg(engine, c("crossmatch", "wublast", "abblast", "ncbi",
                                "hmmer", "decypher"))
  finalOutFn <- paste0(inLstFn, ".out")
  inLftFn <- sub("\\.lst$", ".lft", inLstFn)
  alignOutFn <- paste0(inLstFn, ".align")
  catOutFn <- paste0(inLstFn, ".cat")

  cwd <- getwd()
  tmpDir <- tempdir()
  dir.create(tmpDir)
  setwd(tmpDir)
  # Initialize local library
  cmd <- paste(binary, "-engine", engine, "-pa 1 -species", species, "/dev/null")
  CNEr:::my.system(cmd)

  # do the repeatMasker
  specs <- readLines(inLstFn)
  for(spec in specs){
    # Remove path and .2bit filename to get just the seq:start-end spec:
    # If $spec is the whole sequence, twoBitToFa removes the :start-end part,
    # which causes liftUp to barf later.  So tweak the header back to
    # seq:start-end for liftUp's sake:
    base <- sub("[^:]+:", "", spec)
    cmd <- paste("twoBitToFa", spec, "stdout | sed -e",
                 paste0("\"s/^>.*/>", base,"/\""), ">",
                 paste0(base, ".fa"))
    CNEr:::my.system(cmd)
    
    cmd <- paste(binary, "-engine", engine, " -pa 1 -align -species", 
                 species, paste0(base, ".fa"))
    CNEr:::my.system(cmd)
    if(file.exists(paste0(base, ".fa.cat"))){
      file.rename(paste0(base, ".fa.cat"), catOutFn)
    }
  }
  # Lift up (leave the RepeatMasker header in place because we'll liftUp
  # again later):
  faOutFns <- list.files(path=".", pattern=".*\\.fa\\.out")
  cmd <- paste("liftUp -type=.out stdout",  inLftFn, "error", 
               paste0(faOutFns, collapse=" "), "> tmpOut_out")
  CNEr:::my.system(cmd)
  alignFiles <- list.files(path=".", pattern=".*\\.align")
  if(length(alignFiles) > 0){
    cmd <- paste("liftUp -type=.align stdout", inLftFn, "error",
                 paste0(alignFiles, collapse=" "), "> tmpOut_align")
    CNEr:::my.system(cmd)
  }
  file.rename("tmpOut_out", finalOutFn)
  file.rename("tmpOut_align", alignOutFn)
  setwd(cwd)
  unlink(tmpDir, recursive=TRUE)
  return("success")
}

