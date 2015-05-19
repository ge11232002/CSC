#########################################################################
# File Name: last.R
# Author: Ge Tan
# mail: gtan@me.com
# Created Time: Mon 28 Jul 2014 01:43:58 PM BST
#########################################################################

last <- function(db, queryFn, outputFn,
                 distance="medium", format=c("MAF","tabular"),
                 mc.cores=1L, echoCommand=FALSE){
  format <- match.arg(format)
  # This matrix is taken from http://genomewiki.ucsc.edu/index.php/GorGor3_conservation_lastz_parameters. Default HOXD70 is medium. HoxD55 is far. human-chimp.v2 is close.
  lastzMatrix <- list(medium=matrix(c(91, -114, -31, -123,
                                      -114, 100, -125, -31,
                                      -31, -125, 100,-114,
                                      -123, -31, -114, 91),
                                      nrow=4, ncol=4,
                                      dimnames=list(c("A", "C", "G", "T"),
                                                    c("A", "C", "G", "T"))
                                  ),
                     far=matrix(c(91, -90, -25, -100,
                                    -90, 100, -100, -25,
                                    -25, -100, 100, -90,
                                    -100, -25, -90, 91),
                                    nrow=4, ncol=4,
                                    dimnames=list(c("A", "C", "G", "T"),
                                                  c("A", "C", "G", "T"))
                                  ),
                     near=matrix(c(90, -330, -236, -356,
                                    -330, 100, -318, -236,
                                    -236, -318, 100, -330,
                                    -356, -236, -330, 90),
                                    nrow=4, ncol=4,
                                    dimnames=list(c("A", "C", "G", "T"),
                                                  c("A", "C", "G", "T"))
                                  )
                     )
  matrixFile <- paste(Sys.getpid(), "-lastMatrix.dat", sep="") 
  write.table(lastzMatrix[[distance]], file=matrixFile, quote=FALSE,
              sep=" ", row.names=TRUE, col.names=TRUE)
  ## -a: Gap existence cost.
  ## -b: Gap extension cost.
  ## -e: Minimum alignment score.
  ## -p: Specify a match/mismatch score matrix.  Options -r and -q will be ignored.
  ## -s: Specify which query strand should be used: 0 means reverse only, 1 means forward only, and 2 means both.
  lastOptiosn <- list(near=paste("-a 600 -b 150 -e 3000 -p", matrixFile, "-s 2"),
                      medium=paste("-a 400 -b 30 -e 4500 -p", matrixFile, "-s 2"),
                      far=paste("-a 400 -b 30 -e 6000 -p", matrixFile, "-s 2")
                      )
  formatMapping <- list(MAF=1, tabular=0)
  message("last")
  mc.cores <- as.integer(mc.cores)
  #if(mc.cores != 1L){
  #  stop("The computation in parallel is still not working properly")
  #}
  if(mc.cores == 1L){
    cmd <- paste("lastal", lastOptiosn[[distance]],
               "-f", formatMapping[[format]],
               db, queryFn, ">", outputFn)
  }else{
    cmd <- paste("parallel-fasta", "-j", mc.cores,
                 "\"lastal", lastOptiosn[[distance]],
                 "-f", formatMapping[[format]], db, "\"", "<", queryFn,
                 ">", outputFn)
  }
  if(echoCommand){
    message(cmd)
  }else{
    my.system(cmd)
  }
  unlink(matrixFile)
  return("success")
}



