

readAxt_obsolete = function(axtFile){
  axtContent = readLines(axtFile)
  axtContent = axtContent[grep("^#", axtContent, invert=TRUE)]
  axtContent = axtContent[axtContent != ""]
  axtContent = split(axtContent, c("header", "target", "query"))
  header = strsplit(axtContent$header, split=" ", fixed=TRUE)
  myAxt = axt(targetNames=Rle(sapply(header, "[",2)),
              targetRanges=IRanges(start=as.integer(sapply(header, "[",3)), end=as.integer(sapply(header, "[",4)), names=sapply(header, "[",1)),
              targetSeqs=DNAStringSet(axtContent$target),
              queryNames=Rle(sapply(header, "[",5)),
              queryRanges=IRanges(start=as.integer(sapply(header, "[",6)), end=as.integer(sapply(header, "[",7)), names=sapply(header, "[",1)),
              querySeqs=DNAStringSet(axtContent$query),
              strand=Rle(sapply(header, "[",8)),
              score=Rle(as.integer(sapply(header, "[",9)))
              )
  myAxt2 = axt2(targetRanges=GRanges(seqnames=Rle(sapply(header, "[",2)),
                                     ranges=IRanges(start=as.integer(sapply(header, "[",3)), end=as.integer(sapply(header, "[",4)), names=sapply(header, "[",1)),
                                     strand=Rle("+")),
                targetSeqs=DNAStringSet(axtContent$target),
                queryRanges=GRanges(seqnames=Rle(sapply(header, "[",5)),
                                    ranges=IRanges(start=as.integer(sapply(header, "[",6)), end=as.integer(sapply(header, "[",7)), names=sapply(header, "[",1)),
                                    strand=Rle(sapply(header, "[",8))),
                querySeqs=DNAStringSet(axtContent$query),
                score=as.integer(sapply(header, "[",9))
                )
  rm(axtContent)
  rm(header)
  gc()
  return(myAxt)
}

# tFilterFile = "/export/data/CNEs/danRer7/filters/filter_regions.danRer7.bed"
# qFilterFile = "/export/data/CNEs/hg19/filters/filter_regions.hg19.bed"
readBed = function(bedFile){
## This GRanges is in 1-based.
  require(rtracklayer)
  bed = import(bedFile, asRangedData = FALSE)
  strand(bed) = "+"
  bed = reduce(bed)
}
# system.time(foo<-readBed(bedFile))
#    user  system elapsed
#    51.890   0.757  52.659

#############################C version######################################
readBedToGRanges = function(bedFile=NULL){
  require(GenomicRanges)
  if(is.null(bedFile)){
    return(NULL)
  }
  if(!file.exists(bedFile)){
    stop("No such file ", bedFile) 
  }
  dyn.load("~/Repos/CSC/CNEr/src/io.so")
  bed = .Call("myReadBed", bedFile)
  bed = GRanges(seqnames=Rle(bed[[1]]),
                ranges=IRanges(start=bed[[2]], end=bed[[3]]),
                strand=factor("+"))
  return(bed)
}
# bedHuman = readBedToGRanges("/export/data/CNEs/hg19/filters/filter_regions.hg19.bed")
# bedZebrafish = readBedToGRanges("/export/data/CNEs/danRer7/filters/filter_regions.danRer7.bed")
# system.time(foo<-readBedToGRanges(bedFile))
#    user  system elapsed
#    2.272   0.133   2.414

# axtFiles = c("/export/downloads/ucsc/axtNet/hg19/chr2.hg19.danRer7.net.axt.gz", "/export/downloads/ucsc/axtNet/hg19/chr3.hg19.danRer7.net.axt.gz", "/export/downloads/ucsc/axtNet/hg19/chr4.hg19.danRer7.net.axt.gz")
# axtFiles = list.files(path="/export/downloads/ucsc/axtNet/hg19", pattern=".*hg19\\.mm10*", full.names=TRUE)
# axt = readAxt(axtFiles)
readAxt = function(axtFiles){
  # Read axt files into R axt object
  index_noexists = !file.exists(axtFiles)
  if(any(index_noexists)){
    stop("No such file ", paste(axtFiles[index_noexists], sep=" "))
  }
  require(GenomicRanges)
  require(Biostrings)
  dyn.load("~/Repos/CSC/CNEr/src/io.so")
  myAxt = .Call("myReadAxt", axtFiles)
  axt = axt(targetRanges=GRanges(seqnames=Rle(myAxt[[1]]),
                                 ranges=IRanges(start=myAxt[[2]],
                                                end=myAxt[[3]]),
                                 strand=Rle(myAxt[[4]])),
            targetSeqs=DNAStringSet(myAxt[[5]]),
            queryRanges=GRanges(seqnames=Rle(myAxt[[6]]),
                                ranges=IRanges(start=myAxt[[7]],
                                               end=myAxt[[8]]),
                                strand=Rle(myAxt[[9]])),
            querySeqs=DNAStringSet(myAxt[[10]]),
            score=myAxt[[11]],
            symCount=myAxt[[12]]
            )
  return(axt)
}
##-----------------------------------IO from Biostrings----------------------------
.normargInputFilepath <- function(filepath)
{
    if (!is.character(filepath) || any(is.na(filepath)))
        stop("'filepath' must be a character vector with no NAs")
    ## First pass: expand local paths and download any remote file.
    filepath2 <- character(length(filepath))
    for (i in seq_len(length(filepath))) {
        fp <- filepath[i]
        con <- file(fp)
        con_class <- class(con)[1L]
        close(con)
        if (con_class == "url") {
            filepath2[i] <- tempfile()
            download.file(fp, filepath2[i])
        } else {
            filepath2[i] <- path.expand(fp)
        }
    }
    ## Second pass: check the type of the local files (all files are
    ## now local).
    filetype <- character(length(filepath2))
    for (i in seq_len(length(filepath2))) {
        fp <- filepath2[i]
        con <- file(fp)
        ## Ugly trick to get the type of 'con'. Is there a better way?
        filetype[i] <- showConnections(TRUE)[as.character(con), "class"]
        close(con)
        if (filetype[i] != "file")
            stop("file \"", filepath[i], "\" ",
                 "has unsupported type: ", filetype[i])
    }
    names(filepath2) <- filetype
    filepath2
}

.ExternalFilePtr.close <- function(x)
    .Call("ExternalFilePtr_close", x)

### Returns a list of "external file pointers".
.openInputFiles <- function(filepath)
{
    filepath2 <- .normargInputFilepath(filepath)
    ans <- lapply(filepath2,
           function(fp)
           {
               efp <- .Call("new_input_ExternalFilePtr", fp)
               reg.finalizer(efp, .ExternalFilePtr.close, onexit=TRUE)
               efp
           })
    names(ans) <- filepath
    ans
}

### 'efp_list' must be a list of "external file pointers" returned by
### .openInputFiles() or .openOutputFiles().
.closeFiles <- function(efp_list)
{
    for (efp in efp_list) .ExternalFilePtr.close(efp)
}

.normargNrec <- function(nrec)
{
    if (!isSingleNumber(nrec))
        stop("'nrec' must be a single integer value")
    if (!is.integer(nrec))
        nrec <- as.integer(nrec)
    nrec
}

.normargSkip <- function(skip)
{
    if (!isSingleNumber(skip))
        stop("'skip' must be a single integer value")
    if (!is.integer(skip))
        skip <- as.integer(skip)
    if (skip < 0L)
        stop("'skip' cannot be negative")
    skip
}

axt.info = function(filepath, nrec=-1L, skip=0L, use.names=TRUE, seqtype="B"){
  efp_list <- .openInputFiles(filepath)
  on.exit(.closeFiles(efp_list))
  nrec <- .normargNrec(nrec)
  skip <- .normargSkip(skip)
  use.names <- normargUseNames(use.names)
  seqtype <- match.arg(seqtype, c("B", "DNA", "RNA", "AA"))
  lkup <- get_seqtype_conversion_lookup("B", seqtype)
  .Call("axt_info", efp_list, nrec, skip, use.names, lkup)
}

