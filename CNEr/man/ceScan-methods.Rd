\name{ceScan-methods}
\docType{methods}
\alias{ceScan-methods}
\alias{ceScan}
\alias{ceScan,axt,GRanges,GRanges,Seqinfo-method}
\alias{ceScan,axt,GRanges,missing,missing-method}
\alias{ceScan,axt,missing,GRanges,Seqinfo-method}
\alias{ceScan,axt,missing,missing,missing-method}
\alias{ceScan,character,character,character,Seqinfo-method}
\alias{ceScan,character,character,missing,missing-method}
\alias{ceScan,character,missing,character,Seqinfo-method}
\alias{ceScan,character,missing,missing,missing-method}
\title{ceScan function}
\description{
  This is the main function for conserved noncoding elements (CNEs) identification.
}

\usage{
  ceScan(axts, tFilter, qFilter, qSizes, thresholds="49,50")
}

\arguments{
  \item{axts}{
    A \code{axt} object or \code{character} object with the paths of axt files.
  }
  \item{tFilter}{
    A \code{GRanges} object or \code{character} object with the path of bed file for target genome filter. 
    This argument can also be missing when target filter is not available.
  }
  \item{qFilter}{
    A \code{GRanges} object or \code{character} object with the path of bed file for
    query genome filter.
    This argument can also be missing when query filter is not available.
  }
  \item{qSizes}{
    A \code{Seqinfo} object which contains the \code{seqnames} and \code{seqlengths} for query genome.
    This argument can be missing when \code{qFilter} is missing.
  }
  \item{thresholds}{
    A \code{character} object specifiying the scanning windows and minimal score.
    It can be specified in th form of "45,50" with scanning windows 50 and minial score 45.
    More than one thresholds can be provided.
  }
}
\section{Methods}{
\describe{

\item{\code{signature(axts = "axt", tFilter = "GRanges", qFilter = "GRanges", qSizes = "Seqinfo")}}{

}

\item{\code{signature(axts = "axt", tFilter = "GRanges", qFilter = "missing", qSizes = "missing")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(axts = "axt", tFilter = "missing", qFilter = "GRanges", qSizes = "Seqinfo")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(axts = "axt", tFilter = "missing", qFilter = "missing", qSizes = "missing")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(axts = "character", tFilter = "character", qFilter = "character", qSizes = "Seqinfo")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(axts = "character", tFilter = "character", qFilter = "missing", qSizes = "missing")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(axts = "character", tFilter = "missing", qFilter = "character", qSizes = "Seqinfo")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(axts = "character", tFilter = "missing", qFilter = "missing", qSizes = "missing")}}{
%%  ~~describe this method here~~
}
}}

\details{
  ceScan scan the axts alignmnets and identify the CNEs.
  ceScan can accept axts in \code{axt} object and filter in \code{GRanges} object,
  or directly the axt files and bed files.
  When the axt files and bed files are ready for computation, 
  it is recommended to use them directly rather than read them into R first.
  
  The details of algorithm will given in the vignette.
}

\value{
  A \code{list} of \code{data.frame} is returned. 
  Each element of the list is for one threshold.
}

\author{
  Ge Tan
}

\examples{
  \dontrun{
    axtFileshg19danRer7 = list.files(path="/Users/gtan/CSC/CNEr/axtNet", 
      pattern=".*hg19\\.danRer7\\.*", full.names=TRUE)
    axtshg19danRer7 = readAxt(axtFileshg19danRer7)
    axtFilesdanRer7hg19 = list.files(path="/Users/gtan/CSC/CNEr/axtNet", 
      pattern=".*danRer7\\.hg19\\.*", full.names=TRUE)
    axtsdanRer7hg19 = readAxt(axtFilesdanRer7hg19)
    bedhg19 = readBed("/Users/gtan/CSC/CNEr/filters/filter_regions.hg19.bed")
    beddanRer7 = readBed("/Users/gtan/CSC/CNEr/filters/filter_regions.danRer7.bed")
    library(rtracklayer)
    qSizesdanRer7 = seqinfo(TwoBitFile("/Users/gtan/CSC/CNEr/2bit/danRer7.2bit"))
    qSizeshg19 = seqinfo(TwoBitFile("/Users/gtan/CSC/CNEr/2bit/hg19.2bit"))
  }
}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
