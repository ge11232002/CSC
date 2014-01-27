\name{axt-class}
\Rdversion{1.1}
\docType{class}
\alias{axt}
\alias{axt-class}
\alias{[,axt,ANY,ANY-method}
\alias{c,axt-method}
\alias{length,axt-method}
\alias{queryRanges}
\alias{queryRanges,axt-method}
\alias{querySeqs}
\alias{querySeqs,axt-method}
\alias{score,axt-method}
\alias{subAxt}
\alias{subAxt,axt-method}
\alias{symCount}
\alias{symCount,axt-method}
\alias{targetRanges}
\alias{targetRanges,axt-method}
\alias{targetSeqs}
\alias{targetSeqs,axt-method}

\title{Class \code{"axt"}}
\description{
  The axt S4 object
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("axt", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{targetRanges}:}{Object of class \code{"GRanges"}: The ranges of net alignments on reference genome.}
    \item{\code{targetSeqs}:}{Object of class \code{"DNAStringSet"}: The alignment sequences of reference genome.}
    \item{\code{queryRanges}:}{Object of class \code{"GRanges"}: The ranges of net alignments on query genome.}
    \item{\code{querySeqs}:}{Object of class \code{"DNAStringSet"}: The alignment sequences of query genome.}
    \item{\code{score}:}{Object of class \code{"integer"}: The alignment score.}
    \item{\code{symCount}:}{Object of class \code{"integer"}: The alignment length.}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "axt", i = "ANY", j = "ANY")}: axt getter}
    \item{c}{\code{signature(x = "axt")}: axt concatenator.}
    \item{length}{\code{signature(x = "axt")}: Get the number of alignments.}
    \item{queryRanges}{\code{signature(x = "axt")}: Get the ranges of query genome.}
    \item{querySeqs}{\code{signature(x = "axt")}: Get the alignment sequences of query genome.}
    \item{score}{\code{signature(x = "axt")}: Get the alignment score.}
    \item{subAxt}{\code{signature(x = "axt")}: Get subset of axt alignments based on range.}
    \item{symCount}{\code{signature(x = "axt")}: Get the alignment lengths.}
    \item{targetRanges}{\code{signature(x = "axt")}: Get the ranges of reference genome.}
    \item{targetSeqs}{\code{signature(x = "axt")}: Get the alignment sequences of reference genome.}
	 }
}
\author{
  Ge Tan
}

\seealso{
  \code{\link{readAxt}}
}
\examples{
showClass("axt")
}
\keyword{classes}
