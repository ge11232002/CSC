\name{reverseCigar}
\alias{reverseCigar}

\title{reverse the cigar string}

\description{reverse the cigar string}

\usage{
  reverseCigar(cigar, ops=CIGAR_OPS)
}

\arguments{
  \item{cigar}{
    A character vector of cigar strings.
  }
  \item{ops}{
    A character vector of the extended CIGAR operations. 
    By default, CIGAR_OPS is used.
  }
}

\value{
  A character vector contains the revered cigar strings.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link[GenomicRanges]{cigar-utils}}
}

\examples{
  cigar = c("20M15I10D", "10D15I20M")
  reverseCigar(cigar)
}
