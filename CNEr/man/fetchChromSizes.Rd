\name{fetchChromSizes}
\alias{fetchChromSizes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fetchChromSizes function.
}
\description{
This function tries to automate the fetch of chrom sizes for assembly from UCSC and other sources.
}
\usage{
fetchChromSizes(assembly)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{assembly}{
%%     ~~Describe \code{x} here~~
  A \code{character} object: the canonical name of assembly, i.e., hg19 for UCSC.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function utilises mysql query for UCSC assemblies.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A object of \code{Seqinfo} is returned.
}
\author{
%%  ~~who you are~~
Ge Tan
}
\note{
%%  ~~further notes~~
Currently the assemblies from UCSC are supported.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
  fetchChromSizes("hg19")
}

