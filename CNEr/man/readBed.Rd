\name{readBed}
\alias{readBed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
readBed function
}
\description{
readBed reads bed file with an embeded C IO function.
}
\usage{
readBed(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  The path of bed file.
}
}
\details{
  This function is designed to read the bed file only with three mandatory columns, 
  i.e., "chrom", "chromStart", "chromEnd".
  This function utilises the C interface for speedy import.
  For general bed file import, please use the \code{import.bed} from package \code{rtracklayer}.
  
  In bed file, the "chromStart" is on 0-based coordinate while "chromEnd" is on 1-based coordinate.
  For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
  When it is read into GRanges, both the chromStart and chromEnd are on 1-based coordinate, i.e.,
  chromStart=1 and chromEnd=100.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  A \code{GRanges} is returned. 
  When no strand information is availabe in bed file,
  all the ranges are assumed to be on the positive strand.
}


\author{
Ge Tan
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{import.bed}}
}
\examples{
  \dontrun{
    bedhg19 = readBed("/Users/gtan/CSC/CNEr/filters/filter_regions.hg19.bed")
    beddanRer7 = readBed("/Users/gtan/CSC/CNEr/filters/filter_regions.danRer7.bed")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
