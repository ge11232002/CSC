\name{readAxt}
\alias{readAxt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
readAxt function.
}
\description{
This function reads the axt files into a \code{axt} object.
}
\usage{
readAxt(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Object of \code{character}. The length can be one or more.
}
}
\details{
  This function reads the axt files. 
  The coordinates in \code{axt} object is 1-based.
}
\value{
  A object \code{axt} is returned.
}

\author{
Ge Tan
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{axt}}
}
\examples{
  \dontrun{
    axtFileshg19danRer7 = list.files(path="/Users/gtan/CSC/CNEr/axtNet", 
      pattern=".*hg19\\.danRer7\\.*", full.names=TRUE)
    axtshg19danRer7 = readAxt(axtFileshg19danRer7)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
