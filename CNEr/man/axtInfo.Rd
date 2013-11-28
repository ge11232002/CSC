\name{axtInfo}
\alias{axtInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  axtInfo function
}
\description{
  Given the path of axt file, retrieve the alignments' withs information.
}
\usage{
axtInfo(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Object of \code{character}. The length can be one or more.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  A vector of \code{integer} is returned. 
  It stores the withds of all the alignments.
}

\author{
  Ge Tan
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{link{readAxt}}
}
\examples{
  \dontrun{
    axtFileshg19danRer7 = list.files(path="/Users/gtan/CSC/CNEr/axtNet", 
      pattern=".*hg19\\.danRer7\\.*", full.names=TRUE)
    ans = axtInfo(axtFileshg19danRer7)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
