\name{JASPARDb-class}
\Rdversion{1.1}
\docType{class}
\alias{JASPARDb-class}
\alias{searchDb,JASPARDb-method}

\title{Class \code{"JASPARDb"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("JASPARDb", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{provider}:}{Object of class \code{"character"} ~~ }
    \item{\code{provider_version}:}{Object of class \code{"character"} ~~ }
    \item{\code{release_date}:}{Object of class \code{"character"} ~~ }
    \item{\code{release_name}:}{Object of class \code{"character"} ~~ }
    \item{\code{source_url}:}{Object of class \code{"character"} ~~ }
    \item{\code{metadata_dirpath}:}{Object of class \code{"character"} ~~ }
    \item{\code{conn}:}{Object of class \code{"SQLiteConnection"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{searchDb}{\code{signature(x = "JASPARDb")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("JASPARDb")
}
\keyword{classes}
