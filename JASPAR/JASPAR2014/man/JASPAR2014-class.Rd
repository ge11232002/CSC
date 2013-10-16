\name{JASPAR2014-class}
\Rdversion{1.1}
\docType{class}
\alias{JASPAR2014-class}
\alias{JASPAR2014}

\title{JASPAR2014 object}
\description{
  The JASPAR2014 object class is a thin class for storing the path of JASPAR2014 style SQLite file. 
}

\section{Slots}{
  \describe{
    \item{\code{db}:}{Object of class \code{"character"}: a character string of the path of SQLite file.}
  }
}

\author{
  Ge Tan
}

\seealso{
  \link{JASPAR2014SitesSeqs},
}
\examples{
  \dontrun{
    library(JASPAR2014)
    JASPAR2014
  }
}
\keyword{classes}
