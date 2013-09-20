\name{XMatrixList}
\Rdversion{1.1}
\docType{class}
\alias{XMatrixList}
\alias{bg,XMatrixList-method}
\alias{ID,XMatrixList-method}
\alias{Matrix,XMatrixList-method}
\alias{matrixClass,XMatrixList-method}
\alias{matrixType,XMatrixList-method}
\alias{name,XMatrixList-method}
\alias{strand,XMatrixList-method}

\title{Class \code{"XMatrixList"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
  The XMatrixList virtual class is a container for storing a collection of XMatrix objects.
  No object can be constructed directly from this virtual and it has three subclasses: PFMatrixList, PWMatrixList and ICMatrixList.
  Basically it is a SimpleList and is designed for manipulating the set of XMatrix objects as a whole.
}

\usage{
## Constructors:
PFMatrixList(..., use.names=TRUE)
PWMatrixList(..., use.names=TRUE)
ICMatrixList(..., use.names=TRUE)

## Accessor-like methods:
\S4method{ID}{XMatrixList}(x)
\S4method{name}{XMatrixList}(x)

}

\arguments{
  \item{...}{
    The XMatrix objects are supplied in ....
  }
  \item{use.names}{
    A logical value. When TRUE, the names of the XMatrix will be kept.
  }
}

\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
Ge Tan
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{XMatrix}},
  \code{\link{searchSeq}},
  \code{\link{searchAln}}
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}

\examples{
  \dontrun{
    PFMatrixList(pfm1=pfm, pfm2=pfm, use.names=TRUE)
  }
}
\keyword{classes}

