\name{XMatrixList}
\Rdversion{1.1}
\docType{class}
\alias{XMatrixList}
\alias{PFMatrixList}
\alias{bg}
\alias{ID}
\alias{Matrix}
\alias{matrixClass}
\alias{matrixType}
\alias{name}
\alias{strand}

\title{Class \code{"XMatrixList"}}
\description{
  The XMatrixList virtual class is a container for storing a collection of XMatrix objects.
  No object can be constructed directly from this virtual and it has three subclasses: PFMatrixList, PWMatrixList and ICMatrixList.
  Basically it is a SimpleList and is designed for manipulating the set of XMatrix objects as a whole.
}

\usage{
## Constructors:
PFMatrixList(..., use.names=TRUE)

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

\author{
Ge Tan
}


\seealso{
  \code{\link{XMatrix}},
}

\examples{
  \dontrun{
    PFMatrixList(pfm1=pfm, pfm2=pfm, use.names=TRUE)
  }
}
\keyword{classes}

