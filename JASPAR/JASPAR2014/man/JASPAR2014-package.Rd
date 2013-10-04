\name{JASPAR2014-package}
\alias{JASPAR2014-package}
\docType{package}
\title{
Data package for JASPAR 2014
}
\description{
Data package for JASPAR 2014.
}
\details{
\tabular{ll}{
Package: \tab JASPAR2014\cr
Version: \tab 0.99.1\cr
Date: \tab 2013-10-04\cr
Depends: \tab R (>= 2.9.0), Biostrings(>= 2.26.3)\cr
Suggests: \tab TFBSTools\cr
License: \tab GPL-2 \cr
URL: \tab http://jaspardev.genereg.net/\cr
Type: \tab Package\cr
NeedsCompilation: \tab no\cr
LazyData: \tab yes\cr
}

}

\author{
Ge Tan <ge.tan09@imperial.ac.uk> 

Maintainer: Ge Tan <ge.tan09@imperial.ac.uk>
}
\references{
  See \url{http://jaspardev.genereg.net/} for more details about JASPAR.
}
\keyword{ package }
\examples{
## load the library
library(JASPAR2014)
## list the contents that are loaded into memory
ls("package:JASPAR2014")
## show the db object that is loaded by calling it's name
JASPAR2014
}



