\name{XMatrix-class}
\Rdversion{1.1}
\docType{class}
\alias{XMatrix-class}
\alias{bg,XMatrix-method}
\alias{bg<-,XMatrix-method}
\alias{ID,XMatrix-method}
\alias{ID<-,XMatrix-method}
\alias{length,XMatrix-method}
\alias{Matrix,XMatrix-method}
\alias{Matrix<-,XMatrix-method}
\alias{matrixClass,XMatrix-method}
\alias{matrixClass<-,XMatrix-method}
\alias{name,XMatrix-method}
\alias{name<-,XMatrix-method}
\alias{revcom,XMatrix-method}
\alias{show,XMatrix-method}
\alias{strand,XMatrix-method}
\alias{strand<-,XMatrix-method}
\alias{tags,XMatrix-method}
\alias{update,XMatrix-method}

\title{Class \code{"XMatrix"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
XMatrix is a virtual class. No objects can be created from it directly. Three classes are derived from this class: PFMatrix, PWMatrix and ICMatrix.

PFMatrix is a class whose instances are objects representing raw position frequency matrices (PFMs).

PWMatrix is a class whose instances are objects representing position weight matrices (PWMs).

ICMatrix is a class whose instances are objects representing information content matrices (ICMs).
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Methods}{
  \describe{
    \item{bg}{\code{signature(x = "XMatrix")}: ... }
    \item{bg<-}{\code{signature(x = "XMatrix")}: ... }
    \item{ID}{\code{signature(x = "XMatrix")}: ... }
    \item{ID<-}{\code{signature(x = "XMatrix")}: ... }
    \item{length}{\code{signature(x = "XMatrix")}: ... }
    \item{Matrix}{\code{signature(x = "XMatrix")}: ... }
    \item{Matrix<-}{\code{signature(x = "XMatrix")}: ... }
    \item{matrixClass}{\code{signature(x = "XMatrix")}: ... }
    \item{matrixClass<-}{\code{signature(x = "XMatrix")}: ... }
    \item{name}{\code{signature(x = "XMatrix")}: ... }
    \item{name<-}{\code{signature(x = "XMatrix")}: ... }
    \item{revcom}{\code{signature(x = "XMatrix")}: ... }
    \item{show}{\code{signature(object = "XMatrix")}: ... }
    \item{strand}{\code{signature(x = "XMatrix")}: ... }
    \item{strand<-}{\code{signature(x = "XMatrix")}: ... }
    \item{tags}{\code{signature(x = "XMatrix")}: ... }
    \item{update}{\code{signature(object = "XMatrix")}: ... }
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
showClass("XMatrix")
}
\keyword{classes}
