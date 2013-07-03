\name{JASPAR-class}
\Rdversion{1.1}
\docType{class}
\alias{JASPAR-class}
\alias{[,JASPAR-method}
\alias{FMatrix,JASPAR-method}
\alias{metadata,JASPAR-method}
\alias{names,JASPAR-method}
\alias{seqs,JASPAR-method}
\alias{show,JASPAR-method}

\title{Class \code{"JASPAR"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("JASPAR", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{ID}:}{Object of class \code{"character"} ~~ }
    \item{\code{collection}:}{Object of class \code{"Rle"} ~~ }
    \item{\code{version}:}{Object of class \code{"Rle"} ~~ }
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
    \item{\code{species}:}{Object of class \code{"Rle"} ~~ }
    \item{\code{TF_class}:}{Object of class \code{"Rle"} ~~ }
    \item{\code{medline}:}{Object of class \code{"character"} ~~ }
    \item{\code{family}:}{Object of class \code{"Rle"} ~~ }
    \item{\code{tax_group}:}{Object of class \code{"Rle"} ~~ }
    \item{\code{acc}:}{Object of class \code{"character"} ~~ }
    \item{\code{type}:}{Object of class \code{"Rle"} ~~ }
    \item{\code{pazar_tf_id}:}{Object of class \code{"character"} ~~ }
    \item{\code{comment}:}{Object of class \code{"character"} ~~ }
    \item{\code{FMatrix}:}{Object of class \code{"list"} ~~ }
    \item{\code{seqs}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "JASPAR")}: ... }
    \item{FMatrix}{\code{signature(x = "JASPAR")}: ... }
    \item{metadata}{\code{signature(x = "JASPAR")}: ... }
    \item{names}{\code{signature(x = "JASPAR")}: ... }
    \item{seqs}{\code{signature(x = "JASPAR")}: ... }
    \item{show}{\code{signature(object = "JASPAR")}: ... }
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
showClass("JASPAR")
}
\keyword{classes}
