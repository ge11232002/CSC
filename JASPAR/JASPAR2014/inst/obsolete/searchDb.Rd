\name{searchDb}
\alias{searchDb}
\title{searchDb method}
\description{
  The searchDb method is the major function to query the JASPAR SQLite database.
}
\usage{
  searchDb(x, ID=NULL, name=NULL, species=NULL, class=NULL, type=NULL)
}

\arguments{
  \item{x}{
    A character vector of length 1: The path of SQLite database.
  }
  \item{ID}{
    A character vector: the TFBS IDs, e.g. "MA0004"
  }
  \item{name}{
    A character vector: the TFBS names, e.g. "Arnt"
  }
  \item{species}{
    A character vector: the species, e.g. "Mus musculus"
  }
  \item{class}{
    A character vector: the classes, e.g. "Zipper-Type"
  }
  \item{type}{
    A character vector: the types, e.g. "SELEX"
  }
}

\section{Methods}{
  Search the JASPAR database based on the criteria of ID, name, species, class, type. All the arguments can be NULL except the x. If so, all the records will be retrieved. All the criterias have the "AND" relationship and inside each criteria, all the elements have the "OR" relationship.
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
\link{JASPAR}
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
  searchDb(JASPAR2014, class="Zipper-Type", type="SELEX")
}
\keyword{classes}
