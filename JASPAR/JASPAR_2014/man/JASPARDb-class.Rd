\name{JASPARDb-class}
\Rdversion{1.1}
\docType{class}
\alias{JASPARDb-class}
\alias{searchDb,JASPARDb-method}

\title{JASPARDb object}
\description{
  The JASPARDb object class is a container for the information about one JASPAR database.
}
\usage{
## Constructor
  JASPARDb(provider=character(), provider_version=character(),
           release_date=character(), release_name=character(),
           source_url=character(), db_dirpath=character())
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("JASPARDb", ...)}.
%%  ~~ describe objects here ~~ 
}
\arguments{
  \item{provider}{
    A character vector: the provider of JASPAR database (e.g. SQLite).
  }
  \item{provider_version}{
    A character vector: the version of provider.
  }
  \item{release_date}{
    A character vector: the release date of the database.
  }
  \item{release_name}{
    A character vector: the release name of the database.
  }
  \item{source_url}{
    A character vector: the url of origin.
  }
  \item{db_dirpath}{
    A character vector: the path of the database file.
  }
}

\section{Methods}{
  In the code snippets below, \code{x} is an JASPARDb object.
  \describe{
    \item{}{
      \code{searchDb(x, ID=NULL, name=NULL, species=NULL, class=NULL, type=NULL}:
      Search the JASPAR database based on the criteria of ID, name, species, class, type. These criterias have the "ALL" relationship. If none of them is set, all records will be returned.
    }
    \item{}{
      \code{openDb(x)}:
      Open the SQLite connection of the object \code{JASPARDb}.
    }
    \item{}{
      \code{closeDb(x)}:
      Close the SQLite connection of the object \code{JASPARDb}.
    }
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
\link{JASPAR}
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
  searchDb(JASPAR_2010, class="Zipper-Type", type="SELEX")
  openDb(JASPAR_2010)
  closeDb(JASPAR_2010)
}
\keyword{classes}
