\name{closeDb-methods}
\docType{methods}
\alias{closeDb-methods}
\alias{closeDb,JASPARDb-method}
\title{Close the SQLite connection of the object \code{JASPARDb}}
\description{
Close the SQLite connection of the object \code{JASPARDb}. 
}
\section{Methods}{
\describe{

\item{\code{signature(x = "JASPARDb")}}{
%%  ~~describe this method here~~
}
}}
\seealso{
  \code{\link{openDb}}
}
\examples{
  jasparDb = JASPARDb(release_name="JASPAR_2010", metadata_dirpath="JASPAR_2010.sqlite")
  closeDb(jasparDb)
}
\keyword{methods}
