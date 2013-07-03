\name{openDb-methods}
\docType{methods}
\alias{openDb-methods}
\alias{openDb,JASPARDb-method}
\title{Open the SQLite connection of the object \code{JASPARDb}}.
\description{
Open the SQLite connection of the object \code{JASPARDb}.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "JASPARDb")}}{
%%  ~~describe this method here~~
}
}}
\seealso{
  \code{\link{closeDb}}
}
\examples{
  jasparDb = JASPARDb(release_name="JASPAR_2010", metadata_dirpath="JASPAR_2010.sqlite")
  closeDb(jasparDb)
  openDb(jasparDb)
}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
