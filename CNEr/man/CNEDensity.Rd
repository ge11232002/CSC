\name{CNEDensity}
\alias{CNEDensity}

\title{
CNEDensity function
}
\description{
This function queries the database and generates the CNEs density values.
}
\usage{
CNEDensity(dbName, tableName, whichAssembly=c("L","R"),
           chr, start, end, windowSize,
           minLength=NULL)
}

\arguments{
  \item{dbName}{
    A object of \code{character}, the path of the local SQLite database.
  }
  \item{tableName}{
    A object of \code{character}, the name of table for this CNE data table.
  }
  \item{whichAssembly}{
    A object of \code{character}, the genome to fetch is in the "Left" columns or "Right" columns of the table.
  }
  \item{chr}{
     A object of \code{character}, the chromosome to query.
  }
  \item{start, end}{
    A object of \code{integer}, the start and end coordiante to fetch the CNEs.
  }
  \item{windowSize}{
    A object of \code{integer}, the window size in kb used to smooth the CNEs.
  }
  \item{minLength}{
    A object of \code{integer}, the minimal length of CNEs to fetch.
  }
}

\value{
  A \code{matrix} is returned. 
  The first column is the coordinates and the second column is the density values.
}

\author{
Ge Tan
}


\examples{
    dbName = file.path(system.file("extdata", package=CNEr),
                       "cne.sqlite")
    windowSize= 300L
    minLength = 50L
    cneHg19DanRer7_45_50 = CNEDensity(dbName=dbName, 
                                  tableName="cne_twoWay_danRer7_hg19_len50_id900_v1", 
                                  whichAssembly="R", chr=chr, start=start,
                                  end=end, windowSize=windowSize, 
                                  minLength=minLength)
    cneHg19DanRer7_48_50 = CNEDensity(dbName=dbName, 
                                  tableName="cne_twoWay_danRer7_hg19_len50_id960_v1", 
                                  whichAssembly="R", chr=chr, start=start,
                                  end=end, windowSize=windowSize, 
                                  minLength=minLength)
    cneHg19DanRer7_49_50 = CNEDensity(dbName=dbName, 
                                  tableName="cne_twoWay_danRer7_hg19_len50_id980_v1", 
                                  whichAssembly="R", chr=chr, start=start,
                                  end=end, windowSize=windowSize, 
                                  minLength=minLength)
}


