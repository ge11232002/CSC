\name{JASPAR-class}
\Rdversion{1.1}
\docType{class}
\alias{JASPAR-class}
\alias{[,JASPAR-method}
\alias{FMatrix,JASPAR-method}
\alias{TFBSinfo,JASPAR-method}
\alias{names,JASPAR-method}
\alias{seqs,JASPAR-method}

\title{JASPAR object}
\description{
  The JASPAR object class is a container for storing a set of TFBS records which can be fetched from JASPAR website.
}
\usage{
## Constructors:
JASPAR(ID=character(), collection=Rle(), version=Rle(),
       name=character(), species=Rle(), TF_class=Rle(),
       medline=character(), family=Rle(), tax_group=Rle(),
       acc=character(), type=Rle(), pazar_tf_id=character(),
       comment=character(), FMatrix=list(), seqs=list())
}
\arguments{
  \item{ID}{
    A character vector: a identifier for each model, without version postfix.
  }
  \item{collection}{
    A Rle object: the name of collection.
  }
  \item{version}{
    A Rle object: version information.
  }
  \item{name}{
    A character vector: the name of model.
  }
  \item{species}{
    A Rle object: the species source for the sequences, in Latin. NA is allowed.
  }
  \item{TF_class}{
    A Rle object: structural class of the transcription factor, based on the TFCaT system.
  }
  \item{medline}{
    A character vector: A link to the relevant publication reporting the sites used in the mode building.
  }
  \item{family}{
    A Rle object: structural sub-class of the transcription factor, based on the TFCaT system.
  }
  \item{tax_group}{
    A Rle object: group of species, currently consisting of 4 larger groups: vertebrate, insect, plant, chordate.
  }
  \item{acc}{
    A character vector: a representative protein accession number in Genbank for the transcription factor. Human takes precedence if several exists. 
  }
  \item{type}{
    A Rle object: methodology used for matrix construction.
  }
  \item{pazar_tf_id}{
    A character vector: a link to the PAZAR database.
  }
  \item{comment}{
    A character vector: for some matrices, a curator comment is added.
  }
  \item{FMatrix}{
    A list of matrix: count matrices.
  }
  \item{seqs}{
    A list of DNAStringSet: the binding sites sequences.
  }

}
\section{Methods}{
  In the code snippets below, \code{x} is an JASPAR object.
  \describe{
    \item{}{
      \code{FMatrix(x)}:
      Get the count matrices from a JASPAR object.
    }
    \item{}{
      \code{TFBSinfo(x)}:
      Get the fundamental matrix information.
    }
    \item{}{
      \code{seqs(x)}:
      Get the Binding Sites sequences as a list of DNAStringSet.
    }
    \item{}{
      \code{x[1..n]}:
      Get the subset of JASPAR object.
    }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
  Ge Tan
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \link{JASPARDb},
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
  ans = searchDb(JASPAR2014, class="Zipper-Type", type="SELEX")
  FMatrix(ans)
  TFBSinfo(ans)
  seqs(ans)
  ans[1:5]
}
\keyword{classes}
