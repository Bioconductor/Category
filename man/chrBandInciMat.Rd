\name{chrBandInciMat}

\alias{chrBandInciMat}
\alias{makeChrBandInciMat}

\title{Create an incidence matrix mapping chromosome bands to genes}

\description{

  These functions return a 0/1 incidence matrix with a row for each
  gene ID with at least one chromosme band annotation and a column for
  each chromosome band with at least one gene annotation.
  
}

\usage{
chrBandInciMat(chip, univ = NULL)
makeChrBandInciMat(chrGraph)
}

\arguments{
  \item{chip}{A string giving the annotation source.  For example, \code{"hgu133plus2"}}
  \item{univ}{A vector of gene IDs (these should be Entrez IDs for
  most annotation sources).  The the annotations will be limited
  to those in the set specified by \code{univ}.
  If \code{univ} is \code{NULL} (default), then the gene IDs are those
  found in the annotation data source.}
  \item{chrGraph}{A \code{graph} object as returned by \code{makeChrBandGraph}}
}

\value{

  A (0/1) incidence matrix with gene IDs as rows and chromosome bands
  as columns.  A \code{1} in \code{m[i, j]} indicates that the geneID
  \code{rownames(m)[i]} is annotated at chromosome band
  \code{colnames(m)[j]}.

}

\author{Seth Falcon}

\seealso{\code{\link{makeChrBandGraph}}}
\examples{
\dontrun{
mm <- chrBandInciMat("hgu95av2")
}
}

\keyword{manip}

