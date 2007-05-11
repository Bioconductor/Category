\name{chrBandInciMat}

\alias{chrBandInciMat}
\alias{makeChrBandInciMat}

\title{Create an incidence matrix mapping chromosome bands to genes}

\description{

  These functions return a 0/1 incidence matrix with a row for each
  chromosme band and a column for each gene.  Only those chromosome
  bands with at least one gene annotation will be included.
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

  A (0/1) incidence matrix with chromosome bands as rows and gene IDs
  as columns.  A \code{1} in \code{m[i, j]} indicates that the
  chromosome band \code{rownames(m)[i]} contains the geneID
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

