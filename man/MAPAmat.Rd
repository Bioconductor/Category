\name{MAPAmat}

\alias{MAPAmat}
\alias{chrBandInciMat}
\alias{makeChrBandInciMat}

\title{Create an incidence matrix mapping chromosome bands to genes}

\description{
  These functions return a 0/1 incidence matrix with a row for each
  chromosme band and a column for each gene.  Only those chromosome
  bands with at least one gene annotation will be included.
}

\usage{
MAPAmat(chip, univ = NULL, minCount = 0)
chrBandInciMat(chip, univ = NULL, minCount = 0)
makeChrBandInciMat(chrGraph)
}

\arguments{
  \item{chip}{A string giving the annotation source.  For example,
    \code{"hgu133plus2"}}

  \item{univ}{A vector of gene IDs (these should be Entrez IDs for
    most annotation sources).  The the annotations will be limited to
    those in the set specified by \code{univ}.  If \code{univ} is
    \code{NULL} (default), then the gene IDs are those found in the
    annotation data source.}

  \item{chrGraph}{A \code{graph} object as returned by
    \code{makeChrBandGraph}}

  \item{minCount}{ Bands with less than \code{minCount} genes will be
    excluded from the returned matrix.  If \code{minCount} is \code{0},
    no bands will be removed, this is the default.}
}

\details{
  \code{chrBandInciMat} is a DEPRECATED alias for \code{MAPAmat}.
}

\value{

  A (0/1) incidence matrix with chromosome bands as rows and gene IDs
  as columns.  A \code{1} in \code{m[i, j]} indicates that the
  chromosome band \code{rownames(m)[i]} contains the geneID
  \code{colnames(m)[j]}.

}

\author{Seth Falcon}

\seealso{
  \code{\link{makeChrBandGraph}},
  \code{\link{cateGOry}},
  \code{\link{probes2MAP}}
}

\examples{
have_hgu95av2.db <- suppressWarnings(require("hgu95av2.db"))
if (have_hgu95av2.db)
   mam <- MAPAmat("hgu95av2.db")
}

\keyword{manip}
