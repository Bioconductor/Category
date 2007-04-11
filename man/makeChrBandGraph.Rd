\name{makeChrBandGraph}
\alias{makeChrBandGraph}

\title{Create a graph representing chromosome band annotation data}
\description{
  This function returns a \code{graph} object representing the nested
  structure of chromosome bands (also known as cytogenetic bands).
  The nodes of the graph are band identifiers.  Each node has a
  \code{geneIds} node attribute that lists the gene IDs that are
  annotated at the band (the gene IDs will be Entrez IDs in most
  cases).
}
\usage{
makeChrBandGraph(chip, univ = NULL)
}

\arguments{
  \item{chip}{A string giving the annotation source.  For example, \code{"hgu133plus2"}}
  \item{univ}{A vector of gene IDs (these should be Entrez IDs for
  most annotation sources).  The annotations
  attached to the graph will be limited to those specified by \code{univ}.
  If \code{univ} is \code{NULL} (default), then the gene IDs are those
  found in the annotation data source.}
}
\details{
  This function parses the data stored in the
  \code{<chip>MAP} map from the appropriate annotation data package.
  Although cytogenetic bands are observed in all organisms, currently,
  only human band nomenclature is supported.
}
\value{
  A \code{\link[graph]{graph-class}} instance.  The graph will be a
  tree and the root node is labeled for the organism.
}

\author{Seth Falcon}

\examples{
chrGraph <- makeChrBandGraph("hgu95av2")
chrGraph
}

\keyword{manip}

