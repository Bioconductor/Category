\name{geneGoHyperGeoTest}
\alias{geneGoHyperGeoTest}

\title{Hypergeometric Tests for GO }

\description{
  Given a set of unique Entrez Gene Identifiers, a microarray annotation
  data package name, and the GO ontology of interest, this function will
  compute Hypergeomtric p-values for overrepresentation of each GO term
  in the specified ontology among the GO annotations for the interesting
  genes (as indicated by the Entrez Gene ids).
}

\usage{
geneGoHyperGeoTest(entrezGeneIds, lib, ontology, universe=NULL)
}

\arguments{
  \item{entrezGeneIds}{A vector of Entrez Gene Identifiers}
  \item{lib}{A string giving the name of the annotation data package to
    use.  This must correspond to the microarray chip type that the data
    came from}
  \item{ontology}{One of "BP", "CC", or "MF" used to determine which GO
    ontology to use.}
  \item{universe}{A character vector of unique Entrez Gene identifiers.
    This is the population (the urn) of the Hypergeometric test.  When
    \code{NULL} (default), the population is all Entrez Gene ids in the
    annotation package that have a GO term annotation in the specified
    GO category (see details).}
}

\details{
  The Entrez Gene ids given in \code{entrezGeneIds} define the selected
  set of genes.  The universe of Entrez Gene ids is determined by the
  chip annotation data package (\code{lib}) or specified by the
  \code{universe} argument which must be a subset of the Entrez Gene ids
  represented on the chip.  Both the selected genes and the universe are
  reduced by removing Entrez Gene ids that do not have any annotations
  in the specified GO category.

  For each GO term in the specified category that has at least one
  annotation in the selected gene set (\code{entrezGeneIds}), we
  determine how many of its Entrez Gene annotations are in the universe
  set and how many are in the selected set.  With these counts we
  perform a Hypergeometric test using \code{phyper}.  This is equivalent
  to using Fisher's exact test.

  It is important that the correct chip annotation data package be
  identified as it determines the GO term to Entrez Gene id mapping as
  well as the universe of Entrez Gene ids in the case that the
  'universe' argument is omitted.

  For S. cerevisiae if the 'lib' argument is set to \code{"YEAST"} then
  comparisons and statistics are computed using common names and are
  with respect to all genes annotated in the S. cerevisiae genome not
  with respect to any microarray chip.  This will \bold{not} be the
  right thing to do if you are working with a yeast microarray.
}

\value{
  A \code{\link{GeneCategoryHyperGeoTestResult-class}} instance.
}

\author{S. Falcon}

\seealso{
  \code{\link{GeneCategoryHyperGeoTestResult-class}}
  \code{\link{GeneCategoryHyperGeoTestParams-class}}
  \code{\link{geneKeggHyperGeoTest}}
  \code{\link{geneCategoryHyperGeoTest}}
}

\examples{
library("hgu95av2")
library("GO")
set.seed(123)
probes <- ls(hgu95av2LOCUSID)
probes <- sample(probes, 100)
egIds <- unique(unlist(mget(probes, hgu95av2LOCUSID)))
ans <- geneGoHyperGeoTest(egIds, "hgu95av2", "BP")
print(ans)
}
\keyword{htest}

