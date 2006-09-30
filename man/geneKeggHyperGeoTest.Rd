\name{geneKeggHyperGeoTest}
\alias{geneKeggHyperGeoTest}

\title{Hypergeometric Tests for KEGG (DEFUNCT)}

\description{
        This function is defunct.  Use \code{hyperGTest} instead.
}
\usage{
geneKeggHyperGeoTest(entrezGeneIds, lib, universe=NULL)
}

\arguments{
    \item{entrezGeneIds}{A vector of Entrez Gene Identifiers}
    \item{lib}{A string giving the name of the annotation data package to
      use.  This must correspond to the microarray chip type that the data
      came from}
  \item{universe}{A character vector of unique Entrez Gene identifiers.
    This is the population (the urn) of the Hypergeometric test.  When
    \code{NULL} (default), the population is all Entrez Gene ids in the
    annotation package that have a KEGG pathway annotation (see
    details).}
}
\details{
  The Entrez Gene ids given in \code{entrezGeneIds} define the selected
  set of genes.  The universe of Entrez Gene ids is determined by the
  chip annotation data package (\code{lib}) or specified by the
  \code{universe} argument which must be a subset of the Entrez Gene ids
  represented on the chip.  Both the selected genes and the universe are
  reduced by removing Entrez Gene ids that do not have any KEGG pathway
  annotations.

  For each KEGG pathway that has at least one annotation in the selected
  gene set (\code{entrezGeneIds}), we determine how many of its Entrez
  Gene annotations are in the universe set and how many are in the
  selected set.  With these counts we perform a Hypergeometric test
  using \code{phyper}.  This is equivalent to using Fisher's exact test.

  It is important that the correct chip annotation data package be
  identified as it determines the KEGG pathway to Entrez Gene id mapping
  as well as the universe of Entrez Gene ids in the case that the
  'universe' argument is omitted.

  For S. cerevisiae if the 'lib' argument is set to \code{"YEAST"} then
  comparisons and statistics are computed using common names and are
  with respect to all genes annotated in the S. cerevisiae genome not
  with respect to any microarray chip.  This will \bold{not} be the
  right thing to do if you are working with a yeast microarray.
}
\value{
  A \code{\link{HyperGResult-class}} instance.
}

\author{S. Falcon}

\seealso{
  \code{\link{HyperGResult-class}}
  \code{\link{HyperGParams-class}}
  \code{\link{geneGoHyperGeoTest}}
  \code{\link{hyperGTest}}
}

\keyword{htest}

