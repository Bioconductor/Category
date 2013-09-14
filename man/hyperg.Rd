\name{hyperg}
\alias{hyperg}
\alias{hyperg,character-method}
\alias{hyperg,list-method}

\title{

  Hypergeometric (gene set enrichment) tests on character vectors.

}

\description{

  This function performs a hypergeometric test for over- or
  under-representation of significant \sQuote{genes} amongst those
  assayed in a universe of genes. It provides an interface based on
  character vectors of identifying member of gene sets and the gene
  universe.

}

\usage{
hyperg(assayed, significant, universe,
    representation = c("over", "under"), ...)
}

\arguments{

  \item{assayed}{A vector of assayed genes (or other
    identifiers). \code{assayed} may be a character vector (defining a
    single gene set) or list of character vectors (defining a collection
    of gene sets).}

  \item{significant}{A vector of assayed genes that were
    differentially expressed. If \code{assayed} is a character vector,
    then \code{significant} must also be a character vector; likewise
    when \code{assayed} is a \code{list}.}

  \item{universe}{A character vector defining the universe of
    genes.}

  \item{representation}{Either \dQuote{over} or \dQuote{under}, to
    indicate testing for over- or under-representation, respectively, of
    differentially expressed genes.}

  \item{\dots}{Additional arguments, unused.}

}

\value{

  When invoked with a character vector of \code{assayed} genes, a named
  numeric vector providing the input values, P-value, odds ratio, and
  expected number of significantly expressed genes.

  When invoked with a list of character vectors of \code{assayed} genes,
  a data frame with columns of input values, P-value, odds ratio, and
  expected number of significantly expressed genes.

}

\author{
  Martin Morgan \url{mtmorgan@fhcrc.org} with contributions from Paul
  Shannon.
}

\seealso{

  \code{\link{hyperGTest}} for convenience functions using Bioconductor
  annotation resources such as GO.db.

}

\examples{
set.seed(123)

## artifical sets -- affy probes grouped by protein family
library(hgu95av2.db)
map <- select(hgu95av2.db, keys(hgu95av2.db), "PFAM")
sets <- Filter(function(x) length(x) >= 10, split(map$PROBEID, map$PFAM))

universe <- unlist(sets, use.names=FALSE)
siggenes <- sample(universe, length(universe) / 20)  ## simulate
sigsets <- Map(function(x, y) x[x \%in\% y], sets, MoreArgs=list(y=siggenes))

result <- hyperg(sets, sigsets, universe)
head(result)
}

\keyword{manip}
