\name{gseaperm}
\alias{gseaperm}

\title{Permutation p-values for GSEA}
\description{
  This function performs a GSEA computations and returns p-values
  for each gene set based repeated permutation of the phenotype
  labels.
}
\usage{
gseaperm(eset, fac, mat, nperm)
}

\arguments{
  \item{eset}{An \code{ExpressionSet} object}
  \item{fac}{A \code{factor} identifying the phenotypes in
    \code{eset}.  Usually, this will be one of the columns in the
    phenotype data associated with \code{eset}.}
  \item{mat}{A 0/1 incidence matrix with each row representing a gene
    set and each column representing a gene.  A 1 indicates membership
    of a gene in a gene set.}
  \item{nperm}{Number of permutations to test to build the reference
    distribution.}
}
\details{
  The t-statistic is used (via \code{rowttests}) to test for a
  difference in means between the phenotypes determined by \code{fac}
  within each gene set (given as a row of \code{mat}).

  A reference distribution for these statistics is established by
  permuting \code{fac} and repeating the test \code{B} times.
}
\value{
  A matrix with the same number of rows as \code{mat} and two columns,
  \code{"Lower"} and \code{"Upper"}.  The \code{"Lower"}
  (\code{"Upper"}) column gives the probability of seeing a t-statistic
  smaller (larger) than the observed.
}

\author{Seth Falcon}

\keyword{ }

