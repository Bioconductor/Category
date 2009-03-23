\name{gseattperm}
\alias{gseattperm}

\title{Permutation p-values for GSEA}
\description{
  This function performs GSEA computations and returns p-values
  for each gene set based on repeated permutation of the phenotype
  labels.
}
\usage{
gseattperm(eset, fac, mat, nperm)
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

\examples{
## This example uses a random sample of probesets and a randomly
## generated category matrix.  The results, therefore, are not
## meaningful, but the code demonstrates how to use gseattperm without
## requiring any expensive computations.

## Obtain an ExpressionSet with two types of samples (mol.biol)
haveALL <- require("ALL")
if (haveALL) {
data(ALL)
set.seed(0xabcd)
rndIdx <- sample(1:nrow(ALL), 500)
Bcell <- grep("^B", as.character(ALL$BT))
typeNames <- c("NEG", "BCR/ABL")
bcrAblOrNegIdx <- which(as.character(ALL$mol.biol) \%in\% typeNames)
s <- ALL[rndIdx, intersect(Bcell, bcrAblOrNegIdx)]
s$mol.biol <- factor(s$mol.biol)

## Generate a random category matrix
nCats <- 100
set.seed(0xdcba)
rndCatMat <- matrix(sample(c(0L, 1L), replace=TRUE),
                    nrow=nCats, ncol=nrow(s),
                    dimnames=list(
                      paste("c", 1:nCats, sep=""),
                      featureNames(s)))

## Demonstrate use of gseattperm
N <- 10
pvals <- gseattperm(s, s$mol.biol, rndCatMat, N)
pvals[1:5, ]
}
}


