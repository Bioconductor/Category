\name{NewChrBandTree}
\alias{NewChrBandTree}
\alias{ChrBandTreeFromGraph}

\title{Create a new ChrBandTree object}
\description{
  \code{NewChrBandTree} and \code{ChrBandTreeFromGraph} provide
  constructors for the \code{ChrBandTree} class.
}
\usage{
NewChrBandTree(chip, univ)
ChrBandTreeFromGraph(g)
}

\arguments{
  \item{chip}{The name of an annotation data package}
  \item{univ}{A vector of gene identifiers that defines the universe of
    genes.  Usually, this will be a vector of Entez Gene IDs.  If
    \code{univ} is \code{NULL}, then all genes probed on the specified
    chip will be in the universe.  We strongly recommend using the set
    of genes that remains after applying a non-specific filter as the
    universe.}
  \item{g}{A \code{graph} instance as returned by \code{makeChrBandGraph}}
}
\value{
  A new \code{ChrBandTree} instance.
}

\author{S. Falcon}
\seealso{\code{\link{ChrBandTree-class}}}


