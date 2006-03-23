\name{universeBuilder}
\docType{genericFunction}
\alias{universeBuilder}
\alias{universeBuilder,GeneGoHyperGeoTestParams-method}
\alias{universeBuilder,GeneKeggHyperGeoTestParams-method}

\title{Return a vector of gene identifiers with category annotations}
\description{
  Return all gene ids that are annotated at one or more terms in the
  category.  If the \code{universeGeneIds} slot of \code{p} has length
  greater than zero, then the intersection of the gene ids specified in
  that slot and the normal return value is given.
}
\usage{
universeBuilder(p)
}

\arguments{
  \item{p}{A subclass of \code{GeneCategoryHyperGeoTestParams-class}}
}
\details{
  End users \bold{should not} call this directly.  This method gets
  called from \code{geneCategoryHyperGeoTest}.  To add support for a new
  category, a new method for this generic must be defined.  Its
  signature should match a subclass of
  \code{GeneCategoryHyperGeoTestParams-class} appropriate for the new
  category.
}
\value{
  A vector of gene identifiers.
}

\author{S. Falcon}
\seealso{
  \code{\link{geneCategoryHyperGeoTest}}
  \code{\link{GeneCategoryHyperGeoTestParams-class}}
}

\keyword{manip}
