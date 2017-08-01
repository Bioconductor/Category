\name{universeBuilder}
\docType{methods}
\alias{universeBuilder}
\alias{universeBuilder,GOHyperGParams-method}
\alias{universeBuilder,KEGGHyperGParams-method}
\alias{universeBuilder,PFAMHyperGParams-method}
\alias{universeBuilder,OBOHyperGParams-method}

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
  \item{p}{A subclass of \code{HyperGParams-class}}
}
\details{
  End users \bold{should not} call this directly.  This method gets
  called from \code{hyperGTest}.  To add support for a new
  category, a new method for this generic must be defined.  Its
  signature should match a subclass of
  \code{HyperGParams-class} appropriate for the new
  category.
}
\value{
  A vector of gene identifiers.
}

\author{S. Falcon}
\seealso{
  \code{\link{hyperGTest}}
  \code{\link{HyperGParams-class}}
}

\keyword{manip}
