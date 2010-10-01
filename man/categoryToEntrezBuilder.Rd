\name{categoryToEntrezBuilder}
\docType{methods}

\alias{categoryToEntrezBuilder}
\alias{categoryToEntrezBuilder,GOHyperGParams-method}
\alias{categoryToEntrezBuilder,KEGGHyperGParams-method}
\alias{categoryToEntrezBuilder,PFAMHyperGParams-method}


\title{Return a list mapping category ids to Entrez Gene ids}
\description{
  Return a list mapping category ids to the Entrez Gene ids annotated at
  the category id.  Only those category ids that have at least one
  annotation in the set of Entrez Gene ids specified by the
  \code{geneIds} slot of \code{p} are included.
}
\usage{
categoryToEntrezBuilder(p)
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
  A list mapping category ids to Entrez Gene identifiers.
}

\author{S. Falcon}
\seealso{
  \code{\link{hyperGTest}}
  \code{\link{HyperGParams-class}}
}

\keyword{manip}
