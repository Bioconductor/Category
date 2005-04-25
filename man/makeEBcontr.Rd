\name{makeEBcontr}
\alias{makeEBcontr}
\title{ A function to make the contrast vectors needed for EBarrays }
\description{
 Using EBarrays to detect differential expression requires the
 construction of a set of contrasts. This little helper function
 computes these contrasts for a two level factor.
}
\usage{
makeEBcontr(f1, hival)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f1}{The factor that will define the contrasts.}
  \item{hival}{The \code{level} of the factor to treat as the high level.}
}
\details{
  Not much more to add, see EBarrays for more details. This is used in
  the Category package to let users compute the posterior probability of
  differential expression, and hence to compute expected numbers of
  differentially expressed genes, per category.
}
\value{
An object of class  ``ebarraysPatterns''.
}

\author{R. Gentleman }

\seealso{ \code{\link{ebPatterns}}}
\examples{
 if( require("EBarrays") ) {
  myfac = factor(rep(c("A", "B"), c(12, 24)))
  makeEBcontr(myfac, "B")
 }
}
\keyword{manip}

