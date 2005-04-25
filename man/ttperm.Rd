\name{ttperm}
\alias{ttperm}
\title{A simple function to compute a permutation t-test. }
\description{
  The data matrix, \code{x}, with two-level factor, \code{fac}, is used
  to compute t-tests. The values of \code{fac} are permuted \code{B}
  times and the complete set of t-tests is performed for each
  permutation. 
}
\usage{
ttperm(x, fac, B = 100, tsO = TRUE)
}
\arguments{
  \item{x}{ A data matrix. The number of columns should be the same as
    the length of \code{fac}. }
  \item{fac}{A factor with two levels. }
  \item{B}{An integer specifying the number of permutations. }
  \item{tsO}{A logical indicating whether to compute only the t-test
    statistic for each permuation. If \code{FALSE} then p-values are
    also computed - but this can be very slow. }
}
\details{
 Not much more to say. Probably there is a generic function somewhere,
 but I could not find it.
}
\value{
 A list, the first element is named \code{obs} and contains the true,
 observed, values of the t-statistic. The second element is named
 \code{ans}  and contains a list of length \code{B} containing the
 different permuations. 
}

\author{R. Gentleman }

\seealso{\code{\link{rowttests}} }
\examples{
   x=matrix(rnorm(100), nc=10)
   y = factor(rep(c("A","B"), c(5,5)))
   ttperm(x, y, 10)

}
\keyword{manip}

