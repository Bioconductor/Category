\name{findAMstats}
\alias{findAMstats}
\title{Compute per category summary statistics }
\description{
 For a given incidence matrix, \code{Amat}, compute some per category
 statistics. 
}
\usage{
findAMstats(Amat, tstats)
}
\arguments{
  \item{Amat}{An incidence matrix, with categories as the rows and
    probes as the columns. }
  \item{tstats}{A vector of per probe test statistics (should be the
    same length as \code{ncol(Amat)}. }
}
\details{
  Simple summary statistics are computed, such as the row sums and the
  vector of per category sums of the test statistics, \code{tstats}.
}
\value{
A list with components,
  \item{eDE }{per category sums of the test statistics}
  \item{lens }{row sums of \code{Amat}}
}

\author{R. Gentleman }

\seealso{\code{\link{applyByCategory}}}
\examples{

 ts = rnorm(100)
 Am = matrix(sample(c(0,1), 1000, replace=TRUE), ncol=100)
 findAMstats(Am, ts)
}

\keyword{manip}

