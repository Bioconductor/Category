\name{applyByCategory}
\alias{applyByCategory}
\title{Apply a function to a vector of statistics, by category }
\description{
  The rows of the input matrix, \code{Amat}, define categories. For each
  row we split the values of \code{stats} into sets according to
  the categories and the apply \code{FUN} to each of the sets.
}
\usage{
applyByCategory(stats, Amat, FUN = mean, ...)
}
\arguments{
  \item{stats}{Test statistics of interest. }
  \item{Amat}{A matrix with columns corresponding to values in
    \code{stats} and rows defining categorizations. }
  \item{FUN}{A function to apply to the splits of \code{stats} by
    categories.}
  \item{\dots}{Extra parameters passed to \code{lapply} and used by
    \code{FUN}.}
}
\details{
  Pretty much all in the description field.
}
\value{
 The return value is a list of length equal to \code{nrow(Amat)}. Each
 element of which corresponds to the values obtained by \code{split}ing
 the vector \code{stats} according to the categories defined for that
 row and then \code{apply}ing \code{FUN} to each resultant set.
}

\author{R. Gentleman}


\seealso{ \code{\link{split}}, \code{\link{apply}} }
\examples{
 st1 = rnorm(20)
 catM = matrix(sample(c(0,1), 60, replace=TRUE), nr=3)
 applyByCategory(st1, catM, median)

}

\keyword{manip}

