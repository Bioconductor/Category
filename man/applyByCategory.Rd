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
  \item{stats}{Numeric vector with test statistics of interest. }
  \item{Amat}{A matrix whose columns correspond to values in
    \code{stats}, whose rows correspond to the different
    categorizations. Matrix element \code{[i,j]} codes the
    category membership of test statistic \code{stats[j]} in the
    \code{i}-th categorization.}
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
 element of the list corresponds to the values obtained by \code{split}ing
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

