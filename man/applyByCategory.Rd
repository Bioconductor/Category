\name{applyByCategory}
\alias{applyByCategory}
\title{Apply a function to a vector of statistics, by category }
\description{
  For each category, apply the function \code{FUN} to the set of values
  of \code{stats} belonging to that category.
}

\usage{
applyByCategory(stats, Amat, FUN = mean, ...)
}

\arguments{
  \item{stats}{Numeric vector with test statistics of interest.}
  \item{Amat}{A logical matrix: the adjacency matrix of the
    bipartite genes - category graph.
    Its rows correspond to the categories, columns
    to the genes, and \code{TRUE} indicates membership.
    The columns are assumed to be aligned with the elements of \code{stats}.}
  \item{FUN}{A function to apply to the subsets \code{stats} by
    categories.}
  \item{\dots}{Extra parameters passed to \code{FUN}.}
}

\details{
  For GO categories, the function \code{\link{cateGOry}} might be useful
  for the construction of \code{Amat}.
}

\value{
  The return value is a list or vector of length equal to
  the number of categories. Each element corresponds to the
  values obtained by applying \code{FUN} to the subset of values
  in \code{stats} according to the category defined for that
  row.
}

\author{R. Gentleman, contributions from W. Huber}


\seealso{\code{\link{apply}} }
\examples{
 set.seed(0xabcd)
 st = rnorm(20)
 names(st) = paste("gene", 1:20)

 a = matrix(sample(c(FALSE, TRUE), 60, replace=TRUE), nrow=3,
        dimnames = list(paste("category", LETTERS[1:3]), names(st)))

 applyByCategory(st, a, median)
}

\keyword{manip}

