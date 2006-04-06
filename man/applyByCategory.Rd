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
  \item{stats}{Numeric vector with test statistics of interest. If
    \code{Amat} is an object of class
    \code{\link[graph:graphNEL-class]{graphNEL}}, then \code{stats}
    needs to be named, with names corresponding to those in \code{Amat}.}
  \item{Amat}{A data structure that codes category memberships of the
    elements of \code{stats}. It can be either a matrix or a
    \code{\link[graph:graphNEL-class]{graphNEL}} object.
    See details.}
  \item{FUN}{A function to apply to the subsets \code{stats} by
    categories.}
  \item{\dots}{Extra parameters passed to \code{FUN}.}
}

\details{
  If \code{Amat} is an object of class
  \code{\link[graph:graphNEL-class]{graphNEL}}, the
  categories are taken
  to be those nodes of \code{Amat} whose names do not match the names of
  \code{stats}, and category membership is coded by edges from
  categories to non-category nodes.
  If it is a matrix, its columns correspond to values in
  \code{stats} and its rows to the different
  categories. If matrix element \code{[i,j]!=0}, then
  \code{stats[j]} belongs to the \code{i}-th categorization.

  For situations in which \code{stats} is large, the sparse
  graph representation may be more efficient.
}

\value{
  The return value is a list or vector of length equal to
  the number of categorizations. Each element corresponds to the
  values obtained by applying \code{FUN} to the subset of values
  in \code{stats} according to the category defined for that
  row.
}

\author{R. Gentleman, contributions from W. Huber}


\seealso{ \code{\link{split}}, \code{\link{apply}} }
\examples{
 set.seed(123)
 st = rnorm(20)

 aM = matrix(sample(c(0,1), 60, replace=TRUE), nr=3)
 rownames(aM) = paste("category", 1:3, sep="")
 colnames(aM) = names(st) = paste("g", 1:20, sep="")

 bpg = aM2bpG(aM)

 a1 = applyByCategory(st, aM, median)
 a2 = applyByCategory(st, bpg, median)

 stopifnot( identical(a2, a1) )
}

\keyword{manip}

